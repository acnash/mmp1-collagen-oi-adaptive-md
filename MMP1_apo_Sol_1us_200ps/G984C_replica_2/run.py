#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import math
import os
import signal
import subprocess
import time
from pathlib import Path

import numpy as np
import openmm
from openmm import XmlSerializer, unit
from openmm.app import GromacsGroFile, GromacsTopFile, HBonds, PME, Simulation


JOB_ROOT = Path("/scratch/anash19/MMP1_apo_Sol_1us_200ps/G984C_replica_2")
INPUT_DIR = JOB_ROOT / "input"
OUTPUT_DIR = JOB_ROOT / "output"
VARIANT = "G984C"
REPLICA = 2
SEED = 11520002

TEMPERATURE_K = 310.15
PRESSURE_BAR = 1.0
TIMESTEP_PS = 0.002
LOCAL_PRODUCTION_NS = 5.0
TOTAL_TARGET_NS = 1000.0
REMAINING_NS = TOTAL_TARGET_NS - LOCAL_PRODUCTION_NS
TARGET_STEPS = int(REMAINING_NS * 1000.0 / TIMESTEP_PS)
OUTPUT_INTERVAL_STEPS = int(200.0 / TIMESTEP_PS)
MAX_ALLOCATION_SECONDS = 6 * 24 * 3600 + 23 * 3600

stop_requested = False


def request_stop(signum, _frame):
    global stop_requested
    stop_requested = True
    print(f"Received signal {signum}; finishing the current 200 ps block.", flush=True)


def atomic_bytes(path: Path, data: bytes):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_bytes(data)
    os.replace(temporary, path)


def atomic_text(path: Path, text: str):
    atomic_bytes(path, text.encode("utf-8"))


def enable_periodic_bonded_forces(system):
    names = {
        "HarmonicBondForce", "HarmonicAngleForce", "PeriodicTorsionForce",
        "RBTorsionForce", "CMAPTorsionForce", "CustomBondForce",
        "CustomAngleForce", "CustomTorsionForce",
    }
    for force in system.getForces():
        if force.__class__.__name__ in names and hasattr(force, "setUsesPeriodicBoundaryConditions"):
            force.setUsesPeriodicBoundaryConditions(True)


def calculate_temperature(state, system):
    degrees = 3 * system.getNumParticles() - system.getNumConstraints()
    if any(isinstance(system.getForce(i), openmm.CMMotionRemover) for i in range(system.getNumForces())):
        degrees -= 3
    return (2 * state.getKineticEnergy() / (degrees * unit.MOLAR_GAS_CONSTANT_R)).value_in_unit(unit.kelvin)


def snapshot_path(step):
    return OUTPUT_DIR / "snapshots_200ps" / f"snapshot_step_{step:09d}.npz"


def write_snapshot(state, system, step):
    path = snapshot_path(step)
    temporary = path.with_suffix(".tmp")
    positions = np.asarray(state.getPositions(asNumpy=True).value_in_unit(unit.nanometer), dtype=np.float32)
    velocities = np.asarray(
        state.getVelocities(asNumpy=True).value_in_unit(unit.nanometer / unit.picosecond),
        dtype=np.float32,
    )
    box = np.asarray(state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.nanometer), dtype=np.float64)
    with temporary.open("wb") as handle:
        np.savez_compressed(
            handle,
            production_step=np.int64(step),
            added_production_ns=np.float64(step * TIMESTEP_PS / 1000.0),
            total_unrestrained_ns=np.float64(LOCAL_PRODUCTION_NS + step * TIMESTEP_PS / 1000.0),
            state_time_ps=np.float64(state.getTime().value_in_unit(unit.picosecond)),
            positions_nm=positions,
            velocities_nm_per_ps=velocities,
            box_vectors_nm=box,
            potential_energy_kj_mol=np.float64(state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)),
            kinetic_energy_kj_mol=np.float64(state.getKineticEnergy().value_in_unit(unit.kilojoule_per_mole)),
            temperature_K=np.float64(calculate_temperature(state, system)),
        )
    os.replace(temporary, path)


def reconcile_thermodynamics(path, current_step):
    header = [
        "production_step", "added_production_ns", "total_unrestrained_ns", "state_time_ps",
        "temperature_K", "potential_energy_kj_mol", "kinetic_energy_kj_mol", "volume_nm3",
    ]
    rows = []
    if path.exists():
        with path.open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle):
                if int(row["production_step"]) <= current_step:
                    rows.append(row)
    temporary = path.with_suffix(".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=header)
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, path)
    return {int(row["production_step"]) for row in rows}


def append_thermodynamics(path, state, system, step):
    with path.open("a", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow([
            step,
            step * TIMESTEP_PS / 1000.0,
            LOCAL_PRODUCTION_NS + step * TIMESTEP_PS / 1000.0,
            state.getTime().value_in_unit(unit.picosecond),
            calculate_temperature(state, system),
            state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole),
            state.getKineticEnergy().value_in_unit(unit.kilojoule_per_mole),
            state.getPeriodicBoxVolume().value_in_unit(unit.nanometer**3),
        ])


def write_progress(state, system, simulation, status):
    value = {
        "status": status,
        "variant": VARIANT,
        "replica": REPLICA,
        "production_step": simulation.currentStep,
        "added_production_ns": simulation.currentStep * TIMESTEP_PS / 1000.0,
        "total_unrestrained_ns": LOCAL_PRODUCTION_NS + simulation.currentStep * TIMESTEP_PS / 1000.0,
        "target_unrestrained_ns": TOTAL_TARGET_NS,
        "output_interval_ps": 200.0,
        "temperature_K": calculate_temperature(state, system),
        "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
        "updated_unix": time.time(),
    }
    atomic_text(OUTPUT_DIR / "progress.json", json.dumps(value, indent=2) + "\n")


def request_requeue():
    job_id = os.environ.get("SLURM_JOB_ID")
    if not job_id:
        raise RuntimeError("A restart was requested but SLURM_JOB_ID is unavailable")
    print(f"Checkpoint complete. Requeueing Slurm job {job_id}.", flush=True)
    subprocess.run(["scontrol", "requeue", job_id], check=True)


def main():
    signal.signal(signal.SIGTERM, request_stop)
    if hasattr(signal, "SIGUSR1"):
        signal.signal(signal.SIGUSR1, request_stop)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    (OUTPUT_DIR / "snapshots_200ps").mkdir(exist_ok=True)
    if (OUTPUT_DIR / "complete.json").exists():
        print(f"{VARIANT} replica {REPLICA} is already complete.", flush=True)
        return

    gro = GromacsGroFile(str(INPUT_DIR / "start.gro"))
    topology = GromacsTopFile(
        str(INPUT_DIR / "system.top"),
        periodicBoxVectors=gro.getPeriodicBoxVectors(),
        includeDir=str(INPUT_DIR),
    )
    system = topology.createSystem(
        nonbondedMethod=PME,
        nonbondedCutoff=1.0 * unit.nanometer,
        constraints=HBonds,
        rigidWater=True,
        removeCMMotion=True,
    )
    enable_periodic_bonded_forces(system)
    if any(isinstance(system.getForce(i), openmm.CustomExternalForce) for i in range(system.getNumForces())):
        raise RuntimeError("A positional-restraint force is present")
    barostat = openmm.MonteCarloBarostat(PRESSURE_BAR * unit.bar, TEMPERATURE_K * unit.kelvin, 25)
    barostat.setRandomNumberSeed(SEED + 1_000_000)
    system.addForce(barostat)

    integrator = openmm.LangevinMiddleIntegrator(
        TEMPERATURE_K * unit.kelvin,
        1.0 / unit.picosecond,
        TIMESTEP_PS * unit.picosecond,
    )
    integrator.setConstraintTolerance(1.0e-6)
    integrator.setRandomNumberSeed(SEED)
    platform = openmm.Platform.getPlatformByName("CUDA")
    simulation = Simulation(
        topology.topology,
        system,
        integrator,
        platform,
        {"Precision": "mixed", "DeviceIndex": "0"},
    )

    checkpoint = OUTPUT_DIR / "restart.chk"
    if checkpoint.exists():
        simulation.loadCheckpoint(str(checkpoint))
        continuation = "Sol binary checkpoint"
    else:
        initial = XmlSerializer.deserialize((INPUT_DIR / "initial_state.xml").read_text(encoding="utf-8"))
        simulation.context.setState(initial)
        simulation.context.setStepCount(0)
        continuation = "portable 5 ns endpoint state"

    if simulation.currentStep < 0 or simulation.currentStep > TARGET_STEPS:
        raise RuntimeError(f"Invalid production step {simulation.currentStep}")
    if simulation.currentStep % OUTPUT_INTERVAL_STEPS != 0:
        raise RuntimeError("Restart checkpoint is not aligned to a 200 ps boundary")

    thermo_path = OUTPUT_DIR / "thermodynamics_200ps.csv"
    recorded_steps = reconcile_thermodynamics(thermo_path, simulation.currentStep)
    if simulation.currentStep > 0:
        state = simulation.context.getState(
            getPositions=True, getVelocities=True, getEnergy=True, enforcePeriodicBox=True
        )
        if not snapshot_path(simulation.currentStep).exists():
            write_snapshot(state, system, simulation.currentStep)
        if simulation.currentStep not in recorded_steps:
            append_thermodynamics(thermo_path, state, system, simulation.currentStep)

    current = simulation.context.getState(getEnergy=True)
    print(
        f"Starting {VARIANT} replica {REPLICA}: continuation={continuation}, "
        f"total={LOCAL_PRODUCTION_NS + simulation.currentStep * TIMESTEP_PS / 1000.0:.1f}/1000.0 ns, "
        f"T={calculate_temperature(current, system):.2f} K, restraints=off",
        flush=True,
    )
    allocation_started = time.monotonic()

    while simulation.currentStep < TARGET_STEPS:
        steps = min(OUTPUT_INTERVAL_STEPS, TARGET_STEPS - simulation.currentStep)
        simulation.step(steps)
        state = simulation.context.getState(
            getPositions=True, getVelocities=True, getEnergy=True, enforcePeriodicBox=True
        )
        atomic_bytes(checkpoint, simulation.context.createCheckpoint())
        write_snapshot(state, system, simulation.currentStep)
        append_thermodynamics(thermo_path, state, system, simulation.currentStep)
        status = "complete" if simulation.currentStep == TARGET_STEPS else "running"
        write_progress(state, system, simulation, status)
        print(
            f"{VARIANT} replica {REPLICA}: "
            f"{LOCAL_PRODUCTION_NS + simulation.currentStep * TIMESTEP_PS / 1000.0:.1f}/1000.0 ns",
            flush=True,
        )

        wall_limit = time.monotonic() - allocation_started >= MAX_ALLOCATION_SECONDS
        if simulation.currentStep < TARGET_STEPS and (stop_requested or wall_limit):
            request_requeue()
            return

    final_state = simulation.context.getState(
        getPositions=True, getVelocities=True, getEnergy=True, enforcePeriodicBox=True
    )
    atomic_text(OUTPUT_DIR / "final_1us_state.xml", XmlSerializer.serialize(final_state) + "\n")
    completion = {
        "status": "complete",
        "variant": VARIANT,
        "replica": REPLICA,
        "local_starting_production_ns": LOCAL_PRODUCTION_NS,
        "sol_added_production_ns": REMAINING_NS,
        "total_unrestrained_ns": TOTAL_TARGET_NS,
        "timestep_fs": 2.0,
        "temperature_K": TEMPERATURE_K,
        "pressure_bar": PRESSURE_BAR,
        "position_restraints": False,
        "output_interval_ps": 200.0,
        "snapshot_count": TARGET_STEPS // OUTPUT_INTERVAL_STEPS,
        "coordinates_and_velocities_in_every_snapshot": True,
        "openmm_version": openmm.__version__,
        "completed_unix": time.time(),
    }
    atomic_text(OUTPUT_DIR / "complete.json", json.dumps(completion, indent=2) + "\n")
    print(f"Completed {VARIANT} replica {REPLICA}: 1.000 microsecond total.", flush=True)


if __name__ == "__main__":
    main()
