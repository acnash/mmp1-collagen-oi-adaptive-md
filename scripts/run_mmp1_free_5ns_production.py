#!/usr/bin/env python3
"""Continue all apo mutant-collagen replicas for 5 ns of unrestrained MD."""

from __future__ import annotations

import argparse
import ctypes
import json
import math
import os
import time
import traceback
from datetime import datetime
from pathlib import Path

from openmm import LangevinMiddleIntegrator, MonteCarloBarostat, Platform, XmlSerializer, unit
from openmm.app import (
    CheckpointReporter,
    GromacsGroFile,
    GromacsTopFile,
    HBonds,
    PME,
    PDBFile,
    Simulation,
    StateDataReporter,
    XTCReporter,
)


VARIANTS = ["G978S", "G984C", "G987R"]


def set_windows_sleep_prevention(enabled):
    """Keep Windows awake while an unattended production campaign is active."""
    if os.name != "nt":
        return
    es_continuous = 0x80000000
    es_system_required = 0x00000001
    flags = es_continuous | es_system_required if enabled else es_continuous
    result = ctypes.windll.kernel32.SetThreadExecutionState(flags)
    if result == 0:
        raise OSError("SetThreadExecutionState failed")


def enable_periodic_bonded_forces(system):
    classes = {
        "HarmonicBondForce", "HarmonicAngleForce", "PeriodicTorsionForce",
        "RBTorsionForce", "CMAPTorsionForce", "CustomBondForce",
        "CustomAngleForce", "CustomTorsionForce",
    }
    for force in system.getForces():
        if force.__class__.__name__ in classes and hasattr(force, "setUsesPeriodicBoundaryConditions"):
            force.setUsesPeriodicBoundaryConditions(True)


def finite_state(state):
    energy = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    if not math.isfinite(energy):
        return False
    return all(
        math.isfinite(v)
        for p in state.getPositions().value_in_unit(unit.nanometer)
        for v in (p.x, p.y, p.z)
    )


def write_json(path, value):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2) + "\n")
    temporary.replace(path)


def run_replica(repo, variant, replica, args):
    equil_summary_path = repo / "mmp1_free_conventional_md" / variant / f"replica_{replica}" / "summary.json"
    equil_summary = json.loads(equil_summary_path.read_text())
    equil_dir = Path(equil_summary["output_directory"])
    initial_state_path = equil_dir / "06_npt_unrestrained_state.xml"
    initial_state = XmlSerializer.deserialize(initial_state_path.read_text())

    system_dir = repo / "generated_mutants" / "mmp1_free_150mM_NaCl" / f"collagen_{variant}"
    gro = GromacsGroFile(str(system_dir / f"MMP1_free_{variant}_150mM_NaCl_start.gro"))
    top = GromacsTopFile(
        str(system_dir / "system.top"), periodicBoxVectors=gro.getPeriodicBoxVectors(), includeDir=str(system_dir)
    )
    system = top.createSystem(
        nonbondedMethod=PME,
        nonbondedCutoff=1.0 * unit.nanometer,
        constraints=HBonds,
        rigidWater=True,
        removeCMMotion=True,
    )
    enable_periodic_bonded_forces(system)
    system.addForce(MonteCarloBarostat(args.pressure_bar * unit.bar, args.temperature_k * unit.kelvin, 25))

    seed = args.master_seed + 1000 * VARIANTS.index(variant) + 10 * replica
    integrator = LangevinMiddleIntegrator(
        args.temperature_k * unit.kelvin,
        1.0 / unit.picosecond,
        args.timestep_fs * 0.001 * unit.picosecond,
    )
    integrator.setRandomNumberSeed(seed)
    platform = Platform.getPlatformByName(args.platform)
    properties = {"DeviceIndex": str(args.device_index), "Precision": "mixed", "UseCpuPme": "false"}
    simulation = Simulation(top.topology, system, integrator, platform, properties)

    output_dir = repo / "mmp1_free_5ns_production" / variant / f"replica_{replica}"
    output_dir.mkdir(parents=True, exist_ok=True)
    checkpoint_path = output_dir / "production.chk"
    progress_path = output_dir / "progress.json"
    summary_path = output_dir / "summary.json"
    if summary_path.exists() and not args.force:
        previous = json.loads(summary_path.read_text())
        if previous.get("status") == "completed":
            print(f"Skipping completed {variant} replica {replica}", flush=True)
            return previous

    resumed = checkpoint_path.exists() and not args.force
    if resumed:
        simulation.loadCheckpoint(str(checkpoint_path))
    else:
        simulation.context.setPositions(initial_state.getPositions())
        simulation.context.setVelocities(initial_state.getVelocities())
        simulation.context.setPeriodicBoxVectors(*initial_state.getPeriodicBoxVectors())
        simulation.context.setTime(initial_state.getTime())

    timestep_ps = args.timestep_fs * 0.001
    target_steps = int(round(args.duration_ns * 1000.0 / timestep_ps))
    trajectory_steps = int(round(args.trajectory_interval_ps / timestep_ps))
    log_steps = int(round(args.log_interval_ps / timestep_ps))
    checkpoint_steps = int(round(args.checkpoint_interval_ps / timestep_ps))
    chunk_steps = int(round(args.progress_interval_ps / timestep_ps))
    if simulation.currentStep > target_steps:
        raise RuntimeError("Checkpoint step exceeds requested production length")

    append = resumed and (output_dir / "production.xtc").exists()
    simulation.reporters.append(
        XTCReporter(str(output_dir / "production.xtc"), trajectory_steps, append=append, enforcePeriodicBox=True)
    )
    simulation.reporters.append(
        StateDataReporter(
            str(output_dir / "production.csv"), log_steps, append=append,
            step=True, time=True, potentialEnergy=True, kineticEnergy=True,
            totalEnergy=True, temperature=True, volume=True, density=True,
            speed=True, remainingTime=True, totalSteps=target_steps, separator=",",
        )
    )
    simulation.reporters.append(CheckpointReporter(str(checkpoint_path), checkpoint_steps))

    started = time.time()
    initial_step = simulation.currentStep
    print(
        f"Starting {variant} replica {replica}: step {initial_step}/{target_steps}"
        f" ({args.duration_ns:.3f} ns target)", flush=True
    )
    while simulation.currentStep < target_steps:
        simulation.step(min(chunk_steps, target_steps - simulation.currentStep))
        state = simulation.context.getState(getPositions=True, getEnergy=True)
        if not finite_state(state):
            raise RuntimeError("Non-finite energy or coordinates during production")
        progress = {
            "variant": variant,
            "replica": replica,
            "status": "running",
            "current_step": simulation.currentStep,
            "target_steps": target_steps,
            "completed_ns": simulation.currentStep * timestep_ps / 1000.0,
            "target_ns": args.duration_ns,
            "resumed_from_checkpoint": resumed,
            "updated": datetime.now().isoformat(),
        }
        write_json(progress_path, progress)
        print(
            f"PROGRESS {variant} replica {replica}: {progress['completed_ns']:.2f}/{args.duration_ns:.2f} ns",
            flush=True,
        )

    final_state = simulation.context.getState(
        getPositions=True, getVelocities=True, getEnergy=True, enforcePeriodicBox=True
    )
    if not finite_state(final_state):
        raise RuntimeError("Non-finite final production state")
    (output_dir / "final_state.xml").write_text(XmlSerializer.serialize(final_state) + "\n")
    checkpoint_path.write_bytes(simulation.context.createCheckpoint())
    with (output_dir / "final_structure.pdb").open("w") as handle:
        PDBFile.writeFile(top.topology, final_state.getPositions(), handle, keepIds=True)
    summary = {
        "variant": variant,
        "replica": replica,
        "status": "completed",
        "duration_ns": args.duration_ns,
        "timestep_fs": args.timestep_fs,
        "temperature_k": args.temperature_k,
        "pressure_bar": args.pressure_bar,
        "platform": args.platform,
        "device_index": args.device_index,
        "integrator_seed": seed,
        "initial_state": str(initial_state_path),
        "coordinates_and_velocities_continued": True,
        "positional_restraints_present": False,
        "adaptive_sampling": False,
        "resumed_from_checkpoint": resumed,
        "initial_step": initial_step,
        "final_step": simulation.currentStep,
        "potential_energy_kj_mol": final_state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole),
        "wall_time_seconds_this_invocation": time.time() - started,
        "completed": datetime.now().isoformat(),
    }
    write_json(summary_path, summary)
    write_json(progress_path, {**summary, "target_steps": target_steps})
    print(f"Completed {variant} replica {replica}", flush=True)
    return summary


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--variants", nargs="+", choices=VARIANTS, default=VARIANTS)
    parser.add_argument("--replicas", nargs="+", type=int, default=[1, 2, 3])
    parser.add_argument("--duration-ns", type=float, default=5.0)
    parser.add_argument("--timestep-fs", type=float, default=2.0)
    parser.add_argument("--temperature-k", type=float, default=310.15)
    parser.add_argument("--pressure-bar", type=float, default=1.0)
    parser.add_argument("--platform", choices=["OpenCL", "CUDA"], default="OpenCL")
    parser.add_argument("--device-index", type=int, default=0)
    parser.add_argument("--master-seed", type=int, default=202611150)
    parser.add_argument("--trajectory-interval-ps", type=float, default=10.0)
    parser.add_argument("--log-interval-ps", type=float, default=10.0)
    parser.add_argument("--checkpoint-interval-ps", type=float, default=100.0)
    parser.add_argument("--progress-interval-ps", type=float, default=100.0)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    set_windows_sleep_prevention(True)
    try:
        results = []
        for variant in args.variants:
            for replica in args.replicas:
                try:
                    results.append(run_replica(args.repo_root, variant, replica, args))
                except Exception as exc:
                    failure = {
                        "variant": variant, "replica": replica, "status": "failed",
                        "error": repr(exc), "traceback": traceback.format_exc(),
                        "failed": datetime.now().isoformat(),
                    }
                    failure_dir = args.repo_root / "mmp1_free_5ns_production" / variant / f"replica_{replica}"
                    failure_dir.mkdir(parents=True, exist_ok=True)
                    write_json(failure_dir / "failure.json", failure)
                    print(json.dumps(failure, indent=2), flush=True)
                    raise
        print(json.dumps(results, indent=2))
    finally:
        set_windows_sleep_prevention(False)


if __name__ == "__main__":
    main()
