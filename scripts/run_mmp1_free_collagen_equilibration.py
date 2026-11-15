#!/usr/bin/env python3
"""Conventional equilibration for MMP-1-free periodic mutant collagen.

Protocol per independent replica:
  * restrained energy minimisation
  * 20 ps NVT heating from 10 K to 310.15 K, k=1000 kcal mol-1 A-2
  * 50 ps NPT, k=1000 kcal mol-1 A-2
  * 50 ps NPT, k=500 kcal mol-1 A-2
  * 50 ps NPT, k=100 kcal mol-1 A-2
  * 50 ps NPT, k=10 kcal mol-1 A-2
  * 100 ps NPT, unrestrained

All restraints act on collagen heavy atoms. Coordinates, velocities and box
vectors are transferred without reinitialisation between stages. Replicas use
different initial velocity and stochastic-integrator seeds. This is ordinary
MD equilibration, not adaptive sampling.
"""

from __future__ import annotations

import argparse
import json
import math
import shutil
import time
import traceback
from datetime import datetime
from pathlib import Path

from openmm import (
    CustomExternalForce,
    LangevinMiddleIntegrator,
    MonteCarloBarostat,
    Platform,
    XmlSerializer,
    unit,
)
from openmm.app import (
    GromacsGroFile,
    GromacsTopFile,
    HBonds,
    PME,
    PDBFile,
    Simulation,
    StateDataReporter,
    XTCReporter,
)


SYSTEMS = {
    "G978S": "collagen_G978S",
    "G984C": "collagen_G984C",
    "G987R": "collagen_G987R",
}

KCAL_PER_MOL_A2_TO_KJ_PER_MOL_NM2 = 418.4


def enable_periodic_bonded_forces(system):
    classes = {
        "HarmonicBondForce", "HarmonicAngleForce", "PeriodicTorsionForce",
        "RBTorsionForce", "CMAPTorsionForce", "CustomBondForce",
        "CustomAngleForce", "CustomTorsionForce",
    }
    for force in system.getForces():
        if force.__class__.__name__ in classes and hasattr(force, "setUsesPeriodicBoundaryConditions"):
            force.setUsesPeriodicBoundaryConditions(True)


def collagen_atom_indices(topology):
    indices = []
    heavy = []
    collagen_chains = []
    for chain in topology.chains():
        residues = list(chain.residues())
        if residues and residues[0].name not in {"HOH", "SOL", "WAT", "NA", "CL"}:
            collagen_chains.append(chain.index)
            for atom in chain.atoms():
                indices.append(atom.index)
                if atom.element is not None and atom.element.symbol != "H":
                    heavy.append(atom.index)
    if len(collagen_chains) != 3:
        raise RuntimeError(f"Expected three collagen chains, found {len(collagen_chains)}")
    return indices, heavy


def add_heavy_atom_restraints(system, reference_positions, heavy_indices):
    force = CustomExternalForce("0.5*k_restraint*periodicdistance(x,y,z,x0,y0,z0)^2")
    force.addGlobalParameter("k_restraint", 0.0)
    for name in ("x0", "y0", "z0"):
        force.addPerParticleParameter(name)
    ref_nm = reference_positions.value_in_unit(unit.nanometer)
    for index in heavy_indices:
        p = ref_nm[index]
        force.addParticle(index, [p.x, p.y, p.z])
    system.addForce(force)
    return force


def build_system(system_dir, gro, reference_positions, heavy_indices, barostat, temperature_k):
    top = GromacsTopFile(
        str(system_dir / "system.top"),
        periodicBoxVectors=gro.getPeriodicBoxVectors(),
        includeDir=str(system_dir),
    )
    system = top.createSystem(
        nonbondedMethod=PME,
        nonbondedCutoff=1.0 * unit.nanometer,
        constraints=HBonds,
        rigidWater=True,
        removeCMMotion=True,
    )
    enable_periodic_bonded_forces(system)
    add_heavy_atom_restraints(system, reference_positions, heavy_indices)
    if barostat:
        system.addForce(MonteCarloBarostat(1.0 * unit.bar, temperature_k * unit.kelvin, 25))
    return top, system


def platform_and_properties(name, device_index):
    platform = Platform.getPlatformByName(name)
    if name == "OpenCL":
        return platform, {
            "DeviceIndex": str(device_index),
            "Precision": "mixed",
            "UseCpuPme": "false",
        }
    if name == "CUDA":
        return platform, {"DeviceIndex": str(device_index), "Precision": "mixed"}
    if name == "CPU":
        return platform, {"Threads": "8"}
    return platform, {}


def new_simulation(topology, system, temperature_k, timestep_ps, integrator_seed, platform, properties):
    integrator = LangevinMiddleIntegrator(
        temperature_k * unit.kelvin,
        1.0 / unit.picosecond,
        timestep_ps * unit.picosecond,
    )
    integrator.setRandomNumberSeed(int(integrator_seed))
    return Simulation(topology, system, integrator, platform, properties)


def state_is_finite(state):
    energy = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    if not math.isfinite(energy):
        return False
    for p in state.getPositions().value_in_unit(unit.nanometer):
        if not (math.isfinite(p.x) and math.isfinite(p.y) and math.isfinite(p.z)):
            return False
    return True


def save_state(simulation, topology, output_dir, label):
    state = simulation.context.getState(
        getPositions=True, getVelocities=True, getEnergy=True, enforcePeriodicBox=True
    )
    if not state_is_finite(state):
        raise RuntimeError(f"Non-finite state after {label}")
    with (output_dir / f"{label}.pdb").open("w") as handle:
        PDBFile.writeFile(topology, state.getPositions(), handle, keepIds=True)
    (output_dir / f"{label}_state.xml").write_text(XmlSerializer.serialize(state) + "\n")
    (output_dir / f"{label}.chk").write_bytes(simulation.context.createCheckpoint())
    return state


def add_reporters(simulation, output_dir, label, steps, timestep_ps, trajectory_interval_ps, log_interval_ps):
    trajectory_steps = max(1, int(round(trajectory_interval_ps / timestep_ps)))
    log_steps = max(1, int(round(log_interval_ps / timestep_ps)))
    simulation.reporters.append(XTCReporter(str(output_dir / f"{label}.xtc"), trajectory_steps, enforcePeriodicBox=True))
    simulation.reporters.append(StateDataReporter(
        str(output_dir / f"{label}.csv"),
        log_steps,
        step=True,
        time=True,
        potentialEnergy=True,
        kineticEnergy=True,
        totalEnergy=True,
        temperature=True,
        volume=True,
        density=True,
        speed=True,
        separator=",",
    ))


def run_steps_checked(simulation, steps, chunk_steps):
    remaining = steps
    while remaining:
        take = min(remaining, chunk_steps)
        simulation.step(take)
        state = simulation.context.getState(getPositions=True, getEnergy=True)
        if not state_is_finite(state):
            raise RuntimeError("Non-finite coordinates or potential energy")
        remaining -= take


def run_attempt(repo_root, variant, replica, timestep_ps, args, attempt_dir):
    folder = SYSTEMS[variant]
    system_dir = repo_root / "systems" / "mmp1-free-apo" / folder
    input_gro = system_dir / f"MMP1_free_{variant}_150mM_NaCl_start.gro"
    gro = GromacsGroFile(str(input_gro))
    top_probe = GromacsTopFile(
        str(system_dir / "system.top"), periodicBoxVectors=gro.getPeriodicBoxVectors(), includeDir=str(system_dir)
    )
    _, heavy_indices = collagen_atom_indices(top_probe.topology)
    velocity_seed = args.master_seed + 1000 * list(SYSTEMS).index(variant) + 10 * replica
    integrator_seed_nvt = velocity_seed + 1
    integrator_seed_npt = velocity_seed + 2
    platform, properties = platform_and_properties(args.platform, args.device_index)

    top_nvt, system_nvt = build_system(system_dir, gro, gro.positions, heavy_indices, False, args.temperature_k)
    simulation = new_simulation(
        top_nvt.topology, system_nvt, 10.0, timestep_ps, integrator_seed_nvt, platform, properties
    )
    simulation.context.setPositions(gro.positions)
    simulation.context.setPeriodicBoxVectors(*gro.getPeriodicBoxVectors())
    simulation.context.setParameter("k_restraint", 1000.0 * KCAL_PER_MOL_A2_TO_KJ_PER_MOL_NM2)
    simulation.minimizeEnergy(
        tolerance=args.minimization_tolerance * unit.kilojoule_per_mole / unit.nanometer,
        maxIterations=args.minimization_max_iterations,
    )
    minimised = save_state(simulation, top_nvt.topology, attempt_dir, "00_minimised")

    simulation.context.setVelocitiesToTemperature(10.0 * unit.kelvin, velocity_seed)
    heat_steps = int(round(args.nvt_ps / timestep_ps))
    heating_chunks = max(1, args.heating_chunks)
    steps_per_chunk = heat_steps // heating_chunks
    remainder = heat_steps % heating_chunks
    add_reporters(simulation, attempt_dir, "01_nvt_heating", heat_steps, timestep_ps, args.trajectory_interval_ps, args.log_interval_ps)
    for chunk in range(heating_chunks):
        fraction = (chunk + 1) / heating_chunks
        temperature = 10.0 + fraction * (args.temperature_k - 10.0)
        simulation.integrator.setTemperature(temperature * unit.kelvin)
        simulation.step(steps_per_chunk + (1 if chunk < remainder else 0))
        check = simulation.context.getState(getPositions=True, getEnergy=True)
        if not state_is_finite(check):
            raise RuntimeError(f"Non-finite state during NVT heating chunk {chunk + 1}")
    nvt_state = save_state(simulation, top_nvt.topology, attempt_dir, "01_nvt_heating")

    top_npt, system_npt = build_system(system_dir, gro, gro.positions, heavy_indices, True, args.temperature_k)
    simulation = new_simulation(
        top_npt.topology, system_npt, args.temperature_k, timestep_ps, integrator_seed_npt, platform, properties
    )
    simulation.context.setPositions(nvt_state.getPositions())
    simulation.context.setVelocities(nvt_state.getVelocities())
    simulation.context.setPeriodicBoxVectors(*nvt_state.getPeriodicBoxVectors())

    stages = [
        ("02_npt_k1000", 1000.0, 50.0),
        ("03_npt_k500", 500.0, 50.0),
        ("04_npt_k100", 100.0, 50.0),
        ("05_npt_k10", 10.0, 50.0),
        ("06_npt_unrestrained", 0.0, 100.0),
    ]
    stage_summaries = []
    for label, restraint_k, duration_ps in stages:
        simulation.context.setParameter(
            "k_restraint", restraint_k * KCAL_PER_MOL_A2_TO_KJ_PER_MOL_NM2
        )
        steps = int(round(duration_ps / timestep_ps))
        simulation.reporters = []
        add_reporters(simulation, attempt_dir, label, steps, timestep_ps, args.trajectory_interval_ps, args.log_interval_ps)
        run_steps_checked(simulation, steps, max(1, int(round(5.0 / timestep_ps))))
        state = save_state(simulation, top_npt.topology, attempt_dir, label)
        stage_summaries.append({
            "label": label,
            "duration_ps": duration_ps,
            "restraint_k_kcal_mol_A2": restraint_k,
            "potential_energy_kj_mol": state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole),
            "kinetic_energy_kj_mol": state.getKineticEnergy().value_in_unit(unit.kilojoule_per_mole),
        })

    final_state = simulation.context.getState(
        getPositions=True, getVelocities=True, getEnergy=True, enforcePeriodicBox=True
    )
    final_gro = attempt_dir / f"MMP1_free_{variant}_replica_{replica}_equilibrated.gro"
    # OpenMM has no GRO writer; PDB + XML + checkpoint retain exact final state.
    with (attempt_dir / f"MMP1_free_{variant}_replica_{replica}_equilibrated.pdb").open("w") as handle:
        PDBFile.writeFile(top_npt.topology, final_state.getPositions(), handle, keepIds=True)

    return {
        "variant": variant,
        "replica": replica,
        "status": "completed",
        "input_gro": str(input_gro),
        "platform": args.platform,
        "device_index": args.device_index,
        "timestep_fs": timestep_ps * 1000.0,
        "temperature_k": args.temperature_k,
        "pressure_bar": 1.0,
        "velocity_seed": velocity_seed,
        "nvt_integrator_seed": integrator_seed_nvt,
        "npt_integrator_seed": integrator_seed_npt,
        "collagen_heavy_atoms_restrained": len(heavy_indices),
        "restraint_units": "kcal mol^-1 A^-2",
        "nvt_heating_ps": args.nvt_ps,
        "npt_stages": stage_summaries,
        "velocity_continuity_across_stages": True,
        "adaptive_sampling": False,
        "completed": datetime.now().isoformat(),
    }


def run_replica(repo_root, variant, replica, args):
    base = repo_root / "mmp1_free_conventional_md" / variant / f"replica_{replica}"
    base.mkdir(parents=True, exist_ok=True)
    completed_summary = base / "summary.json"
    if completed_summary.exists() and not args.force:
        previous = json.loads(completed_summary.read_text())
        if previous.get("status") == "completed":
            print(f"Skipping completed {variant} replica {replica}", flush=True)
            return previous

    errors = []
    for attempt, timestep_ps in enumerate(args.timestep_ps, 1):
        attempt_dir = base / f"attempt_{attempt}_{timestep_ps * 1000:g}fs"
        if attempt_dir.exists():
            shutil.rmtree(attempt_dir)
        attempt_dir.mkdir(parents=True)
        started = time.time()
        try:
            print(f"Running {variant} replica {replica} at {timestep_ps * 1000:g} fs", flush=True)
            summary = run_attempt(repo_root, variant, replica, timestep_ps, args, attempt_dir)
            summary["wall_time_seconds"] = time.time() - started
            summary["attempt"] = attempt
            summary["output_directory"] = str(attempt_dir)
            completed_summary.write_text(json.dumps(summary, indent=2) + "\n")
            return summary
        except Exception as exc:
            error = {
                "attempt": attempt,
                "timestep_fs": timestep_ps * 1000.0,
                "error": repr(exc),
                "traceback": traceback.format_exc(),
            }
            errors.append(error)
            (attempt_dir / "failure.json").write_text(json.dumps(error, indent=2) + "\n")
            print(f"Attempt failed: {exc}", flush=True)
    failed = {
        "variant": variant,
        "replica": replica,
        "status": "failed",
        "errors": errors,
        "completed": datetime.now().isoformat(),
    }
    completed_summary.write_text(json.dumps(failed, indent=2) + "\n")
    raise RuntimeError(f"All timestep attempts failed for {variant} replica {replica}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--variants", nargs="+", choices=list(SYSTEMS), default=list(SYSTEMS))
    parser.add_argument("--replicas", nargs="+", type=int, default=[1, 2, 3])
    parser.add_argument("--platform", choices=["OpenCL", "CUDA", "CPU"], default="OpenCL")
    parser.add_argument("--device-index", type=int, default=0)
    parser.add_argument("--temperature-k", type=float, default=310.15)
    parser.add_argument("--nvt-ps", type=float, default=20.0)
    parser.add_argument("--heating-chunks", type=int, default=20)
    parser.add_argument("--timestep-ps", nargs="+", type=float, default=[0.002, 0.001, 0.0005])
    parser.add_argument("--master-seed", type=int, default=202611140)
    parser.add_argument("--minimization-tolerance", type=float, default=10.0)
    parser.add_argument("--minimization-max-iterations", type=int, default=5000)
    parser.add_argument("--trajectory-interval-ps", type=float, default=10.0)
    parser.add_argument("--log-interval-ps", type=float, default=1.0)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    summaries = []
    for variant in args.variants:
        for replica in args.replicas:
            summaries.append(run_replica(args.repo_root, variant, replica, args))
    print(json.dumps(summaries, indent=2))


if __name__ == "__main__":
    main()
