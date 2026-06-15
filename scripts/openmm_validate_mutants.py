#!/usr/bin/env python3
"""Run short OpenMM minimization/NVT/NPT validation for prepared mutant systems."""

from __future__ import annotations

import argparse
import json
from datetime import datetime
from pathlib import Path

from openmm import LangevinMiddleIntegrator, MonteCarloBarostat, Platform, unit
from openmm.app import (
    GromacsGroFile,
    GromacsTopFile,
    HBonds,
    PME,
    Simulation,
    StateDataReporter,
)


MUTANT_SYSTEMS = {
    "G978S": ("collagen_G978S", "NPT_eq_collagen_G978S_150mM_NaCl.gro"),
    "G984C": ("collagen_G984C", "NPT_eq_collagen_G984C_150mM_NaCl.gro"),
    "G987R": ("collagen_G987R", "NPT_eq_collagen_G987R_150mM_NaCl.gro"),
}


def enable_periodic_bonded_forces(system):
    periodic_force_classes = {
        "HarmonicBondForce",
        "HarmonicAngleForce",
        "PeriodicTorsionForce",
        "RBTorsionForce",
        "CMAPTorsionForce",
        "CustomBondForce",
        "CustomAngleForce",
        "CustomTorsionForce",
    }
    for force in system.getForces():
        if force.__class__.__name__ in periodic_force_classes and hasattr(force, "setUsesPeriodicBoundaryConditions"):
            force.setUsesPeriodicBoundaryConditions(True)


def write_gro(path: Path, topology, positions, box_vectors, title: str) -> None:
    positions_nm = positions.value_in_unit(unit.nanometer)
    try:
        box_nm = box_vectors.value_in_unit(unit.nanometer)
    except AttributeError:
        box_nm = box_vectors
    with path.open("w") as handle:
        handle.write(f"{title}\n")
        handle.write(f"{len(list(topology.atoms())):5d}\n")
        for atom in topology.atoms():
            try:
                resid = int(atom.residue.id)
            except Exception:
                resid = atom.residue.index + 1
            pos = positions_nm[atom.index]
            handle.write(
                f"{resid % 100000:5d}"
                f"{atom.residue.name[:5]:<5}"
                f"{atom.name[:5]:>5}"
                f"{(atom.index + 1) % 100000:5d}"
                f"{pos.x:8.3f}{pos.y:8.3f}{pos.z:8.3f}\n"
            )
        a, b, c = box_nm
        handle.write(
            f"{float(a.x):10.5f}"
            f"{float(b.y):10.5f}"
            f"{float(c.z):10.5f}\n"
        )


def create_simulation(topology, system, threads: int, temperature_k: float, timestep_ps: float, seed: int, barostat: bool):
    if barostat:
        system.addForce(MonteCarloBarostat(1.0 * unit.bar, temperature_k * unit.kelvin, 25))
    integrator = LangevinMiddleIntegrator(temperature_k * unit.kelvin, 1.0 / unit.picosecond, timestep_ps * unit.picosecond)
    integrator.setRandomNumberSeed(seed)
    platform = Platform.getPlatformByName("CPU")
    properties = {"Threads": str(threads)}
    return Simulation(topology, system, integrator, platform, properties)


def build_system(system_dir: Path, gro_name: str):
    gro = GromacsGroFile(str(system_dir / gro_name))
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
    )
    enable_periodic_bonded_forces(system)
    return gro, top, system


def run_one(repo_root: Path, variant: str, threads: int, nvt_ps: float, npt_ps: float, timestep_ps: float, temperature_k: float) -> dict:
    folder, gro_name = MUTANT_SYSTEMS[variant]
    system_dir = repo_root / "generated_mutants" / "salted_150mM_NaCl" / folder
    validation_dir = system_dir / "openmm_validation"
    validation_dir.mkdir(exist_ok=True)

    steps_nvt = int(round(nvt_ps / timestep_ps))
    steps_npt = int(round(npt_ps / timestep_ps))
    started = datetime.now().isoformat()

    gro, top, system_nvt = build_system(system_dir, gro_name)
    sim_nvt = create_simulation(top.topology, system_nvt, threads, temperature_k, timestep_ps, seed=1101, barostat=False)
    sim_nvt.context.setPositions(gro.positions)
    sim_nvt.context.setPeriodicBoxVectors(*gro.getPeriodicBoxVectors())
    sim_nvt.context.setVelocitiesToTemperature(temperature_k * unit.kelvin, 1101)
    sim_nvt.reporters.append(
        StateDataReporter(
            str(validation_dir / "nvt.log"),
            max(1, steps_nvt // 10),
            step=True,
            time=True,
            potentialEnergy=True,
            temperature=True,
            speed=True,
        )
    )

    sim_nvt.minimizeEnergy(maxIterations=500)
    state_min = sim_nvt.context.getState(getPositions=True, getEnergy=True, enforcePeriodicBox=True)
    write_gro(validation_dir / f"{variant}_after_minimization.gro", top.topology, state_min.getPositions(), state_min.getPeriodicBoxVectors(), f"{variant} after OpenMM minimization")

    sim_nvt.step(steps_nvt)
    state_nvt = sim_nvt.context.getState(getPositions=True, getVelocities=True, getEnergy=True, enforcePeriodicBox=True)
    write_gro(validation_dir / f"{variant}_after_10ps_nvt.gro", top.topology, state_nvt.getPositions(), state_nvt.getPeriodicBoxVectors(), f"{variant} after 10 ps NVT")

    _, _, system_npt = build_system(system_dir, gro_name)
    sim_npt = create_simulation(top.topology, system_npt, threads, temperature_k, timestep_ps, seed=2202, barostat=True)
    sim_npt.context.setPositions(state_nvt.getPositions())
    sim_npt.context.setPeriodicBoxVectors(*state_nvt.getPeriodicBoxVectors())
    sim_npt.context.setVelocities(state_nvt.getVelocities())
    sim_npt.reporters.append(
        StateDataReporter(
            str(validation_dir / "npt.log"),
            max(1, steps_npt // 10),
            step=True,
            time=True,
            potentialEnergy=True,
            kineticEnergy=True,
            totalEnergy=True,
            temperature=True,
            volume=True,
            density=True,
            speed=True,
        )
    )

    sim_npt.step(steps_npt)
    state_npt = sim_npt.context.getState(getPositions=True, getEnergy=True, enforcePeriodicBox=True)
    final_name = f"NPT_eq_collagen_{variant}_150mM_NaCl_openmm_10ps_npt.gro"
    final_path = system_dir / final_name
    write_gro(final_path, top.topology, state_npt.getPositions(), state_npt.getPeriodicBoxVectors(), f"{variant} after OpenMM minimization, 10 ps NVT, 10 ps NPT")
    write_gro(validation_dir / final_name, top.topology, state_npt.getPositions(), state_npt.getPeriodicBoxVectors(), f"{variant} after OpenMM minimization, 10 ps NVT, 10 ps NPT")

    summary = {
        "variant": variant,
        "started": started,
        "completed": datetime.now().isoformat(),
        "threads": threads,
        "platform": "CPU",
        "temperature_k": temperature_k,
        "timestep_ps": timestep_ps,
        "nvt_ps": nvt_ps,
        "npt_ps": npt_ps,
        "input_gro": gro_name,
        "final_gro": final_name,
        "minimized_potential_kj_per_mol": state_min.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole),
        "post_nvt_potential_kj_per_mol": state_nvt.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole),
        "post_npt_potential_kj_per_mol": state_npt.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole),
    }
    (validation_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description="Validate prepared mutant systems with short OpenMM CPU minimization/NVT/NPT runs.")
    parser.add_argument("--repo_root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--variants", nargs="+", choices=sorted(MUTANT_SYSTEMS), default=sorted(MUTANT_SYSTEMS))
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--nvt_ps", type=float, default=10.0)
    parser.add_argument("--npt_ps", type=float, default=10.0)
    parser.add_argument("--timestep_ps", type=float, default=0.002)
    parser.add_argument("--temperature_k", type=float, default=310.15)
    args = parser.parse_args()

    summaries = []
    for variant in args.variants:
        print(f"Validating {variant}...", flush=True)
        summary = run_one(args.repo_root, variant, args.threads, args.nvt_ps, args.npt_ps, args.timestep_ps, args.temperature_k)
        summaries.append(summary)
        print(f"Completed {variant}: {summary['final_gro']}", flush=True)

    print(json.dumps(summaries, indent=2))


if __name__ == "__main__":
    main()
