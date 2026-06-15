#!/usr/bin/env python3
"""Run final mutant NVT/NPT equilibration and save coordinates with velocities."""

from __future__ import annotations

import argparse
import json
from datetime import datetime
from pathlib import Path

from openmm import LangevinMiddleIntegrator, MonteCarloBarostat, Platform, XmlSerializer, unit
from openmm.app import GromacsGroFile, GromacsTopFile, HBonds, PME, Simulation, StateDataReporter


SYSTEMS = {
    "G978S": (
        "collagen_G978S",
        "NPT_eq_collagen_G978S_150mM_NaCl_openmm_10ps_npt.gro",
    ),
    "G984C": (
        "collagen_G984C",
        "NPT_eq_collagen_G984C_150mM_NaCl_openmm_1ps_nvt_1ps_npt.gro",
    ),
    "G987R": (
        "collagen_G987R",
        "NPT_eq_collagen_G987R_150mM_NaCl_openmm_1ps_nvt_1ps_npt.gro",
    ),
}


def enable_periodic_bonded_forces(system) -> None:
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


def write_gro_with_velocities(path: Path, topology, positions, velocities, box_vectors, title: str) -> None:
    positions_nm = positions.value_in_unit(unit.nanometer)
    velocities_nm_ps = velocities.value_in_unit(unit.nanometer / unit.picosecond)
    box_nm = box_vectors.value_in_unit(unit.nanometer)
    atoms = list(topology.atoms())
    with path.open("w") as handle:
        handle.write(f"{title}\n")
        handle.write(f"{len(atoms):5d}\n")
        for atom in atoms:
            try:
                resid = int(atom.residue.id)
            except Exception:
                resid = atom.residue.index + 1
            pos = positions_nm[atom.index]
            vel = velocities_nm_ps[atom.index]
            handle.write(
                f"{resid % 100000:5d}"
                f"{atom.residue.name[:5]:<5}"
                f"{atom.name[:5]:>5}"
                f"{(atom.index + 1) % 100000:5d}"
                f"{pos.x:8.3f}{pos.y:8.3f}{pos.z:8.3f}"
                f"{vel.x:8.4f}{vel.y:8.4f}{vel.z:8.4f}\n"
            )
        a, b, c = box_nm
        handle.write(f"{float(a.x):10.5f}{float(b.y):10.5f}{float(c.z):10.5f}\n")


def build_system(system_dir: Path, gro_name: str, constraints: str):
    gro = GromacsGroFile(str(system_dir / gro_name))
    top = GromacsTopFile(
        str(system_dir / "system.top"),
        periodicBoxVectors=gro.getPeriodicBoxVectors(),
        includeDir=str(system_dir),
    )
    system = top.createSystem(
        nonbondedMethod=PME,
        nonbondedCutoff=1.0 * unit.nanometer,
        constraints=None if constraints == "None" else HBonds,
        rigidWater=constraints != "None",
    )
    enable_periodic_bonded_forces(system)
    return gro, top, system


def create_simulation(topology, system, threads: int, temperature_k: float, timestep_ps: float, seed: int) -> Simulation:
    integrator = LangevinMiddleIntegrator(
        temperature_k * unit.kelvin,
        1.0 / unit.picosecond,
        timestep_ps * unit.picosecond,
    )
    integrator.setRandomNumberSeed(seed)
    return Simulation(topology, system, integrator, Platform.getPlatformByName("CPU"), {"Threads": str(threads)})


def write_failure(output_dir: Path, variant: str, stage: str, exc: Exception, summary: dict) -> None:
    summary.update(
        {
            "variant": variant,
            "status": "failed",
            "failed_stage": stage,
            "error": repr(exc),
            "completed": datetime.now().isoformat(),
        }
    )
    (output_dir / "failure_summary.json").write_text(json.dumps(summary, indent=2) + "\n")


def run_variant(repo_root: Path, variant: str, threads: int, timestep_ps: float, temperature_k: float, nvt_ps: float, npt_ps: float, constraints: str) -> dict:
    folder, gro_name = SYSTEMS[variant]
    system_dir = repo_root / "generated_mutants" / "salted_150mM_NaCl" / folder
    output_dir = system_dir / "openmm_final_10ps_10ps_with_velocities"
    output_dir.mkdir(exist_ok=True)
    steps_nvt = int(round(nvt_ps / timestep_ps))
    steps_npt = int(round(npt_ps / timestep_ps))

    summary = {
        "variant": variant,
        "started": datetime.now().isoformat(),
        "input_gro": gro_name,
        "threads": threads,
        "platform": "CPU",
        "temperature_k": temperature_k,
        "timestep_ps": timestep_ps,
        "nvt_ps": nvt_ps,
        "npt_ps": npt_ps,
        "steps_nvt": steps_nvt,
        "steps_npt": steps_npt,
        "constraints": constraints,
        "rigid_water": constraints != "None",
        "outputs_include_velocities": True,
    }

    gro, top, system_nvt = build_system(system_dir, gro_name, constraints)
    sim_nvt = create_simulation(top.topology, system_nvt, threads, temperature_k, timestep_ps, seed=7101)
    sim_nvt.context.setPositions(gro.positions)
    sim_nvt.context.setPeriodicBoxVectors(*gro.getPeriodicBoxVectors())
    sim_nvt.context.setVelocitiesToTemperature(temperature_k * unit.kelvin, 7101)
    sim_nvt.reporters.append(
        StateDataReporter(
            str(output_dir / "nvt.log"),
            max(1, steps_nvt // 10),
            step=True,
            time=True,
            potentialEnergy=True,
            kineticEnergy=True,
            totalEnergy=True,
            temperature=True,
            speed=True,
        )
    )
    try:
        sim_nvt.step(steps_nvt)
    except Exception as exc:
        write_failure(output_dir, variant, "NVT", exc, summary)
        raise

    state_nvt = sim_nvt.context.getState(getPositions=True, getVelocities=True, getEnergy=True, enforcePeriodicBox=True)
    write_gro_with_velocities(
        output_dir / f"{variant}_after_10ps_nvt_with_velocities.gro",
        top.topology,
        state_nvt.getPositions(),
        state_nvt.getVelocities(),
        state_nvt.getPeriodicBoxVectors(),
        f"{variant} after 10 ps NVT with velocities",
    )

    _, _, system_npt = build_system(system_dir, gro_name, constraints)
    system_npt.addForce(MonteCarloBarostat(1.0 * unit.bar, temperature_k * unit.kelvin, 25))
    sim_npt = create_simulation(top.topology, system_npt, threads, temperature_k, timestep_ps, seed=7202)
    sim_npt.context.setPositions(state_nvt.getPositions())
    sim_npt.context.setPeriodicBoxVectors(*state_nvt.getPeriodicBoxVectors())
    sim_npt.context.setVelocities(state_nvt.getVelocities())
    sim_npt.reporters.append(
        StateDataReporter(
            str(output_dir / "npt.log"),
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
    try:
        sim_npt.step(steps_npt)
    except Exception as exc:
        write_failure(output_dir, variant, "NPT", exc, summary)
        raise

    state_npt = sim_npt.context.getState(getPositions=True, getVelocities=True, getEnergy=True, enforcePeriodicBox=True)
    final_name = f"NPT_eq_collagen_{variant}_150mM_NaCl_openmm_10ps_nvt_10ps_npt_with_velocities.gro"
    write_gro_with_velocities(
        system_dir / final_name,
        top.topology,
        state_npt.getPositions(),
        state_npt.getVelocities(),
        state_npt.getPeriodicBoxVectors(),
        f"{variant} after 10 ps NVT and 10 ps NPT with velocities",
    )
    write_gro_with_velocities(
        output_dir / final_name,
        top.topology,
        state_npt.getPositions(),
        state_npt.getVelocities(),
        state_npt.getPeriodicBoxVectors(),
        f"{variant} after 10 ps NVT and 10 ps NPT with velocities",
    )
    (output_dir / f"{variant}_final_state.xml").write_text(XmlSerializer.serialize(state_npt) + "\n")
    (output_dir / f"{variant}_final_checkpoint.chk").write_bytes(sim_npt.context.createCheckpoint())

    summary.update(
        {
            "status": "completed",
            "completed": datetime.now().isoformat(),
            "final_gro": final_name,
            "post_nvt_potential_kj_per_mol": state_nvt.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole),
            "post_npt_potential_kj_per_mol": state_npt.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole),
        }
    )
    (output_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description="Run final 10 ps NVT + 10 ps NPT and save velocities.")
    parser.add_argument("--repo_root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--variants", nargs="+", choices=sorted(SYSTEMS), default=sorted(SYSTEMS))
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--timestep_ps", type=float, default=0.002)
    parser.add_argument("--temperature_k", type=float, default=312.5)
    parser.add_argument("--nvt_ps", type=float, default=10.0)
    parser.add_argument("--npt_ps", type=float, default=10.0)
    parser.add_argument("--constraints", choices=["HBonds", "None"], default="HBonds")
    args = parser.parse_args()

    summaries = []
    for variant in args.variants:
        print(f"Final equilibration for {variant}...", flush=True)
        summaries.append(
            run_variant(
                args.repo_root,
                variant,
                args.threads,
                args.timestep_ps,
                args.temperature_k,
                args.nvt_ps,
                args.npt_ps,
                args.constraints,
            )
        )
        print(f"Completed {variant}", flush=True)
    print(json.dumps(summaries, indent=2))


if __name__ == "__main__":
    main()
