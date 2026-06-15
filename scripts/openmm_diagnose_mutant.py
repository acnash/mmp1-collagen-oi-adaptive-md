#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
from openmm import LangevinMiddleIntegrator, Platform, unit
from openmm.app import GromacsGroFile, GromacsTopFile, HBonds, PME, Simulation


SYSTEMS = {
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


def atom_label(topology, index: int) -> str:
    atom = list(topology.atoms())[index]
    return f"{index} 1based={index + 1} {atom.residue.name}{atom.residue.id}:{atom.name}"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("variant", choices=sorted(SYSTEMS))
    parser.add_argument("--repo_root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--steps", type=int, default=1000)
    parser.add_argument("--timestep_ps", type=float, default=0.002)
    parser.add_argument("--temperature_k", type=float, default=310.15)
    parser.add_argument("--constraints", choices=["HBonds", "None"], default="HBonds")
    args = parser.parse_args()

    folder, gro_name = SYSTEMS[args.variant]
    system_dir = args.repo_root / "generated_mutants" / "salted_150mM_NaCl" / folder
    gro = GromacsGroFile(str(system_dir / gro_name))
    top = GromacsTopFile(
        str(system_dir / "system.top"),
        periodicBoxVectors=gro.getPeriodicBoxVectors(),
        includeDir=str(system_dir),
    )
    system = top.createSystem(
        nonbondedMethod=PME,
        nonbondedCutoff=1.0 * unit.nanometer,
        constraints=None if args.constraints == "None" else HBonds,
        rigidWater=args.constraints != "None",
    )
    enable_periodic_bonded_forces(system)
    integrator = LangevinMiddleIntegrator(
        args.temperature_k * unit.kelvin,
        1.0 / unit.picosecond,
        args.timestep_ps * unit.picosecond,
    )
    integrator.setRandomNumberSeed(999)
    sim = Simulation(top.topology, system, integrator, Platform.getPlatformByName("CPU"), {"Threads": str(args.threads)})
    sim.context.setPositions(gro.positions)
    sim.context.setPeriodicBoxVectors(*gro.getPeriodicBoxVectors())
    sim.context.setVelocitiesToTemperature(args.temperature_k * unit.kelvin, 999)
    sim.minimizeEnergy(maxIterations=500)

    state = sim.context.getState(getPositions=True, getForces=True, getEnergy=True, enforcePeriodicBox=True)
    frc = state.getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole / unit.nanometer)
    norms = np.linalg.norm(frc, axis=1)
    print("after_min_potential_kj_mol", state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole))
    print("after_min_max_force", float(np.nanmax(norms)), atom_label(top.topology, int(np.nanargmax(norms))))
    print("nonfinite_forces", np.where(~np.isfinite(norms))[0].tolist())

    for i in range(args.steps):
        try:
            sim.step(1)
        except Exception as exc:
            print("exception_at_step", i + 1, repr(exc))
            state = sim.context.getState(getPositions=True, getVelocities=True, getForces=True, getEnergy=True, enforcePeriodicBox=True)
            pos = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
            vel = state.getVelocities(asNumpy=True).value_in_unit(unit.nanometer / unit.picosecond)
            frc = state.getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole / unit.nanometer)
            print("bad_pos", np.where(~np.isfinite(pos).all(axis=1))[0].tolist())
            print("bad_vel", np.where(~np.isfinite(vel).all(axis=1))[0].tolist())
            bad_frc = np.where(~np.isfinite(frc).all(axis=1))[0].tolist()
            print("bad_frc", bad_frc)
            norms = np.linalg.norm(frc, axis=1)
            if np.isfinite(norms).any():
                idx = int(np.nanargmax(norms))
                print("max_force_after_exception", float(np.nanmax(norms)), atom_label(top.topology, idx))
            raise
        if (i + 1) % 10 == 0:
            state = sim.context.getState(getForces=True, getEnergy=True)
            frc = state.getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole / unit.nanometer)
            norms = np.linalg.norm(frc, axis=1)
            idx = int(np.nanargmax(norms))
            print("step", i + 1, "energy", state.getPotentialEnergy(), "max_force", float(np.nanmax(norms)), atom_label(top.topology, idx), flush=True)


if __name__ == "__main__":
    main()
