#!/usr/bin/env python3
"""Build MMP-1-free, periodically bonded collagen systems at 150 mM NaCl.

The starting collagen coordinates and force-field parameters are taken from the
compact, pre-adaptive MMP-1/collagen systems.  MMP-1 and all original solvent
are discarded, the unchanged periodic collagen topology is retained, and the
original periodic box is freshly solvated with OpenMM's equilibrated TIP3P
water box.  Ions are then introduced by replacing water molecules to give
approximately 150 mM NaCl and neutralise the collagen mutation.
"""

from __future__ import annotations

import argparse
import json
import math
import random
import shutil
from pathlib import Path

import numpy as np
from openmm import NonbondedForce, System, unit
from openmm.app import GromacsGroFile, GromacsTopFile, Modeller, PME, Topology


AVOGADRO_PER_NM3_MOLAR = 0.602214076
TARGET_SALT_MOLAR = 0.150

SYSTEMS = {
    "G978S": ("collagen_G978S", "NPT_eq_collagen_G978S_150mM_NaCl.gro", "collagen_G978S"),
    "G984C": ("collagen_G984C", "NPT_eq_collagen_G984C_150mM_NaCl.gro", "collagen_G984C"),
    "G987R": ("collagen_G987R", "NPT_eq_collagen_G987R_150mM_NaCl.gro", "collagen_G987R"),
}


class CollisionForceField:
    """Minimal ForceField facade used only by Modeller.addSolvent()."""

    def __init__(self, particle_parameters):
        self.particle_parameters = particle_parameters

    def createSystem(self, topology, **kwargs):
        if topology.getNumAtoms() != len(self.particle_parameters):
            raise ValueError("Collision system particle count does not match collagen topology")
        system = System()
        nonbonded = NonbondedForce()
        for charge, sigma, epsilon in self.particle_parameters:
            system.addParticle(12.0)
            nonbonded.addParticle(charge, sigma, epsilon)
        system.addForce(nonbonded)
        return system


def copy_collagen_topology(source_topology, source_positions, selected_chain_indices):
    selected_chains = {i for i in selected_chain_indices}
    topology = Topology()
    topology.setPeriodicBoxVectors(source_topology.getPeriodicBoxVectors())
    atom_map = {}
    selected_old_atoms = []
    positions = []
    for chain in source_topology.chains():
        if chain.index not in selected_chains:
            continue
        new_chain = topology.addChain(chain.id)
        for residue in chain.residues():
            new_residue = topology.addResidue(residue.name, new_chain, residue.id, residue.insertionCode)
            for atom in residue.atoms():
                new_atom = topology.addAtom(atom.name, atom.element, new_residue, atom.id, atom.formalCharge)
                atom_map[atom] = new_atom
                selected_old_atoms.append(atom)
                positions.append(source_positions[atom.index])
    for bond in source_topology.bonds():
        atom1, atom2 = bond[0], bond[1]
        if atom1 in atom_map and atom2 in atom_map:
            topology.addBond(atom_map[atom1], atom_map[atom2])
    return topology, unit.Quantity(positions), selected_old_atoms


def box_matrix_nm(box_vectors):
    return np.array([[v.x, v.y, v.z] for v in box_vectors.value_in_unit(unit.nanometer)], dtype=float)


def minimum_image_distances(points, point, box, inverse_box):
    delta = points - point
    fractional = delta @ inverse_box
    fractional -= np.rint(fractional)
    delta = fractional @ box
    return np.sqrt(np.sum(delta * delta, axis=1))


def choose_ion_waters(water_oxygens, solute_heavy, count, box, seed, minimum_distance_nm=0.5):
    rng = random.Random(seed)
    inverse_box = np.linalg.inv(box)
    candidates = list(range(len(water_oxygens)))
    rng.shuffle(candidates)
    selected = []
    selected_positions = []
    for index in candidates:
        pos = water_oxygens[index]
        if len(solute_heavy) and np.min(minimum_image_distances(solute_heavy, pos, box, inverse_box)) < minimum_distance_nm:
            continue
        if selected_positions:
            chosen = np.asarray(selected_positions)
            if np.min(minimum_image_distances(chosen, pos, box, inverse_box)) < minimum_distance_nm:
                continue
        selected.append(index)
        selected_positions.append(pos)
        if len(selected) == count:
            return selected
    raise RuntimeError(f"Could select only {len(selected)} of {count} ion sites")


def write_gro(path, title, collagen_atoms, collagen_positions, waters, sodium_positions, chloride_positions, box_vectors):
    records = []
    residue_number = 1
    for atom, position in zip(collagen_atoms, collagen_positions.value_in_unit(unit.nanometer)):
        records.append((residue_number, atom.residue.name, atom.name, position))
        if atom.index == list(atom.residue.atoms())[-1].index:
            residue_number += 1
    for water in waters:
        for name, position in water:
            records.append((residue_number, "SOL", name, position))
        residue_number += 1
    for position in sodium_positions:
        records.append((residue_number, "NA", "NA", position))
        residue_number += 1
    for position in chloride_positions:
        records.append((residue_number, "CL", "CL", position))
        residue_number += 1

    with path.open("w") as handle:
        handle.write(title + "\n")
        handle.write(f"{len(records):5d}\n")
        for atom_number, (resnr, resname, atomname, position) in enumerate(records, 1):
            handle.write(
                f"{resnr % 100000:5d}{resname[:5]:<5}{atomname[:5]:>5}{atom_number % 100000:5d}"
                f"{position[0]:8.3f}{position[1]:8.3f}{position[2]:8.3f}\n"
            )
        box = box_matrix_nm(box_vectors)
        if np.max(np.abs(box - np.diag(np.diag(box)))) < 1e-10:
            handle.write(f"{box[0, 0]:10.5f}{box[1, 1]:10.5f}{box[2, 2]:10.5f}\n")
        else:
            handle.write(
                f"{box[0,0]:10.5f}{box[1,1]:10.5f}{box[2,2]:10.5f}"
                f"{box[0,1]:10.5f}{box[0,2]:10.5f}{box[1,0]:10.5f}"
                f"{box[1,2]:10.5f}{box[2,0]:10.5f}{box[2,1]:10.5f}\n"
            )


def write_topology(path, mutant_molecule, n_waters, n_na, n_cl):
    path.write_text(
        '#include "forcefield.itp"\n\n'
        f'#include "{mutant_molecule}.itp"\n'
        '#include "collagen.itp"\n\n'
        '#include "tip3p.itp"\n'
        '#include "ions.itp"\n\n'
        '[ system ]\n'
        'MMP-1-free periodic mutant collagen in 150 mM NaCl\n\n'
        '[ molecules ]\n'
        f'{mutant_molecule:<20} 1\n'
        f'{"collagen":<20} 2\n'
        f'{"SOL":<20} {n_waters}\n'
        f'{"NA":<20} {n_na}\n'
        f'{"CL":<20} {n_cl}\n'
    )


def prepare_variant(repo_root, variant, seed):
    folder, gro_name, mutant_molecule = SYSTEMS[variant]
    source_dir = repo_root / "generated_mutants" / "salted_150mM_NaCl" / folder
    output_dir = repo_root / "generated_mutants" / "mmp1_free_150mM_NaCl" / folder
    output_dir.mkdir(parents=True, exist_ok=True)

    gro = GromacsGroFile(str(source_dir / gro_name))
    top = GromacsTopFile(
        str(source_dir / "system.top"),
        periodicBoxVectors=gro.getPeriodicBoxVectors(),
        includeDir=str(source_dir),
    )
    full_system = top.createSystem(nonbondedMethod=PME, nonbondedCutoff=1.0 * unit.nanometer)
    nonbonded = next(force for force in full_system.getForces() if isinstance(force, NonbondedForce))

    # Chain 0 is MMP-1. Chains 1-3 are the three collagen strands/molecules.
    collagen_topology, collagen_positions, old_atoms = copy_collagen_topology(
        top.topology, gro.positions, [1, 2, 3]
    )
    particle_parameters = [nonbonded.getParticleParameters(atom.index) for atom in old_atoms]
    collagen_charge = sum(
        q.value_in_unit(unit.elementary_charge) for q, sigma, epsilon in particle_parameters
    )
    rounded_charge = int(round(collagen_charge))
    if abs(collagen_charge - rounded_charge) > 1e-3:
        raise RuntimeError(f"Non-integral collagen charge for {variant}: {collagen_charge}")

    modeller = Modeller(collagen_topology, collagen_positions)
    modeller.addSolvent(
        CollisionForceField(particle_parameters),
        model="tip3p",
        boxVectors=gro.getPeriodicBoxVectors(),
        ionicStrength=0 * unit.molar,
        neutralize=False,
    )

    all_atoms = list(modeller.topology.atoms())
    positions_nm = np.array(
        [[p.x, p.y, p.z] for p in modeller.positions.value_in_unit(unit.nanometer)], dtype=float
    )
    collagen_atom_count = len(old_atoms)
    water_residues = [r for r in modeller.topology.residues() if r.name in {"HOH", "SOL", "WAT"}]
    waters = []
    water_oxygens = []
    for residue in water_residues:
        residue_atoms = list(residue.atoms())
        oxygen = next(a for a in residue_atoms if a.element.symbol == "O")
        hydrogens = [a for a in residue_atoms if a.element.symbol == "H"]
        waters.append([
            ("OW", positions_nm[oxygen.index]),
            ("HW1", positions_nm[hydrogens[0].index]),
            ("HW2", positions_nm[hydrogens[1].index]),
        ])
        water_oxygens.append(positions_nm[oxygen.index])
    water_oxygens = np.asarray(water_oxygens)

    box = box_matrix_nm(gro.getPeriodicBoxVectors())
    volume_nm3 = abs(np.linalg.det(box))
    salt_pairs = int(round(TARGET_SALT_MOLAR * volume_nm3 * AVOGADRO_PER_NM3_MOLAR))
    n_na = salt_pairs + max(0, -rounded_charge)
    n_cl = salt_pairs + max(0, rounded_charge)
    n_ions = n_na + n_cl
    solute_heavy = np.array([
        positions_nm[a.index] for a in all_atoms[:collagen_atom_count] if a.element.symbol != "H"
    ])
    selected = choose_ion_waters(water_oxygens, solute_heavy, n_ions, box, seed)
    ion_set = set(selected)
    sodium_sites = selected[:n_na]
    chloride_sites = selected[n_na:]
    remaining_waters = [water for i, water in enumerate(waters) if i not in ion_set]

    for filename in ["forcefield.itp", "ffbonded.itp", "ffnonbonded.itp", "tip3p.itp", "ions.itp", "collagen.itp", f"{mutant_molecule}.itp"]:
        shutil.copy2(source_dir / filename, output_dir / filename)
    write_topology(output_dir / "system.top", mutant_molecule, len(remaining_waters), n_na, n_cl)
    output_gro = output_dir / f"MMP1_free_{variant}_150mM_NaCl_start.gro"
    write_gro(
        output_gro,
        f"MMP-1-free {variant}, freshly solvated at 150 mM NaCl",
        old_atoms,
        collagen_positions,
        remaining_waters,
        [water_oxygens[i] for i in sodium_sites],
        [water_oxygens[i] for i in chloride_sites],
        gro.getPeriodicBoxVectors(),
    )

    summary = {
        "variant": variant,
        "source_gro": str(source_dir / gro_name),
        "output_gro": str(output_gro),
        "mmp1_atoms_removed": 5825,
        "collagen_atoms": collagen_atom_count,
        "collagen_charge_e": rounded_charge,
        "box_volume_nm3": volume_nm3,
        "target_nacl_molar": TARGET_SALT_MOLAR,
        "nacl_pairs": salt_pairs,
        "water_molecules": len(remaining_waters),
        "sodium_ions": n_na,
        "chloride_ions": n_cl,
        "net_charge_e": rounded_charge + n_na - n_cl,
        "periodic_collagen_topology_retained": True,
        "fresh_tip3p_resolvation": True,
        "ion_site_seed": seed,
    }
    (output_dir / "preparation_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--variants", nargs="+", choices=sorted(SYSTEMS), default=sorted(SYSTEMS))
    parser.add_argument("--seed", type=int, default=20261114)
    args = parser.parse_args()
    summaries = [prepare_variant(args.repo_root, variant, args.seed + i) for i, variant in enumerate(args.variants)]
    print(json.dumps(summaries, indent=2))


if __name__ == "__main__":
    main()
