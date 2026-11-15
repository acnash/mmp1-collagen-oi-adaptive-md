#!/usr/bin/env python3
"""Audit completion and gross collagen geometry for the nine apo replicas."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import numpy as np
from openmm import XmlSerializer, unit
from openmm.app import GromacsGroFile, GromacsTopFile


VARIANTS = ["G978S", "G984C", "G987R"]


def positions_nm(quantity):
    return np.array([[p.x, p.y, p.z] for p in quantity.value_in_unit(unit.nanometer)])


def box_nm(vectors):
    return np.array([[v.x, v.y, v.z] for v in vectors.value_in_unit(unit.nanometer)])


def minimum_image_distance(a, b, box):
    delta = a - b
    fractional = delta @ np.linalg.inv(box)
    fractional -= np.rint(fractional)
    delta = fractional @ box
    return float(np.linalg.norm(delta))


def opening_metric(topology, positions, box):
    chains = []
    for chain in topology.chains():
        residues = list(chain.residues())
        if len(residues) == 33:
            chains.append(chain)
    if len(chains) != 3:
        raise RuntimeError(f"Expected three 33-residue collagen chains, found {len(chains)}")
    distances = []
    per_residue = []
    # Collagen topology spans model residue 963-995; 977-987 maps to indices 14-24.
    for residue_index in range(14, 25):
        ca = []
        for chain in chains:
            residue = list(chain.residues())[residue_index]
            atom = next(atom for atom in residue.atoms() if atom.name == "CA")
            ca.append(positions[atom.index])
        local = [
            minimum_image_distance(ca[0], ca[1], box),
            minimum_image_distance(ca[0], ca[2], box),
            minimum_image_distance(ca[1], ca[2], box),
        ]
        distances.extend(local)
        per_residue.append(max(local))
    return float(np.mean(distances)), float(np.max(per_residue))


def read_thermodynamics(path):
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    temperature_key = next(key for key in rows[0] if "Temperature" in key)
    density_key = next(key for key in rows[0] if "Density" in key)
    potential_key = next(key for key in rows[0] if "Potential Energy" in key)
    temperatures = np.array([float(row[temperature_key]) for row in rows])
    densities = np.array([float(row[density_key]) for row in rows])
    potentials = np.array([float(row[potential_key]) for row in rows])
    tail = slice(max(0, len(rows) // 2), None)
    return {
        "mean_temperature_k_last_50ps": float(np.mean(temperatures[tail])),
        "mean_density_g_ml_last_50ps": float(np.mean(densities[tail])),
        "potential_energy_finite": bool(np.all(np.isfinite(potentials))),
        "log_rows": len(rows),
    }


def main():
    repo = Path(__file__).resolve().parents[1]
    records = []
    for variant in VARIANTS:
        system_dir = repo / "systems" / "mmp1-free-apo" / f"collagen_{variant}"
        prep = json.loads((system_dir / "preparation_summary.json").read_text())
        gro = GromacsGroFile(str(system_dir / f"MMP1_free_{variant}_150mM_NaCl_start.gro"))
        top = GromacsTopFile(
            str(system_dir / "system.top"), periodicBoxVectors=gro.getPeriodicBoxVectors(), includeDir=str(system_dir)
        )
        start_mean, start_max = opening_metric(
            top.topology, positions_nm(gro.positions), box_nm(gro.getPeriodicBoxVectors())
        )
        for replica in (1, 2, 3):
            replica_dir = repo / "mmp1_free_conventional_md" / variant / f"replica_{replica}"
            summary = json.loads((replica_dir / "summary.json").read_text())
            attempt_dir = Path(summary["output_directory"])
            state = XmlSerializer.deserialize((attempt_dir / "06_npt_unrestrained_state.xml").read_text())
            final_mean, final_max = opening_metric(
                top.topology, positions_nm(state.getPositions()), box_nm(state.getPeriodicBoxVectors())
            )
            thermo = read_thermodynamics(attempt_dir / "06_npt_unrestrained.csv")
            record = {
                "variant": variant,
                "replica": replica,
                "status": summary["status"],
                "timestep_fs": summary["timestep_fs"],
                "velocity_seed": summary["velocity_seed"],
                "net_charge_e": prep["net_charge_e"],
                "nacl_pairs": prep["nacl_pairs"],
                "water_molecules": prep["water_molecules"],
                "start_opening_mean_nm": start_mean,
                "final_opening_mean_nm": final_mean,
                "opening_mean_change_nm": final_mean - start_mean,
                "start_max_local_pair_nm": start_max,
                "final_max_local_pair_nm": final_max,
                **thermo,
            }
            record["audit_pass"] = bool(
                record["status"] == "completed"
                and record["timestep_fs"] == 2.0
                and record["net_charge_e"] == 0
                and record["potential_energy_finite"]
                and 300.0 <= record["mean_temperature_k_last_50ps"] <= 320.0
                and 0.94 <= record["mean_density_g_ml_last_50ps"] <= 1.05
            )
            records.append(record)

    output_dir = repo / "mmp1_free_conventional_md"
    (output_dir / "audit_summary.json").write_text(json.dumps(records, indent=2) + "\n")
    with (output_dir / "audit_summary.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)
    print(json.dumps(records, indent=2))


if __name__ == "__main__":
    main()
