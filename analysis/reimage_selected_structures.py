#!/usr/bin/env python3
"""Make selected MMP-1/collagen complexes whole across periodic boundaries."""

from __future__ import annotations

import argparse
from collections import defaultdict, deque
from pathlib import Path

import numpy as np

from run_peer_review_trajectory_analysis_standalone import (
    SYSTEMS,
    build_system_topology,
    pdb_element,
)


def read_gro(path: Path, expected_atoms: int) -> tuple[np.ndarray, np.ndarray]:
    lines = path.read_text(errors="ignore").splitlines()
    count = int(lines[1].strip())
    if count != expected_atoms:
        raise ValueError(f"{path} contains {count} atoms; expected {expected_atoms}")
    xyz = np.array(
        [[float(line[20:28]), float(line[28:36]), float(line[36:44])] for line in lines[2 : 2 + count]],
        dtype=float,
    )
    lengths = np.array([float(value) for value in lines[2 + count].split()[:3]], dtype=float)
    return xyz, lengths


def minimum_image(vector: np.ndarray, lengths: np.ndarray) -> np.ndarray:
    return vector - lengths * np.round(vector / lengths)


def unwrap_group(
    raw: np.ndarray,
    group_indices: set[int],
    bonds: list[tuple[int, int]],
    lengths: np.ndarray,
) -> np.ndarray:
    adjacency: dict[int, list[int]] = defaultdict(list)
    for first, second in bonds:
        if first in group_indices and second in group_indices:
            adjacency[first].append(second)
            adjacency[second].append(first)
    unwrapped = raw.copy()
    visited: set[int] = set()
    for start in sorted(group_indices):
        if start in visited:
            continue
        visited.add(start)
        queue = deque([start])
        while queue:
            current = queue.popleft()
            for neighbor in adjacency[current]:
                if neighbor in visited:
                    continue
                delta = minimum_image(raw[neighbor] - raw[current], lengths)
                unwrapped[neighbor] = unwrapped[current] + delta
                visited.add(neighbor)
                queue.append(neighbor)
    return unwrapped


def stitch_chain_by_residue(topology, raw: np.ndarray, corrected: np.ndarray, chain: str, lengths: np.ndarray) -> None:
    """Join any disconnected bonded components using consecutive collagen C-alpha atoms."""
    residues: dict[int, list[int]] = defaultdict(list)
    ca_index: dict[int, int] = {}
    for atom in topology.collagen_atoms:
        if atom.chain != chain:
            continue
        residues[int(atom.resid)].append(atom.index)
        if atom.name == "CA":
            ca_index[int(atom.resid)] = atom.index
    ordered = [resid for resid in sorted(residues) if resid in ca_index]
    for previous, current in zip(ordered, ordered[1:]):
        previous_ca = ca_index[previous]
        current_ca = ca_index[current]
        step = minimum_image(raw[current_ca] - raw[previous_ca], lengths)
        target = corrected[previous_ca] + step
        corrected[residues[current]] += target - corrected[current_ca]
    for resid in ordered:
        anchor = ca_index[resid]
        for atom_index in residues[resid]:
            corrected[atom_index] = corrected[anchor] + minimum_image(
                raw[atom_index] - raw[anchor], lengths
            )


def write_pdb(path: Path, topology, xyz_nm: np.ndarray, source: Path) -> None:
    chain_ids = {"mmp1": "A", "chain_1": "B", "chain_2": "C", "chain_3": "D"}
    lines = [f"REMARK Periodic-boundary-corrected selected structure from {source}"]
    serial = 1
    previous_chain = None
    for atom in sorted(topology.atoms, key=lambda item: item.index):
        chain_id = chain_ids[atom.chain]
        if previous_chain is not None and chain_id != previous_chain:
            lines.append("TER")
        x, y, z = xyz_nm[atom.index] * 10.0
        record = "HETATM" if atom.resname.upper() in {"ZN", "CA", "MG"} else "ATOM  "
        lines.append(
            f"{record}{serial:5d} {atom.name[:4].rjust(4)} {atom.resname[:3].rjust(3)} "
            f"{chain_id}{atom.resid:4d}    {x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{0.00:6.2f}"
            f"          {pdb_element(atom.name, atom.resname).rjust(2)}"
        )
        serial += 1
        previous_chain = chain_id
    lines.extend(["TER", "END"])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    for system in SYSTEMS:
        topology = build_system_topology(args.root.resolve(), system)
        source = args.input_dir / f"{system}_generation_028_selected.gro"
        raw, lengths = read_gro(source, topology.total_atoms)
        corrected = raw.copy()
        groups = {
            "mmp1": {atom.index for atom in topology.mmp_atoms},
            **{
                chain: {atom.index for atom in topology.collagen_atoms if atom.chain == chain}
                for chain in ("chain_1", "chain_2", "chain_3")
            },
        }
        for indices in groups.values():
            group_unwrapped = unwrap_group(raw, indices, topology.protein_bonds, lengths)
            corrected[list(indices)] = group_unwrapped[list(indices)]

        for chain in ("chain_1", "chain_2", "chain_3"):
            stitch_chain_by_residue(topology, raw, corrected, chain, lengths)

        mmp_center = corrected[list(groups["mmp1"])].mean(axis=0)
        for chain in ("chain_1", "chain_2", "chain_3"):
            indices = sorted(groups[chain])
            chain_center = corrected[indices].mean(axis=0)
            shift = lengths * np.round((mmp_center - chain_center) / lengths)
            corrected[indices] += shift

        write_pdb(
            args.output_dir / f"{system}_generation_028_selected.pdb",
            topology,
            corrected,
            source,
        )


if __name__ == "__main__":
    main()
