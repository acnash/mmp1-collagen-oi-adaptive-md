#!/usr/bin/env python3
"""Extract compact structural data from selected adaptive unwinding frames.

The script analyses the best_opening_frame.gro selected for each completed adaptive
generation. It does not read or copy trajectories. Outputs are CSV files containing
residue-resolved collagen opening geometry, MMP-1 contacts, and frame summaries.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd


SYSTEMS = {
    "wild_type": {
        "label": "Wild type",
        "run_prefix": "adaptive_run_wild_type_",
        "relative_dir": Path("systems/mmp1-bound-adaptive/wild_type"),
    },
    "G978S": {
        "label": "G978S",
        "run_prefix": "adaptive_run_G978S_",
        "relative_dir": Path("systems/mmp1-bound-adaptive/collagen_G978S"),
    },
    "G984C": {
        "label": "G984C",
        "run_prefix": "adaptive_run_G984C_",
        "relative_dir": Path("systems/mmp1-bound-adaptive/collagen_G984C"),
    },
    "G987R": {
        "label": "G987R",
        "run_prefix": "adaptive_run_G987R_",
        "relative_dir": Path("systems/mmp1-bound-adaptive/collagen_G987R"),
    },
}

WINDOW_RESNAMES = {
    977: "ALA",
    978: "GLY",
    979: "GLN",
    980: "ARG",
    981: "GLY",
    982: "ILE",
    983: "VAL",
    984: "GLY",
    985: "LEU",
    986: "HYP",
    987: "GLY",
}

MUTATIONS = {
    "G978S": (978, "SER"),
    "G984C": (984, "CYS"),
    "G987R": (987, "ARG"),
}


@dataclass
class Atom:
    index0: int
    atom_number: int
    resid: int
    resname: str
    name: str
    xyz_a: np.ndarray
    chain: Optional[str] = None


@dataclass
class FrameRecord:
    system: str
    system_label: str
    run_dir: Path
    generation_index: int
    opening_score_angstrom: float
    best_gro: Path
    worker_dir: str


def read_json(path: Path) -> dict:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {}


def find_latest_run(system_dir: Path, prefix: str) -> Optional[Path]:
    runs = sorted(path for path in system_dir.glob(f"{prefix}*") if path.is_dir())
    return runs[-1].resolve() if runs else None


def resolve_moved_path(path_like: Optional[str], run_dir: Path) -> Optional[Path]:
    if not path_like:
        return None
    path = Path(path_like)
    if path.exists():
        return path.resolve()
    text = str(path_like)
    marker = "adaptive_run_"
    marker_index = text.find(marker)
    if marker_index >= 0:
        original = Path(text[marker_index:])
        repaired = run_dir.joinpath(*original.parts[1:])
        if repaired.exists():
            return repaired.resolve()
    return None


def collect_frames(system: str, run_dir: Path) -> List[FrameRecord]:
    label = SYSTEMS[system]["label"]
    context = read_json(run_dir / "context.json")
    rows: List[FrameRecord] = []
    for generation in context.get("selected_frames", []):
        generation_index = generation.get("generation")
        candidates = generation.get("frames", [])
        if generation_index is None or not candidates:
            continue
        selected = sorted(candidates, key=lambda row: row.get("rank", 999))[0]
        best_gro = resolve_moved_path(selected.get("best_gro"), run_dir)
        if best_gro is None:
            continue
        rows.append(
            FrameRecord(
                system=system,
                system_label=label,
                run_dir=run_dir,
                generation_index=int(generation_index),
                opening_score_angstrom=float(selected.get("best_score_angstrom", np.nan)),
                best_gro=best_gro,
                worker_dir=str(selected.get("worker_dir", "")),
            )
        )
    return sorted(rows, key=lambda row: row.generation_index)


def atom_count_from_itp(path: Path) -> int:
    in_atoms = False
    count = 0
    for line in path.read_text(errors="ignore").splitlines():
        text = line.strip()
        if not text or text.startswith(";"):
            continue
        if text.startswith("["):
            in_atoms = text.lower().startswith("[ atoms")
            continue
        if in_atoms and text.split()[0].isdigit():
            count += 1
    return count


def molecule_name_from_itp(path: Path) -> Optional[str]:
    in_moleculetype = False
    for line in path.read_text(errors="ignore").splitlines():
        text = line.strip()
        if not text or text.startswith(";"):
            continue
        if text.startswith("["):
            in_moleculetype = text.lower().startswith("[ moleculetype")
            continue
        if in_moleculetype:
            return text.split()[0]
    return None


def molecule_lines(system_top: Path) -> List[Tuple[str, int]]:
    rows: List[Tuple[str, int]] = []
    in_molecules = False
    for line in system_top.read_text(errors="ignore").splitlines():
        text = line.strip()
        if not text or text.startswith(";"):
            continue
        if text.startswith("["):
            in_molecules = text.lower().startswith("[ molecules")
            continue
        if in_molecules:
            fields = text.split()
            if len(fields) >= 2 and fields[1].isdigit():
                rows.append((fields[0], int(fields[1])))
    return rows


def molecule_ranges(input_dir: Path) -> Dict[str, List[Tuple[int, int, str]]]:
    counts: Dict[str, int] = {}
    for itp in input_dir.glob("*.itp"):
        if itp.name in {"forcefield.itp", "ffbonded.itp", "ffnonbonded.itp", "tip3p.itp", "ions.itp"}:
            continue
        molecule_name = molecule_name_from_itp(itp)
        if molecule_name:
            counts[molecule_name] = atom_count_from_itp(itp)

    current = 1
    collagen_index = 0
    ranges: Dict[str, List[Tuple[int, int, str]]] = {"mmp1": [], "collagen": []}
    for molecule, copies in molecule_lines(input_dir / "system.top"):
        atom_count = counts.get(molecule)
        if atom_count is None:
            continue
        for _ in range(copies):
            start = current
            end = current + atom_count - 1
            lower = molecule.lower()
            if lower.startswith("mmp1"):
                ranges["mmp1"].append((start, end, molecule))
            elif lower.startswith("collagen"):
                collagen_index += 1
                ranges["collagen"].append((start, end, f"chain_{collagen_index}"))
            current = end + 1
    if len(ranges["mmp1"]) != 1 or len(ranges["collagen"]) != 3:
        raise RuntimeError(f"Could not identify one MMP-1 and three collagen molecules in {input_dir}")
    return ranges


def parse_gro(path: Path) -> Tuple[List[Atom], np.ndarray]:
    lines = path.read_text(errors="ignore").splitlines()
    atom_count = int(lines[1].strip())
    atoms: List[Atom] = []
    for index, line in enumerate(lines[2 : 2 + atom_count]):
        atoms.append(
            Atom(
                index0=index,
                atom_number=int(line[15:20]),
                resid=int(line[0:5]),
                resname=line[5:10].strip(),
                name=line[10:15].strip(),
                xyz_a=np.array(
                    [float(line[20:28]), float(line[28:36]), float(line[36:44])],
                    dtype=float,
                )
                * 10.0,
            )
        )
    values = [float(value) * 10.0 for value in lines[2 + atom_count].split()]
    if len(values) == 3:
        box = np.diag(values)
    elif len(values) == 9:
        box = np.array(
            [
                [values[0], values[3], values[4]],
                [values[5], values[1], values[6]],
                [values[7], values[8], values[2]],
            ],
            dtype=float,
        )
    else:
        raise ValueError(f"Unsupported GRO box with {len(values)} values in {path}")
    return atoms, box


def assign_molecules(atoms: List[Atom], ranges: Dict[str, List[Tuple[int, int, str]]]) -> Tuple[List[Atom], List[Atom]]:
    mmp_range = ranges["mmp1"][0]
    mmp_atoms = [atom for atom in atoms if mmp_range[0] <= atom.atom_number <= mmp_range[1]]
    collagen_atoms: List[Atom] = []
    for atom in atoms:
        for start, end, chain in ranges["collagen"]:
            if start <= atom.atom_number <= end:
                atom.chain = chain
                collagen_atoms.append(atom)
                break
    return mmp_atoms, collagen_atoms


def group_residues(atoms: Iterable[Atom]) -> Dict[str, List[List[Atom]]]:
    result: Dict[str, List[List[Atom]]] = {}
    for chain in ["chain_1", "chain_2", "chain_3"]:
        chain_atoms = [atom for atom in atoms if atom.chain == chain]
        residues: List[List[Atom]] = []
        current: List[Atom] = []
        current_key = None
        for atom in chain_atoms:
            key = (atom.resid, atom.resname)
            if current_key is None:
                current_key = key
            if key != current_key:
                residues.append(current)
                current = []
                current_key = key
            current.append(atom)
        if current:
            residues.append(current)
        result[chain] = residues
    return result


def expected_names(system: str, chain: str) -> List[str]:
    names = [WINDOW_RESNAMES[resid] for resid in sorted(WINDOW_RESNAMES)]
    if system in MUTATIONS and chain == "chain_1":
        resid, replacement = MUTATIONS[system]
        names[resid - min(WINDOW_RESNAMES)] = replacement
    return names


def find_window(residues: Sequence[List[Atom]], expected: Sequence[str]) -> List[List[Atom]]:
    length = len(expected)
    for start in range(len(residues) - length + 1):
        candidate = residues[start : start + length]
        if [residue[0].resname for residue in candidate] == list(expected):
            return list(candidate)
    by_id = [residue for residue in residues if min(WINDOW_RESNAMES) <= residue[0].resid <= max(WINDOW_RESNAMES)]
    if len(by_id) == length:
        return by_id
    raise RuntimeError(f"Could not identify collagen residues 977-987 from sequence {expected}")


def minimum_image(delta: np.ndarray, box: np.ndarray) -> np.ndarray:
    fractional = np.asarray(delta, dtype=float) @ np.linalg.inv(box)
    fractional -= np.round(fractional)
    return fractional @ box


def pbc_distance(a: np.ndarray, b: np.ndarray, box: np.ndarray) -> float:
    return float(np.linalg.norm(minimum_image(b - a, box)))


def atom_by_name(residue: Sequence[Atom], name: str) -> Atom:
    for atom in residue:
        if atom.name == name:
            return atom
    raise RuntimeError(f"Residue {residue[0].resid} {residue[0].resname} has no {name} atom")


def unwrap_ca_window(windows: Dict[str, List[List[Atom]]], box: np.ndarray) -> np.ndarray:
    chains: List[np.ndarray] = []
    for chain in ["chain_1", "chain_2", "chain_3"]:
        raw = np.array([atom_by_name(residue, "CA").xyz_a for residue in windows[chain]])
        unwrapped = [raw[0]]
        for point in raw[1:]:
            unwrapped.append(unwrapped[-1] + minimum_image(point - unwrapped[-1], box))
        chains.append(np.array(unwrapped))
    reference = chains[0].mean(axis=0)
    for index in [1, 2]:
        shift = minimum_image(chains[index].mean(axis=0) - reference, box)
        chains[index] = chains[index] + (reference + shift - chains[index].mean(axis=0))
    return np.stack(chains, axis=0)


def triangle_area(a: np.ndarray, b: np.ndarray, c: np.ndarray, box: np.ndarray) -> float:
    ab = minimum_image(b - a, box)
    ac = minimum_image(c - a, box)
    return float(0.5 * np.linalg.norm(np.cross(ab, ac)))


def kabsch_rmsd(reference: np.ndarray, mobile: np.ndarray) -> float:
    reference_centered = reference - reference.mean(axis=0)
    mobile_centered = mobile - mobile.mean(axis=0)
    covariance = mobile_centered.T @ reference_centered
    u, _, vt = np.linalg.svd(covariance)
    correction = np.eye(3)
    correction[-1, -1] = np.sign(np.linalg.det(u @ vt))
    rotation = u @ correction @ vt
    aligned = mobile_centered @ rotation
    return float(np.sqrt(np.mean(np.sum((aligned - reference_centered) ** 2, axis=1))))


def heavy(atoms: Iterable[Atom]) -> List[Atom]:
    return [atom for atom in atoms if not atom.name.upper().startswith("H")]


def residue_mmp_contacts(
    residue: Sequence[Atom],
    mmp_xyz: np.ndarray,
    box: np.ndarray,
    cutoff_a: float,
) -> Tuple[int, float]:
    residue_xyz = np.array([atom.xyz_a for atom in heavy(residue)], dtype=float)
    pair_count = 0
    minimum = float("inf")
    inverse_box = np.linalg.inv(box)
    for point in residue_xyz:
        delta = mmp_xyz - point
        fractional = delta @ inverse_box
        fractional -= np.round(fractional)
        distances = np.linalg.norm(fractional @ box, axis=1)
        pair_count += int(np.count_nonzero(distances <= cutoff_a))
        minimum = min(minimum, float(distances.min()))
    return pair_count, minimum


def analyse_frame(
    frame: FrameRecord,
    ranges: Dict[str, List[Tuple[int, int, str]]],
    reference_ca: Optional[np.ndarray],
    contact_cutoff_a: float,
) -> Tuple[List[dict], List[dict], dict, np.ndarray]:
    atoms, box = parse_gro(frame.best_gro)
    mmp_atoms, collagen_atoms = assign_molecules(atoms, ranges)
    residues = group_residues(collagen_atoms)
    windows = {
        chain: find_window(residues[chain], expected_names(frame.system, chain))
        for chain in ["chain_1", "chain_2", "chain_3"]
    }
    ca = unwrap_ca_window(windows, box)
    mmp_xyz = np.array([atom.xyz_a for atom in heavy(mmp_atoms)], dtype=float)

    geometry_rows: List[dict] = []
    contact_rows: List[dict] = []
    all_distances: List[float] = []
    all_areas: List[float] = []
    all_contacts = 0
    min_mmp_distance = float("inf")

    residue_numbers = list(range(min(WINDOW_RESNAMES), max(WINDOW_RESNAMES) + 1))
    for position, canonical_resid in enumerate(residue_numbers):
        points = [atom_by_name(windows[chain][position], "CA").xyz_a for chain in ["chain_1", "chain_2", "chain_3"]]
        distances = [
            pbc_distance(points[0], points[1], box),
            pbc_distance(points[0], points[2], box),
            pbc_distance(points[1], points[2], box),
        ]
        area = triangle_area(points[0], points[1], points[2], box)
        all_distances.extend(distances)
        all_areas.append(area)
        geometry_rows.append(
            {
                "system": frame.system,
                "system_label": frame.system_label,
                "generation_index": frame.generation_index,
                "opening_score_angstrom": frame.opening_score_angstrom,
                "residue_position": canonical_resid,
                "chain_1_resname": windows["chain_1"][position][0].resname,
                "chain_2_resname": windows["chain_2"][position][0].resname,
                "chain_3_resname": windows["chain_3"][position][0].resname,
                "chain_1_chain_2_ca_distance_angstrom": distances[0],
                "chain_1_chain_3_ca_distance_angstrom": distances[1],
                "chain_2_chain_3_ca_distance_angstrom": distances[2],
                "mean_ca_distance_angstrom": float(np.mean(distances)),
                "maximum_ca_distance_angstrom": float(np.max(distances)),
                "distance_asymmetry_angstrom": float(np.max(distances) - np.min(distances)),
                "ca_triangle_area_angstrom2": area,
                "best_gro": str(frame.best_gro),
            }
        )
        for chain in ["chain_1", "chain_2", "chain_3"]:
            residue = windows[chain][position]
            contact_count, minimum = residue_mmp_contacts(residue, mmp_xyz, box, contact_cutoff_a)
            all_contacts += contact_count
            min_mmp_distance = min(min_mmp_distance, minimum)
            contact_rows.append(
                {
                    "system": frame.system,
                    "system_label": frame.system_label,
                    "generation_index": frame.generation_index,
                    "opening_score_angstrom": frame.opening_score_angstrom,
                    "chain": chain,
                    "residue_position": canonical_resid,
                    "resname": residue[0].resname,
                    "heavy_atom_contact_pairs": contact_count,
                    "contact_cutoff_angstrom": contact_cutoff_a,
                    "minimum_mmp1_distance_angstrom": minimum,
                    "best_gro": str(frame.best_gro),
                }
            )

    flattened_ca = ca.reshape(-1, 3)
    radius_of_gyration = float(np.sqrt(np.mean(np.sum((flattened_ca - flattened_ca.mean(axis=0)) ** 2, axis=1))))
    rmsd = 0.0 if reference_ca is None else kabsch_rmsd(reference_ca.reshape(-1, 3), flattened_ca)
    summary = {
        "system": frame.system,
        "system_label": frame.system_label,
        "run_dir": str(frame.run_dir),
        "generation_index": frame.generation_index,
        "opening_score_angstrom": frame.opening_score_angstrom,
        "mean_window_ca_distance_angstrom": float(np.mean(all_distances)),
        "maximum_window_ca_distance_angstrom": float(np.max(all_distances)),
        "mean_window_triangle_area_angstrom2": float(np.mean(all_areas)),
        "maximum_window_triangle_area_angstrom2": float(np.max(all_areas)),
        "window_ca_radius_of_gyration_angstrom": radius_of_gyration,
        "window_ca_rmsd_from_system_generation_0_angstrom": rmsd,
        "window_mmp1_heavy_atom_contact_pairs": all_contacts,
        "minimum_window_mmp1_distance_angstrom": min_mmp_distance,
        "contact_cutoff_angstrom": contact_cutoff_a,
        "best_gro": str(frame.best_gro),
        "worker_dir": frame.worker_dir,
    }
    return geometry_rows, contact_rows, summary, ca


def resolve_runs(args: argparse.Namespace) -> Dict[str, Path]:
    root = args.root.resolve()
    explicit = {
        "wild_type": args.wild_type_run,
        "G978S": args.g978s_run,
        "G984C": args.g984c_run,
        "G987R": args.g987r_run,
    }
    runs: Dict[str, Path] = {}
    for system, meta in SYSTEMS.items():
        if explicit[system]:
            run_dir = explicit[system].expanduser().resolve()
        else:
            run_dir = find_latest_run(root / meta["relative_dir"], meta["run_prefix"])
        if run_dir is not None and run_dir.exists():
            runs[system] = run_dir
    if not runs:
        raise SystemExit("No adaptive run directories were found")
    return runs


def input_dir_for_run(run_dir: Path) -> Path:
    context = read_json(run_dir / "context.json")
    configured = context.get("input_dir")
    if configured and Path(configured).exists():
        return Path(configured).resolve()
    return run_dir.parent.resolve()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument("--wild_type_run", type=Path)
    parser.add_argument("--g978s_run", type=Path)
    parser.add_argument("--g984c_run", type=Path)
    parser.add_argument("--g987r_run", type=Path)
    parser.add_argument("--analysis_dir", type=Path, default=Path("compact_structural_analysis"))
    parser.add_argument("--contact_cutoff_angstrom", type=float, default=4.5)
    args = parser.parse_args()

    output_dir = args.analysis_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    geometry_rows: List[dict] = []
    contact_rows: List[dict] = []
    summary_rows: List[dict] = []

    for system, run_dir in resolve_runs(args).items():
        frames = collect_frames(system, run_dir)
        if not frames:
            print(f"No selected frames found for {system}: {run_dir}")
            continue
        ranges = molecule_ranges(input_dir_for_run(run_dir))
        reference_ca: Optional[np.ndarray] = None
        for frame in frames:
            geometry, contacts, summary, ca = analyse_frame(
                frame,
                ranges,
                reference_ca,
                args.contact_cutoff_angstrom,
            )
            if reference_ca is None:
                reference_ca = ca.copy()
            geometry_rows.extend(geometry)
            contact_rows.extend(contacts)
            summary_rows.append(summary)
            print(f"{system} generation {frame.generation_index:03d} analysed")

    geometry_df = pd.DataFrame(geometry_rows)
    contacts_df = pd.DataFrame(contact_rows)
    summary_df = pd.DataFrame(summary_rows)
    geometry_df.to_csv(output_dir / "residue_resolved_opening_geometry.csv", index=False)
    contacts_df.to_csv(output_dir / "residue_resolved_mmp1_contacts.csv", index=False)
    summary_df.to_csv(output_dir / "structural_frame_summary.csv", index=False)

    if not summary_df.empty:
        baseline = summary_df[summary_df["generation_index"] == 0][
            [
                "system",
                "window_mmp1_heavy_atom_contact_pairs",
                "mean_window_triangle_area_angstrom2",
                "window_ca_radius_of_gyration_angstrom",
            ]
        ].rename(
            columns={
                "window_mmp1_heavy_atom_contact_pairs": "baseline_mmp1_contact_pairs",
                "mean_window_triangle_area_angstrom2": "baseline_mean_triangle_area_angstrom2",
                "window_ca_radius_of_gyration_angstrom": "baseline_window_ca_radius_of_gyration_angstrom",
            }
        )
        changes = summary_df.merge(baseline, on="system", how="left")
        changes["change_in_mmp1_contact_pairs"] = (
            changes["window_mmp1_heavy_atom_contact_pairs"] - changes["baseline_mmp1_contact_pairs"]
        )
        changes["change_in_mean_triangle_area_angstrom2"] = (
            changes["mean_window_triangle_area_angstrom2"]
            - changes["baseline_mean_triangle_area_angstrom2"]
        )
        changes["change_in_window_ca_radius_of_gyration_angstrom"] = (
            changes["window_ca_radius_of_gyration_angstrom"]
            - changes["baseline_window_ca_radius_of_gyration_angstrom"]
        )
        changes.to_csv(output_dir / "structural_changes_from_generation_zero.csv", index=False)

    print(f"CSV outputs written to {output_dir}")


if __name__ == "__main__":
    main()
