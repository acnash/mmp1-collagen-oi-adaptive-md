#!/usr/bin/env python3
"""Compare MMP-1 interdomain and catalytic-opening distances in adaptive MD.

This standalone reader is intended to be placed in ``scripts/`` and launched
from the repository root.  It analyses the wild-type, G978S, G984C and G987R
adaptive-MD XTC trajectories without modifying them.

The three reported coordinates reproduce the definitions established in the
MMP1 restraints project:

* mass-weighted catalytic-domain to HPX-domain centre-of-mass distance;
* SER123 CA to SER347 CA distance;
* catalytic-opening proxy: biological N171 CA to T230 CA, represented in the
  coordinate topology by ASN152 CA and THR211 CA (a +19 biological offset).

Domain definitions are coordinate residues 81-250 (catalytic) and 266-447
(HPX).  The finite MMP-1 molecule is reconstructed through covalent bonds and
minimum-image bond vectors before its domain COM and SER distance are measured.
The catalytic-opening vector is measured directly by the periodic minimum
image.  Distances are written in both nm and Angstrom.

Only NumPy and OpenMM are required for the CSV analysis.  Matplotlib is
optional; when available, comparison plots are also written.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
import subprocess
import sys
from collections import defaultdict, deque
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Optional

import numpy as np

OPENMM_IMPORT_ERROR = None
try:
    import openmm  # noqa: F401 - import verifies the analysis environment
    from openmm.app.internal.xtc_utils import get_xtc_natoms, read_xtc
except Exception as exc:  # pragma: no cover - server dependency check
    openmm = None
    get_xtc_natoms = None
    read_xtc = None
    OPENMM_IMPORT_ERROR = exc


SYSTEMS = {
    "wild_type": {
        "label": "Wild type",
        "relative_dir": Path("generated_mutants/salted_150mM_NaCl/wild_type"),
        "run_prefix": "adaptive_run_wild_type_",
    },
    "G978S": {
        "label": "G978S",
        "relative_dir": Path("generated_mutants/salted_150mM_NaCl/collagen_G978S"),
        "run_prefix": "adaptive_run_G978S_",
    },
    "G984C": {
        "label": "G984C",
        "relative_dir": Path("generated_mutants/salted_150mM_NaCl/collagen_G984C"),
        "run_prefix": "adaptive_run_G984C_",
    },
    "G987R": {
        "label": "G987R",
        "relative_dir": Path("generated_mutants/salted_150mM_NaCl/collagen_G987R"),
        "run_prefix": "adaptive_run_G987R_",
    },
}

CATALYTIC_RANGE = (81, 250)
HPX_RANGE = (266, 447)
LANDMARKS = {
    "ser123": (123, "SER", "CA"),
    "ser347": (347, "SER", "CA"),
    # Published/biological N171 and T230 correspond to coordinate IDs -19.
    "n171": (152, "ASN", "CA"),
    "t230": (211, "THR", "CA"),
}

METRICS_NM = (
    "catalytic_hpx_com_nm",
    "ser123_ser347_ca_nm",
    "catalytic_opening_n171_t230_ca_nm",
)

METRIC_LABELS = {
    "catalytic_hpx_com_nm": "Catalytic-HPX COM",
    "ser123_ser347_ca_nm": "SER123-SER347 CA",
    "catalytic_opening_n171_t230_ca_nm": "N171-T230 catalytic opening",
}


@dataclass(frozen=True)
class LocalAtom:
    nr: int
    atom_type: str
    resid: int
    resname: str
    name: str
    mass: float


@dataclass(frozen=True)
class MoleculeType:
    name: str
    atoms: tuple[LocalAtom, ...]
    bonds: tuple[tuple[int, int], ...]


@dataclass(frozen=True)
class AtomMeta:
    index: int
    molecule: str
    molecule_copy: int
    category: str
    resid: int
    resname: str
    name: str
    atom_type: str
    mass: float


@dataclass
class SystemTopology:
    system: str
    input_dir: Path
    input_gro: Path
    mmp_atoms: list[AtomMeta]
    mmp_bonds: list[tuple[int, int]]
    total_atoms: int
    catalytic_indices: np.ndarray
    catalytic_masses: np.ndarray
    hpx_indices: np.ndarray
    hpx_masses: np.ndarray
    landmark_indices: dict[str, int]
    unwrap_parent_order: list[tuple[int, int]]


class RunningStats:
    """Numerically stable population statistics."""

    def __init__(self) -> None:
        self.n = 0
        self.mean = 0.0
        self.m2 = 0.0
        self.minimum = float("inf")
        self.maximum = float("-inf")

    def add(self, value: float) -> None:
        value = float(value)
        if not np.isfinite(value):
            raise ValueError(f"Non-finite metric value: {value}")
        self.n += 1
        delta = value - self.mean
        self.mean += delta / self.n
        self.m2 += delta * (value - self.mean)
        self.minimum = min(self.minimum, value)
        self.maximum = max(self.maximum, value)

    @property
    def sd(self) -> float:
        return float(np.sqrt(self.m2 / self.n)) if self.n else float("nan")

    def as_fields(self, prefix: str) -> dict[str, float | int]:
        return {
            f"{prefix}_n": self.n,
            f"{prefix}_mean": self.mean if self.n else float("nan"),
            f"{prefix}_sd": self.sd,
            f"{prefix}_min": self.minimum if self.n else float("nan"),
            f"{prefix}_max": self.maximum if self.n else float("nan"),
        }


def parse_int_list(text: str) -> list[int]:
    values: list[int] = []
    for part in text.split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part:
            first, last = part.split("-", 1)
            start, stop = int(first), int(last)
            if stop < start:
                raise ValueError(f"Invalid descending generation range: {part}")
            values.extend(range(start, stop + 1))
        else:
            values.append(int(part))
    if not values:
        raise ValueError("At least one generation must be supplied")
    return sorted(set(values))


def sha256(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(chunk_size), b""):
            digest.update(block)
    return digest.hexdigest()


def git_commit(root: Path) -> Optional[str]:
    result = subprocess.run(
        ["git", "-C", str(root), "rev-parse", "HEAD"],
        check=False,
        capture_output=True,
        text=True,
    )
    return result.stdout.strip() if result.returncode == 0 else None


def read_json(path: Path) -> dict:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {}


def strip_comment(line: str) -> str:
    return line.split(";", 1)[0].strip()


def parse_itp_molecule_types(path: Path) -> dict[str, MoleculeType]:
    raw: dict[str, dict] = {}
    section = ""
    current: Optional[str] = None
    awaiting_molecule_name = False
    for original in path.read_text(errors="ignore").splitlines():
        line = strip_comment(original)
        if not line or line.startswith("#"):
            continue
        if line.startswith("[") and line.endswith("]"):
            section = line.strip("[] ").lower()
            awaiting_molecule_name = section == "moleculetype"
            continue
        fields = line.split()
        if awaiting_molecule_name:
            current = fields[0]
            raw[current] = {"atoms": [], "bonds": []}
            awaiting_molecule_name = False
            continue
        if current is None:
            continue
        if section == "atoms" and len(fields) >= 7 and fields[0].isdigit():
            try:
                atom_name = fields[4]
                inferred_mass = 1.008 if atom_name.upper().startswith("H") else 12.0
                raw[current]["atoms"].append(
                    LocalAtom(
                        nr=int(fields[0]),
                        atom_type=fields[1],
                        resid=int(fields[2]),
                        resname=fields[3].upper(),
                        name=atom_name.upper(),
                        mass=float(fields[7]) if len(fields) >= 8 else inferred_mass,
                    )
                )
            except ValueError:
                continue
        elif section == "bonds" and len(fields) >= 2 and fields[0].isdigit() and fields[1].isdigit():
            raw[current]["bonds"].append((int(fields[0]), int(fields[1])))
    result = {}
    for name, values in raw.items():
        atoms = tuple(sorted(values["atoms"], key=lambda atom: atom.nr))
        if atoms:
            result[name] = MoleculeType(name, atoms, tuple(values["bonds"]))
    return result


def parse_molecule_table(system_top: Path) -> list[tuple[str, int]]:
    rows: list[tuple[str, int]] = []
    in_molecules = False
    for original in system_top.read_text(errors="ignore").splitlines():
        line = strip_comment(original)
        if not line:
            continue
        if line.startswith("[") and line.endswith("]"):
            in_molecules = line.strip("[] ").lower() == "molecules"
            continue
        if in_molecules:
            fields = line.split()
            if len(fields) >= 2:
                try:
                    rows.append((fields[0], int(fields[1])))
                except ValueError:
                    pass
    return rows


def choose_input_gro(input_dir: Path) -> Path:
    candidates = sorted(input_dir.glob("*.gro"))
    preferred = [p for p in candidates if "150mM_NaCl" in p.name and "openmm" not in p.name]
    if preferred:
        return preferred[0]
    if candidates:
        return candidates[0]
    raise FileNotFoundError(f"No input GRO file found in {input_dir}")


def unique_atom_index(atoms: list[AtomMeta], resid: int, resname: str, name: str, label: str) -> int:
    matches = [a.index for a in atoms if a.resid == resid and a.resname == resname and a.name == name]
    if len(matches) != 1:
        observed = sorted({(a.resid, a.resname, a.name) for a in atoms if a.resid == resid and a.name == name})
        raise RuntimeError(
            f"{label}: expected exactly one {resname}{resid} {name}; found {len(matches)}. "
            f"Observed CA identity at that coordinate ID: {observed}"
        )
    return matches[0]


def build_unwrap_order(mmp_atoms: list[AtomMeta], bonds: list[tuple[int, int]], required: set[int]) -> list[tuple[int, int]]:
    indices = {atom.index for atom in mmp_atoms}
    adjacency: dict[int, list[int]] = defaultdict(list)
    for first, second in bonds:
        if first in indices and second in indices:
            adjacency[first].append(second)
            adjacency[second].append(first)
    if not required:
        raise RuntimeError("No required MMP-1 atoms were selected")
    anchor = min(required)
    visited = {anchor}
    queue = deque([anchor])
    order: list[tuple[int, int]] = []
    while queue:
        parent = queue.popleft()
        for child in adjacency[parent]:
            if child not in visited:
                visited.add(child)
                order.append((parent, child))
                queue.append(child)
    missing = sorted(required - visited)
    if missing:
        raise RuntimeError(
            f"Covalent MMP-1 graph does not connect {len(missing)} required domain/landmark atoms; "
            f"first missing indices: {missing[:10]}"
        )
    return order


def build_system_topology(root: Path, system: str) -> SystemTopology:
    input_dir = (root / SYSTEMS[system]["relative_dir"]).resolve()
    system_top = input_dir / "system.top"
    if not system_top.exists():
        raise FileNotFoundError(system_top)
    molecule_types: dict[str, MoleculeType] = {}
    for itp in sorted(input_dir.glob("*.itp")):
        molecule_types.update(parse_itp_molecule_types(itp))

    mmp_atoms: list[AtomMeta] = []
    mmp_bonds: list[tuple[int, int]] = []
    offset = 0
    counters: dict[str, int] = defaultdict(int)
    mmp_copies = 0
    for molecule_name, copies in parse_molecule_table(system_top):
        if molecule_name not in molecule_types:
            raise RuntimeError(f"No ITP atom definition found for {molecule_name} in {input_dir}")
        molecule = molecule_types[molecule_name]
        if [a.nr for a in molecule.atoms] != list(range(1, len(molecule.atoms) + 1)):
            raise RuntimeError(f"Non-contiguous atom numbering in molecule {molecule_name}")
        for _ in range(copies):
            counters[molecule_name] += 1
            is_mmp = molecule_name.lower().startswith("mmp1")
            if is_mmp:
                mmp_copies += 1
                for local in molecule.atoms:
                    mmp_atoms.append(
                        AtomMeta(
                            index=offset + local.nr - 1,
                            molecule=molecule_name,
                            molecule_copy=counters[molecule_name],
                            category="mmp1",
                            resid=local.resid,
                            resname=local.resname,
                            name=local.name,
                            atom_type=local.atom_type,
                            mass=local.mass,
                        )
                    )
                mmp_bonds.extend((offset + a - 1, offset + b - 1) for a, b in molecule.bonds)
            offset += len(molecule.atoms)
    if mmp_copies != 1 or not mmp_atoms:
        raise RuntimeError(f"Expected one MMP-1 molecule in {system_top}; found {mmp_copies}")

    input_gro = choose_input_gro(input_dir)
    gro_lines = input_gro.read_text(errors="ignore").splitlines()
    gro_count = int(gro_lines[1].strip())
    if gro_count != offset:
        raise RuntimeError(f"Topology has {offset} atoms but {input_gro} has {gro_count}")

    catalytic = [a for a in mmp_atoms if CATALYTIC_RANGE[0] <= a.resid <= CATALYTIC_RANGE[1]]
    hpx = [a for a in mmp_atoms if HPX_RANGE[0] <= a.resid <= HPX_RANGE[1]]
    if not catalytic or not hpx:
        raise RuntimeError(f"Could not construct catalytic and HPX selections for {system}")
    if any(a.mass <= 0 for a in catalytic + hpx):
        raise RuntimeError(f"Non-positive atom mass in a domain selection for {system}")

    landmarks = {
        label: unique_atom_index(mmp_atoms, resid, resname, name, f"{system} {label}")
        for label, (resid, resname, name) in LANDMARKS.items()
    }
    catalytic_indices = np.asarray([a.index for a in catalytic], dtype=int)
    hpx_indices = np.asarray([a.index for a in hpx], dtype=int)
    required = set(catalytic_indices) | set(hpx_indices) | set(landmarks.values())
    unwrap_order = build_unwrap_order(mmp_atoms, mmp_bonds, required)
    return SystemTopology(
        system=system,
        input_dir=input_dir,
        input_gro=input_gro,
        mmp_atoms=mmp_atoms,
        mmp_bonds=mmp_bonds,
        total_atoms=offset,
        catalytic_indices=catalytic_indices,
        catalytic_masses=np.asarray([a.mass for a in catalytic], dtype=float),
        hpx_indices=hpx_indices,
        hpx_masses=np.asarray([a.mass for a in hpx], dtype=float),
        landmark_indices=landmarks,
        unwrap_parent_order=unwrap_order,
    )


def find_latest_run(input_dir: Path, prefix: str, generations: Iterable[int]) -> Path:
    runs = sorted(path for path in input_dir.glob(f"{prefix}*") if path.is_dir())
    if not runs:
        raise FileNotFoundError(f"No adaptive run beginning {prefix!r} found in {input_dir}")
    required = list(generations)
    complete = [r for r in runs if all((r / f"generation_{g:03d}").is_dir() for g in required)]
    return (complete[-1] if complete else runs[-1]).resolve()


def resolve_recorded_path(path_like: str, run_dir: Path) -> Optional[Path]:
    path = Path(path_like)
    if path.exists():
        return path.resolve()
    text = str(path_like)
    position = text.find("adaptive_run_")
    if position >= 0:
        original = Path(text[position:])
        repaired = run_dir.joinpath(*original.parts[1:])
        if repaired.exists():
            return repaired.resolve()
    return None


def worker_xtc_segments(worker_dir: Path) -> list[Path]:
    meta = read_json(worker_dir / "worker.json")
    recorded = []
    for value in meta.get("archived_xtc_segments", []):
        repaired = resolve_recorded_path(str(value), worker_dir.parents[1])
        if repaired is not None and repaired.suffix == ".xtc":
            recorded.append(repaired)
    fallback = sorted((worker_dir / "previous_xtc_segments").glob("*.xtc"))
    segments: list[Path] = []
    seen: set[Path] = set()
    for path in recorded + fallback:
        resolved = path.resolve()
        if resolved not in seen:
            segments.append(resolved)
            seen.add(resolved)
    current = worker_dir / "production.xtc"
    if current.exists() and current.resolve() not in seen:
        segments.append(current.resolve())
    return segments


def normalize_xtc_arrays(coords, box, times, steps, natoms: int):
    coords, box = np.asarray(coords), np.asarray(box)
    times, steps = np.asarray(times), np.asarray(steps)
    nframes = len(times)
    if coords.shape == (natoms, 3, nframes):
        coords = np.moveaxis(coords, 2, 0)
    elif coords.shape != (nframes, natoms, 3):
        raise RuntimeError(f"Unexpected XTC coordinate shape {coords.shape}")
    if box.shape == (3, 3, nframes):
        box = np.moveaxis(box, 2, 0)
    elif box.shape != (nframes, 3, 3):
        raise RuntimeError(f"Unexpected XTC box shape {box.shape}")
    return coords, box, times, steps


def load_sampled_worker(worker_dir: Path, natoms: int, stride: int):
    frames, boxes, times, steps, sources = [], [], [], [], []
    for segment in worker_xtc_segments(worker_dir):
        segment_natoms = int(get_xtc_natoms(str(segment).encode("utf-8")))
        if segment_natoms != natoms:
            raise RuntimeError(f"{segment} has {segment_natoms} atoms; expected {natoms}")
        raw = read_xtc(str(segment).encode("utf-8"))
        coords, box, seg_times, seg_steps = normalize_xtc_arrays(*raw, natoms)
        selection = np.arange(0, len(seg_times), stride, dtype=int)
        frames.append(np.asarray(coords[selection], dtype=np.float32))
        boxes.append(np.asarray(box[selection], dtype=np.float32))
        times.append(np.asarray(seg_times[selection], dtype=float))
        steps.append(np.asarray(seg_steps[selection], dtype=np.int64))
        sources.extend([segment.name] * len(selection))
        del raw, coords, box
    if not frames:
        raise FileNotFoundError(f"No XTC segments found in {worker_dir}")
    xyz = np.concatenate(frames, axis=0)
    unitcell = np.concatenate(boxes, axis=0)
    time_ps = np.concatenate(times)
    step = np.concatenate(steps)
    source_array = np.asarray(sources, dtype=object)
    order = np.lexsort((time_ps, step))
    xyz, unitcell = xyz[order], unitcell[order]
    time_ps, step, source_array = time_ps[order], step[order], source_array[order]
    keep = np.ones(len(step), dtype=bool)
    if len(step) > 1:
        keep[1:] = step[1:] != step[:-1]
    return xyz[keep], unitcell[keep], time_ps[keep], step[keep], source_array[keep]


def minimum_image(delta: np.ndarray, box: np.ndarray, inverse_box: Optional[np.ndarray] = None) -> np.ndarray:
    box = np.asarray(box, dtype=float)
    inverse_box = np.linalg.inv(box) if inverse_box is None else inverse_box
    fractional = np.asarray(delta, dtype=float) @ inverse_box
    fractional -= np.round(fractional)
    return fractional @ box


def make_mmp_whole(xyz: np.ndarray, box: np.ndarray, topology: SystemTopology) -> np.ndarray:
    inverse = np.linalg.inv(np.asarray(box, dtype=float))
    whole = np.asarray(xyz, dtype=float).copy()
    for parent, child in topology.unwrap_parent_order:
        whole[child] = whole[parent] + minimum_image(xyz[child] - xyz[parent], box, inverse)
    return whole


def weighted_com(xyz: np.ndarray, indices: np.ndarray, masses: np.ndarray) -> np.ndarray:
    return np.average(xyz[indices], axis=0, weights=masses)


def measure_frame(xyz: np.ndarray, box: np.ndarray, topology: SystemTopology) -> dict[str, float]:
    if xyz.shape != (topology.total_atoms, 3):
        raise RuntimeError(f"Coordinate shape {xyz.shape} does not match {topology.total_atoms} atoms")
    whole = make_mmp_whole(xyz, box, topology)
    cat_com = weighted_com(whole, topology.catalytic_indices, topology.catalytic_masses)
    hpx_com = weighted_com(whole, topology.hpx_indices, topology.hpx_masses)
    ser_vector = whole[topology.landmark_indices["ser347"]] - whole[topology.landmark_indices["ser123"]]
    inverse = np.linalg.inv(np.asarray(box, dtype=float))
    opening_vector = minimum_image(
        xyz[topology.landmark_indices["t230"]] - xyz[topology.landmark_indices["n171"]],
        box,
        inverse,
    )
    result = {
        "catalytic_hpx_com_nm": float(np.linalg.norm(hpx_com - cat_com)),
        "ser123_ser347_ca_nm": float(np.linalg.norm(ser_vector)),
        "catalytic_opening_n171_t230_ca_nm": float(np.linalg.norm(opening_vector)),
    }
    if not all(np.isfinite(v) and v > 0 for v in result.values()):
        raise RuntimeError(f"Invalid distance measurement: {result}")
    for name in METRICS_NM:
        result[name.replace("_nm", "_angstrom")] = 10.0 * result[name]
    return result


def new_metric_stats() -> dict[str, RunningStats]:
    return {metric: RunningStats() for metric in METRICS_NM}


def update_metric_stats(stats: dict[str, RunningStats], measurements: dict[str, float]) -> None:
    for metric in METRICS_NM:
        stats[metric].add(measurements[metric])


def summary_row(prefix: dict, stats: dict[str, RunningStats]) -> dict:
    row = dict(prefix)
    for metric in METRICS_NM:
        row.update(stats[metric].as_fields(metric))
        # Add Angstrom means/SDs for convenient plotting and reporting.
        row[f"{metric.replace('_nm', '_angstrom')}_mean"] = 10.0 * stats[metric].mean
        row[f"{metric.replace('_nm', '_angstrom')}_sd"] = 10.0 * stats[metric].sd
    return row


def write_csv(path: Path, rows: list[dict], fieldnames: Optional[list[str]] = None) -> None:
    if not rows and fieldnames is None:
        raise ValueError(f"Cannot infer fields for empty CSV: {path}")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames or list(rows[0]), extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def comparison_rows(generation_rows: list[dict]) -> list[dict]:
    lookup = {(r["system"], int(r["generation_index"])): r for r in generation_rows}
    rows = []
    generations = sorted({int(r["generation_index"]) for r in generation_rows})
    for generation in generations:
        wt = lookup.get(("wild_type", generation))
        if wt is None:
            continue
        for mutant in ("G978S", "G984C", "G987R"):
            current = lookup.get((mutant, generation))
            if current is None:
                continue
            for metric in METRICS_NM:
                wt_mean = float(wt[f"{metric}_mean"])
                mutant_mean = float(current[f"{metric}_mean"])
                rows.append({
                    "generation_index": generation,
                    "mutant": mutant,
                    "metric": metric,
                    "metric_label": METRIC_LABELS[metric],
                    "wild_type_mean_nm": wt_mean,
                    "mutant_mean_nm": mutant_mean,
                    "mutant_minus_wild_type_nm": mutant_mean - wt_mean,
                    "mutant_minus_wild_type_angstrom": 10.0 * (mutant_mean - wt_mean),
                    "percent_difference_from_wild_type": 100.0 * (mutant_mean - wt_mean) / wt_mean,
                    "wild_type_population_sd_nm": wt[f"{metric}_sd"],
                    "mutant_population_sd_nm": current[f"{metric}_sd"],
                    "wild_type_sampled_frames": wt[f"{metric}_n"],
                    "mutant_sampled_frames": current[f"{metric}_n"],
                })
    return rows


def make_plots(generation_rows: list[dict], output_dir: Path) -> Optional[str]:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:  # pragma: no cover - optional server dependency
        return f"Matplotlib unavailable; CSV outputs are complete but plots were skipped: {exc}"
    # Manuscript-wide system palette, established in Figure 3.
    colours = {"wild_type": "#303030", "G978S": "#2574A9", "G984C": "#2E8B57", "G987R": "#C0392B"}
    for metric in METRICS_NM:
        fig, ax = plt.subplots(figsize=(10.5, 6.2))
        for system in SYSTEMS:
            selected = sorted(
                (r for r in generation_rows if r["system"] == system),
                key=lambda r: int(r["generation_index"]),
            )
            if not selected:
                continue
            x = np.asarray([int(r["generation_index"]) for r in selected])
            y = np.asarray([10.0 * float(r[f"{metric}_mean"]) for r in selected])
            sd = np.asarray([10.0 * float(r[f"{metric}_sd"]) for r in selected])
            ax.plot(x, y, marker="o", markersize=3.5, linewidth=1.7, color=colours[system], label=SYSTEMS[system]["label"])
            ax.fill_between(x, y - sd, y + sd, color=colours[system], alpha=0.13, linewidth=0)
        ax.set_xlabel("Adaptive generation")
        ax.set_ylabel("Distance (Angstrom)")
        ax.set_title(METRIC_LABELS[metric])
        ax.grid(alpha=0.2)
        ax.legend(frameon=False, ncol=2)
        fig.tight_layout()
        fig.savefig(output_dir / f"{metric}_by_generation.png", dpi=300, facecolor="white")
        plt.close(fig)
    return None


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compare MMP-1 domain and catalytic-opening distances across wild type and three OI mutants.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--root", type=Path, default=Path("."), help="Adaptive-MD repository root")
    parser.add_argument("--output_dir", type=Path, default=Path("mmp1_domain_opening_analysis"))
    parser.add_argument("--generations", default="0-29", help="Comma-separated generations and inclusive ranges")
    parser.add_argument("--frame_stride", type=int, default=20, help="Read every Nth stored XTC frame")
    parser.add_argument("--wild_type_run", type=Path, help="Explicit wild-type adaptive_run_* directory")
    parser.add_argument("--g978s_run", type=Path, help="Explicit G978S adaptive_run_* directory")
    parser.add_argument("--g984c_run", type=Path, help="Explicit G984C adaptive_run_* directory")
    parser.add_argument("--g987r_run", type=Path, help="Explicit G987R adaptive_run_* directory")
    parser.add_argument("--validate_only", action="store_true", help="Validate topology selections without reading XTC files")
    args = parser.parse_args()
    if args.frame_stride < 1:
        raise SystemExit("--frame_stride must be at least 1")

    root = args.root.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()
    generations = parse_int_list(args.generations)
    explicit_runs = {
        "wild_type": args.wild_type_run,
        "G978S": args.g978s_run,
        "G984C": args.g984c_run,
        "G987R": args.g987r_run,
    }

    topologies = {system: build_system_topology(root, system) for system in SYSTEMS}
    for system, topology in topologies.items():
        print(
            f"VALID {system}: {topology.total_atoms} total atoms; "
            f"CAT={len(topology.catalytic_indices)} atoms; HPX={len(topology.hpx_indices)} atoms; "
            f"landmarks={topology.landmark_indices}",
            flush=True,
        )
    if args.validate_only:
        print("Topology, domain and residue-identity validation passed for all four systems")
        return
    if read_xtc is None or get_xtc_natoms is None or openmm is None:
        raise SystemExit(
            "OpenMM is required for XTC analysis. Activate the simulation environment and rerun. "
            f"Import error: {OPENMM_IMPORT_ERROR}"
        )

    output_dir.mkdir(parents=True, exist_ok=True)
    script_path = Path(__file__).resolve()
    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "repository_root": str(root),
        "repository_git_commit": git_commit(root),
        "script": str(script_path),
        "script_sha256": sha256(script_path),
        "generations": generations,
        "frame_stride": args.frame_stride,
        "definitions": {
            "catalytic_domain_coordinate_residues": list(CATALYTIC_RANGE),
            "hpx_domain_coordinate_residues": list(HPX_RANGE),
            "domain_com": "all selected atoms, actual ITP masses, covalently unwrapped MMP-1",
            "serine_landmarks": "SER123 CA to SER347 CA in covalently unwrapped MMP-1",
            "catalytic_opening": "minimum-image CA distance, biological N171/T230 = coordinate ASN152/THR211",
            "sd": "population SD across sampled frames; not a confidence interval",
        },
        "systems": {},
        "warnings": [],
    }

    time_fields = [
        "system", "system_label", "generation_index", "worker", "sampled_frame_index",
        "step", "time_ps", "trajectory_segment",
        "catalytic_hpx_com_nm", "catalytic_hpx_com_angstrom",
        "ser123_ser347_ca_nm", "ser123_ser347_ca_angstrom",
        "catalytic_opening_n171_t230_ca_nm", "catalytic_opening_n171_t230_ca_angstrom",
    ]
    worker_rows: list[dict] = []
    generation_stats: dict[tuple[str, int], dict[str, RunningStats]] = defaultdict(new_metric_stats)
    system_stats: dict[str, dict[str, RunningStats]] = defaultdict(new_metric_stats)
    time_path = output_dir / "mmp1_domain_opening_timeseries.csv"

    with time_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=time_fields)
        writer.writeheader()
        for system, topology in topologies.items():
            explicit = explicit_runs[system]
            if explicit is not None:
                run_dir = explicit.expanduser().resolve()
                if not run_dir.is_dir():
                    raise FileNotFoundError(f"Explicit run directory does not exist: {run_dir}")
            else:
                run_dir = find_latest_run(topology.input_dir, SYSTEMS[system]["run_prefix"], generations)
            manifest["systems"][system] = {
                "input_dir": str(topology.input_dir),
                "input_gro": str(topology.input_gro),
                "run_dir": str(run_dir),
                "total_atoms": topology.total_atoms,
                "landmark_indices_zero_based": topology.landmark_indices,
                "workers": [],
            }
            for generation in generations:
                generation_dir = run_dir / f"generation_{generation:03d}"
                if not generation_dir.is_dir():
                    raise FileNotFoundError(generation_dir)
                workers = sorted(path for path in generation_dir.glob("worker_*") if path.is_dir())
                if not workers:
                    raise FileNotFoundError(f"No worker directories in {generation_dir}")
                for worker_dir in workers:
                    print(f"[{datetime.now().isoformat(timespec='seconds')}] {system} generation {generation} {worker_dir.name}", flush=True)
                    xyz, boxes, times_ps, steps, sources = load_sampled_worker(
                        worker_dir, topology.total_atoms, args.frame_stride
                    )
                    stats = new_metric_stats()
                    for frame_index, (frame, box, time_ps, step, source) in enumerate(
                        zip(xyz, boxes, times_ps, steps, sources)
                    ):
                        measurements = measure_frame(frame, box, topology)
                        update_metric_stats(stats, measurements)
                        update_metric_stats(generation_stats[(system, generation)], measurements)
                        update_metric_stats(system_stats[system], measurements)
                        writer.writerow({
                            "system": system,
                            "system_label": SYSTEMS[system]["label"],
                            "generation_index": generation,
                            "worker": worker_dir.name,
                            "sampled_frame_index": frame_index,
                            "step": int(step),
                            "time_ps": float(time_ps),
                            "trajectory_segment": str(source),
                            **measurements,
                        })
                    row = summary_row({
                        "system": system,
                        "system_label": SYSTEMS[system]["label"],
                        "generation_index": generation,
                        "worker": worker_dir.name,
                    }, stats)
                    worker_rows.append(row)
                    manifest["systems"][system]["workers"].append({
                        "generation": generation,
                        "worker": worker_dir.name,
                        "sampled_frames": len(xyz),
                        "segments": [str(p) for p in worker_xtc_segments(worker_dir)],
                    })
                    handle.flush()
                    del xyz, boxes

    generation_rows = [
        summary_row({
            "system": system,
            "system_label": SYSTEMS[system]["label"],
            "generation_index": generation,
        }, stats)
        for (system, generation), stats in sorted(
            generation_stats.items(), key=lambda item: (list(SYSTEMS).index(item[0][0]), item[0][1])
        )
    ]
    overall_rows = [
        summary_row({"system": system, "system_label": SYSTEMS[system]["label"]}, system_stats[system])
        for system in SYSTEMS
    ]
    comparisons = comparison_rows(generation_rows)
    write_csv(output_dir / "worker_summary.csv", worker_rows)
    write_csv(output_dir / "system_generation_summary.csv", generation_rows)
    write_csv(output_dir / "overall_system_summary.csv", overall_rows)
    write_csv(output_dir / "mutant_vs_wild_type_by_generation.csv", comparisons)
    plot_warning = make_plots(generation_rows, output_dir)
    if plot_warning:
        manifest["warnings"].append(plot_warning)
        print(f"WARNING: {plot_warning}", file=sys.stderr)
    manifest["outputs"] = sorted(str(p.relative_to(output_dir)) for p in output_dir.iterdir() if p.is_file())
    (output_dir / "analysis_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    archive = shutil.make_archive(str(output_dir), "zip", root_dir=output_dir)
    print(f"Analysis complete: {output_dir}")
    print(f"Archive: {archive}")


if __name__ == "__main__":
    main()
