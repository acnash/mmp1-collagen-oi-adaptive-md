#!/usr/bin/env python3
"""Trajectory-wide peer-review analysis for the MMP-1/collagen campaigns.

Upload this single file to the repository root on SOL and run it from there.
It reads existing production XTC files; it does not run molecular dynamics.

The default analysis samples every 20th stored frame (200 ps because the
campaign wrote frames every 10 ps) from generations 0, 10, 20 and 28 for all
four systems.  It writes compact CSV tables and a ZIP archive suitable for
download through the web portal.  It also exports the generation-28 worker
maxima and the selected wild-type and mutant structures as PyMOL-ready PDBs.

Only NumPy and OpenMM are required.  The XTC reader bundled with OpenMM is
used directly, so MDTraj and MDAnalysis are not required.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import shutil
import subprocess
import sys
import time
from collections import defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Optional

import numpy as np

try:
    import openmm
    from openmm.app.internal.xtc_utils import get_xtc_natoms, read_xtc
except Exception as exc:  # pragma: no cover - server dependency check
    raise SystemExit(
        "OpenMM is required. Activate the same openmm environment used for the "
        f"simulations and rerun this file. Import error: {exc}"
    )


SYSTEMS = {
    "wild_type": {
        "label": "Wild type",
        "relative_dir": Path("systems/mmp1-bound-adaptive/wild_type"),
        "run_prefix": "adaptive_run_wild_type_",
    },
    "G978S": {
        "label": "G978S",
        "relative_dir": Path("systems/mmp1-bound-adaptive/collagen_G978S"),
        "run_prefix": "adaptive_run_G978S_",
    },
    "G984C": {
        "label": "G984C",
        "relative_dir": Path("systems/mmp1-bound-adaptive/collagen_G984C"),
        "run_prefix": "adaptive_run_G984C_",
    },
    "G987R": {
        "label": "G987R",
        "relative_dir": Path("systems/mmp1-bound-adaptive/collagen_G987R"),
        "run_prefix": "adaptive_run_G987R_",
    },
}

WINDOW_POSITIONS = tuple(range(977, 988))
CHAIN_NAMES = ("chain_1", "chain_2", "chain_3")
CHAIN_PAIRS = (("chain_1", "chain_2"), ("chain_1", "chain_3"), ("chain_2", "chain_3"))
HBOND_CRITERIA = ((0.32, 120.0), (0.35, 120.0), (0.35, 135.0), (0.35, 150.0))
CONTACT_CUTOFF_NM = 0.45
CONTACT_CANDIDATE_CUTOFF_NM = 1.20


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
    chain: str
    resid: int
    resname: str
    name: str
    atom_type: str
    mass: float

    @property
    def is_hydrogen(self) -> bool:
        return self.mass < 2.2 or self.name.upper().startswith("H")

    @property
    def is_heavy(self) -> bool:
        return not self.is_hydrogen


@dataclass
class SystemTopology:
    system: str
    input_dir: Path
    input_gro: Path
    atoms: list[AtomMeta]
    collagen_atoms: list[AtomMeta]
    mmp_atoms: list[AtomMeta]
    collagen_bonds: list[tuple[int, int]]
    protein_bonds: list[tuple[int, int]]
    total_atoms: int


def parse_int_list(text: str) -> list[int]:
    values = []
    for part in text.split(","):
        part = part.strip()
        if part:
            values.append(int(part))
    if not values:
        raise ValueError("At least one generation must be supplied")
    return values


def sha256(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            block = handle.read(chunk_size)
            if not block:
                break
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
    """Parse every moleculetype in a GROMACS ITP file."""
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
                atom = LocalAtom(
                    nr=int(fields[0]),
                    atom_type=fields[1],
                    resid=int(fields[2]),
                    resname=fields[3],
                    name=atom_name,
                    mass=float(fields[7]) if len(fields) >= 8 else inferred_mass,
                )
            except ValueError:
                continue
            raw[current]["atoms"].append(atom)
        elif section == "bonds" and len(fields) >= 2 and fields[0].isdigit() and fields[1].isdigit():
            raw[current]["bonds"].append((int(fields[0]), int(fields[1])))

    result = {}
    for name, values in raw.items():
        atoms = tuple(sorted(values["atoms"], key=lambda atom: atom.nr))
        if atoms:
            result[name] = MoleculeType(name=name, atoms=atoms, bonds=tuple(values["bonds"]))
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
    preferred = [path for path in candidates if "150mM_NaCl" in path.name and "openmm" not in path.name]
    if preferred:
        return preferred[0]
    if candidates:
        return candidates[0]
    raise FileNotFoundError(f"No input GRO file found in {input_dir}")


def build_system_topology(root: Path, system: str) -> SystemTopology:
    input_dir = (root / SYSTEMS[system]["relative_dir"]).resolve()
    system_top = input_dir / "system.top"
    if not system_top.exists():
        raise FileNotFoundError(system_top)

    molecule_types: dict[str, MoleculeType] = {}
    for itp in sorted(input_dir.glob("*.itp")):
        molecule_types.update(parse_itp_molecule_types(itp))

    atoms: list[AtomMeta] = []
    collagen_atoms: list[AtomMeta] = []
    mmp_atoms: list[AtomMeta] = []
    collagen_bonds: list[tuple[int, int]] = []
    protein_bonds: list[tuple[int, int]] = []
    offset = 0
    collagen_copy = 0
    molecule_copy_counter: dict[str, int] = defaultdict(int)

    for molecule_name, copies in parse_molecule_table(system_top):
        if molecule_name not in molecule_types:
            raise RuntimeError(f"No ITP atom definition found for molecule {molecule_name} in {input_dir}")
        molecule = molecule_types[molecule_name]
        if [atom.nr for atom in molecule.atoms] != list(range(1, len(molecule.atoms) + 1)):
            raise RuntimeError(f"Non-contiguous atom numbering in molecule {molecule_name}")

        for _ in range(copies):
            molecule_copy_counter[molecule_name] += 1
            molecule_copy = molecule_copy_counter[molecule_name]
            lower = molecule_name.lower()
            category = "other"
            chain = ""
            if lower.startswith("mmp1"):
                category = "mmp1"
                chain = "mmp1"
            elif lower.startswith("collagen"):
                category = "collagen"
                collagen_copy += 1
                chain = f"chain_{collagen_copy}"

            if category != "other":
                for local in molecule.atoms:
                    meta = AtomMeta(
                        index=offset + local.nr - 1,
                        molecule=molecule_name,
                        molecule_copy=molecule_copy,
                        category=category,
                        chain=chain,
                        resid=local.resid,
                        resname=local.resname,
                        name=local.name,
                        atom_type=local.atom_type,
                        mass=local.mass,
                    )
                    atoms.append(meta)
                    if category == "collagen":
                        collagen_atoms.append(meta)
                    else:
                        mmp_atoms.append(meta)
                if category == "collagen":
                    for first, second in molecule.bonds:
                        collagen_bonds.append((offset + first - 1, offset + second - 1))
                if category in {"mmp1", "collagen"}:
                    for first, second in molecule.bonds:
                        protein_bonds.append((offset + first - 1, offset + second - 1))
            offset += len(molecule.atoms)

    if collagen_copy != 3:
        raise RuntimeError(f"Expected three collagen molecules in {system_top}; found {collagen_copy}")
    if not mmp_atoms:
        raise RuntimeError(f"No MMP-1 molecule identified in {system_top}")

    input_gro = choose_input_gro(input_dir)
    gro_lines = input_gro.read_text(errors="ignore").splitlines()
    gro_count = int(gro_lines[1].strip())
    if gro_count != offset:
        raise RuntimeError(f"Topology has {offset} atoms but {input_gro} has {gro_count}")

    return SystemTopology(
        system=system,
        input_dir=input_dir,
        input_gro=input_gro,
        atoms=atoms,
        collagen_atoms=collagen_atoms,
        mmp_atoms=mmp_atoms,
        collagen_bonds=collagen_bonds,
        protein_bonds=protein_bonds,
        total_atoms=offset,
    )


def find_latest_run(input_dir: Path, prefix: str, required_generations: Optional[Iterable[int]] = None) -> Path:
    runs = sorted(path for path in input_dir.glob(f"{prefix}*") if path.is_dir())
    if not runs:
        raise FileNotFoundError(f"No adaptive run beginning {prefix!r} found in {input_dir}")
    required = list(required_generations or [])
    complete_candidates = [
        run
        for run in runs
        if all((run / f"generation_{generation:03d}").is_dir() for generation in required)
    ]
    if complete_candidates:
        return complete_candidates[-1].resolve()
    return runs[-1].resolve()


def resolve_recorded_path(path_like: str, run_dir: Path) -> Optional[Path]:
    path = Path(path_like)
    if path.exists():
        return path.resolve()
    text = str(path_like)
    marker = "adaptive_run_"
    position = text.find(marker)
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
    segments = []
    seen = set()
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
    coords = np.asarray(coords)
    box = np.asarray(box)
    times = np.asarray(times)
    steps = np.asarray(steps)
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
    frames = []
    boxes = []
    times = []
    steps = []
    sources = []
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

    order = np.lexsort((time_ps, step))
    xyz = xyz[order]
    unitcell = unitcell[order]
    time_ps = time_ps[order]
    step = step[order]
    source_array = np.asarray(sources, dtype=object)[order]

    # Resumed segments can contain the checkpoint frame twice.
    keep = np.ones(len(step), dtype=bool)
    if len(step) > 1:
        keep[1:] = step[1:] != step[:-1]
    return xyz[keep], unitcell[keep], time_ps[keep], step[keep], source_array[keep]


def minimum_image(delta: np.ndarray, box: np.ndarray, inverse_box: Optional[np.ndarray] = None) -> np.ndarray:
    if inverse_box is None:
        inverse_box = np.linalg.inv(np.asarray(box, dtype=float))
    fractional = np.asarray(delta, dtype=float) @ inverse_box
    fractional -= np.round(fractional)
    return fractional @ np.asarray(box, dtype=float)


def pbc_distance(first: np.ndarray, second: np.ndarray, box: np.ndarray, inverse_box: Optional[np.ndarray] = None) -> float:
    return float(np.linalg.norm(minimum_image(second - first, box, inverse_box)))


def pbc_distances(reference: np.ndarray, points: np.ndarray, box: np.ndarray, inverse_box: Optional[np.ndarray] = None) -> np.ndarray:
    delta = minimum_image(np.asarray(points) - np.asarray(reference), box, inverse_box)
    return np.linalg.norm(delta, axis=-1)


def angle_dha(donor: np.ndarray, hydrogen: np.ndarray, acceptor: np.ndarray, box: np.ndarray, inverse_box: Optional[np.ndarray] = None) -> float:
    hd = minimum_image(donor - hydrogen, box, inverse_box)
    ha = minimum_image(acceptor - hydrogen, box, inverse_box)
    denominator = np.linalg.norm(hd) * np.linalg.norm(ha)
    if denominator == 0:
        return float("nan")
    cosine = float(np.clip(np.dot(hd, ha) / denominator, -1.0, 1.0))
    return float(np.degrees(np.arccos(cosine)))


def pbc_dihedral(a: np.ndarray, b: np.ndarray, c: np.ndarray, d: np.ndarray, box: np.ndarray, inverse_box: Optional[np.ndarray] = None) -> float:
    p0 = np.asarray(a, dtype=float)
    p1 = p0 + minimum_image(b - a, box, inverse_box)
    p2 = p1 + minimum_image(c - b, box, inverse_box)
    p3 = p2 + minimum_image(d - c, box, inverse_box)
    b0 = -(p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2
    norm = np.linalg.norm(b1)
    if norm == 0:
        return float("nan")
    b1 = b1 / norm
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return float(np.degrees(np.arctan2(y, x)))


def atom_lookup(topology: SystemTopology) -> dict[tuple[str, int, str], AtomMeta]:
    return {(atom.chain, atom.resid, atom.name): atom for atom in topology.collagen_atoms}


def collagen_residue_atoms(topology: SystemTopology) -> dict[tuple[str, int], list[AtomMeta]]:
    groups: dict[tuple[str, int], list[AtomMeta]] = defaultdict(list)
    for atom in topology.collagen_atoms:
        if atom.resid in WINDOW_POSITIONS:
            groups[(atom.chain, atom.resid)].append(atom)
    missing = [key for key in ((chain, resid) for chain in CHAIN_NAMES for resid in WINDOW_POSITIONS) if key not in groups]
    if missing:
        raise RuntimeError(f"Could not identify the complete collagen window; missing {missing}")
    return groups


def acceptor_atom(atom: AtomMeta) -> bool:
    name = atom.name.upper()
    if name in {"O", "OXT", "OD1", "OD2", "OE1", "OE2", "OG", "OG1", "OH", "SG"}:
        return True
    return False


def hbond_candidates(topology: SystemTopology):
    meta_by_index = {atom.index: atom for atom in topology.collagen_atoms}
    bonded: dict[int, list[int]] = defaultdict(list)
    for first, second in topology.collagen_bonds:
        bonded[first].append(second)
        bonded[second].append(first)

    donors = []
    acceptors = []
    for atom in topology.collagen_atoms:
        if atom.resid not in WINDOW_POSITIONS:
            continue
        if acceptor_atom(atom):
            acceptors.append(atom)
        if atom.is_heavy and atom.name.upper().startswith(("N", "O", "S")):
            for neighbor in bonded.get(atom.index, []):
                hydrogen = meta_by_index.get(neighbor)
                if hydrogen is not None and hydrogen.is_hydrogen:
                    donors.append((atom, hydrogen))

    candidates = []
    for donor, hydrogen in donors:
        for acceptor in acceptors:
            if donor.chain == acceptor.chain:
                continue
            if donor.index == acceptor.index:
                continue
            backbone = donor.name.upper() == "N" and acceptor.name.upper() in {"O", "OXT"}
            includes_sulfur = donor.name.upper().startswith("S") or acceptor.name.upper().startswith("S")
            candidates.append((donor, hydrogen, acceptor, backbone, includes_sulfur))
    return candidates


def dihedral_definitions(topology: SystemTopology):
    lookup = atom_lookup(topology)
    definitions = []
    for chain in CHAIN_NAMES:
        for resid in WINDOW_POSITIONS:
            if resid > min(WINDOW_POSITIONS):
                keys = ((chain, resid - 1, "C"), (chain, resid, "N"), (chain, resid, "CA"), (chain, resid, "C"))
                if all(key in lookup for key in keys):
                    definitions.append((chain, resid, "phi", tuple(lookup[key].index for key in keys)))
            if resid < max(WINDOW_POSITIONS):
                keys = ((chain, resid, "N"), (chain, resid, "CA"), (chain, resid, "C"), (chain, resid + 1, "N"))
                if all(key in lookup for key in keys):
                    definitions.append((chain, resid, "psi", tuple(lookup[key].index for key in keys)))
                keys = ((chain, resid, "CA"), (chain, resid, "C"), (chain, resid + 1, "N"), (chain, resid + 1, "CA"))
                if all(key in lookup for key in keys):
                    definitions.append((chain, resid, "omega", tuple(lookup[key].index for key in keys)))
    return definitions


def mmp_domain(resid: int) -> str:
    if resid <= 260:
        return "catalytic"
    if resid >= 275:
        return "hemopexin"
    return "linker"


def identify_active_site(topology: SystemTopology, reference_xyz: np.ndarray, box: np.ndarray):
    glutamate = [atom for atom in topology.mmp_atoms if atom.resid == 200 and atom.resname.upper() == "GLU"]
    glu_oxygen = [atom for atom in glutamate if atom.name.upper() in {"OE1", "OE2"}]
    if len(glu_oxygen) != 2:
        raise RuntimeError("Could not identify both Glu200 carboxylate oxygens")

    zinc = [atom for atom in topology.mmp_atoms if atom.name.upper() == "ZN"]
    if not zinc:
        raise RuntimeError("Could not identify MMP-1 zinc atoms")
    catalytic_zinc = min(
        zinc,
        key=lambda atom: min(pbc_distance(reference_xyz[atom.index], reference_xyz[oxygen.index], box) for oxygen in glu_oxygen),
    )

    histidines = {}
    for resid in (199, 203, 209):
        matches = [atom for atom in topology.mmp_atoms if atom.resid == resid and atom.name.upper() == "NE2"]
        if matches:
            histidines[resid] = matches[0]
    return catalytic_zinc, glu_oxygen, histidines


def choose_contact_candidates(
    xyz: np.ndarray,
    boxes: np.ndarray,
    window_heavy: list[AtomMeta],
    mmp_heavy: list[AtomMeta],
) -> list[AtomMeta]:
    query_indices = [atom.index for atom in window_heavy]
    mmp_indices = np.array([atom.index for atom in mmp_heavy], dtype=int)
    selected = set()
    probe_frames = sorted(set(np.linspace(0, len(xyz) - 1, min(11, len(xyz)), dtype=int).tolist()))
    for frame_index in probe_frames:
        frame = xyz[frame_index]
        box = boxes[frame_index]
        inverse_box = np.linalg.inv(np.asarray(box, dtype=float))
        mmp_xyz = frame[mmp_indices]
        minimum = np.full(len(mmp_indices), np.inf, dtype=float)
        for query in query_indices:
            minimum = np.minimum(minimum, pbc_distances(frame[query], mmp_xyz, box, inverse_box))
        for local_index in np.flatnonzero(minimum <= CONTACT_CANDIDATE_CUTOFF_NM):
            selected.add(int(local_index))
    return [mmp_heavy[index] for index in sorted(selected)]


def circular_summary(sum_sin: float, sum_cos: float, count: int):
    if count == 0:
        return float("nan"), float("nan")
    mean = math.atan2(sum_sin / count, sum_cos / count)
    radius = min(max(math.hypot(sum_sin / count, sum_cos / count), 1.0e-12), 1.0)
    std = math.sqrt(-2.0 * math.log(radius))
    return math.degrees(mean), math.degrees(std)


def opening_geometry(
    frame: np.ndarray,
    box: np.ndarray,
    lookup: dict[tuple[str, int, str], AtomMeta],
    inverse_box: Optional[np.ndarray] = None,
):
    pair_values: dict[tuple[str, str], list[float]] = {pair: [] for pair in CHAIN_PAIRS}
    areas = []
    asymmetries = []
    for resid in WINDOW_POSITIONS:
        points = {chain: frame[lookup[(chain, resid, "CA")].index] for chain in CHAIN_NAMES}
        d12 = pbc_distance(points["chain_1"], points["chain_2"], box, inverse_box)
        d13 = pbc_distance(points["chain_1"], points["chain_3"], box, inverse_box)
        d23 = pbc_distance(points["chain_2"], points["chain_3"], box, inverse_box)
        pair_values[("chain_1", "chain_2")].append(d12)
        pair_values[("chain_1", "chain_3")].append(d13)
        pair_values[("chain_2", "chain_3")].append(d23)
        v12 = minimum_image(points["chain_2"] - points["chain_1"], box, inverse_box)
        v13 = minimum_image(points["chain_3"] - points["chain_1"], box, inverse_box)
        areas.append(0.5 * np.linalg.norm(np.cross(v12, v13)))
        asymmetries.append(max(d12, d13, d23) - min(d12, d13, d23))
    pair_means = {pair: float(np.mean(values)) for pair, values in pair_values.items()}
    overall = float(np.mean([value for values in pair_values.values() for value in values]))
    return overall, pair_means, float(np.mean(areas)), float(np.mean(asymmetries))


def analyse_worker(
    system: str,
    generation: int,
    worker_dir: Path,
    topology: SystemTopology,
    stride: int,
    frame_writer,
    hbond_counts: dict,
    contact_counts: dict,
    dihedral_accumulator: dict,
):
    started = time.time()
    xyz, boxes, times_ps, steps, sources = load_sampled_worker(worker_dir, topology.total_atoms, stride)
    lookup = atom_lookup(topology)
    residue_groups = collagen_residue_atoms(topology)
    hbonds = hbond_candidates(topology)
    dihedrals = dihedral_definitions(topology)
    window_heavy = [atom for atom in topology.collagen_atoms if atom.resid in WINDOW_POSITIONS and atom.is_heavy]
    mmp_heavy = [atom for atom in topology.mmp_atoms if atom.is_heavy]
    candidate_mmp = choose_contact_candidates(xyz, boxes, window_heavy, mmp_heavy)
    reference_xyz = xyz[0]
    catalytic_zinc, glu_oxygen, histidines = identify_active_site(topology, reference_xyz, boxes[0])
    frame_rows = []

    hbond_donor_indices = np.asarray([item[0].index for item in hbonds], dtype=int)
    hbond_hydrogen_indices = np.asarray([item[1].index for item in hbonds], dtype=int)
    hbond_acceptor_indices = np.asarray([item[2].index for item in hbonds], dtype=int)
    candidate_mmp_indices = np.asarray([atom.index for atom in candidate_mmp], dtype=int)

    for frame_number, (frame, box, time_ps, step, source) in enumerate(zip(xyz, boxes, times_ps, steps, sources)):
        inverse_box = np.linalg.inv(np.asarray(box, dtype=float))
        opening, pair_means, triangle_area, asymmetry = opening_geometry(frame, box, lookup, inverse_box)

        hbond_metrics = {}
        if len(hbonds):
            donor_xyz = frame[hbond_donor_indices]
            hydrogen_xyz = frame[hbond_hydrogen_indices]
            acceptor_xyz = frame[hbond_acceptor_indices]
            da_vectors = minimum_image(acceptor_xyz - donor_xyz, box, inverse_box)
            da_distances = np.linalg.norm(da_vectors, axis=1)
            hd_vectors = minimum_image(donor_xyz - hydrogen_xyz, box, inverse_box)
            ha_vectors = minimum_image(acceptor_xyz - hydrogen_xyz, box, inverse_box)
            denominators = np.linalg.norm(hd_vectors, axis=1) * np.linalg.norm(ha_vectors, axis=1)
            cosines = np.divide(
                np.sum(hd_vectors * ha_vectors, axis=1),
                denominators,
                out=np.full(len(hbonds), np.nan, dtype=float),
                where=denominators > 0,
            )
            dha_angles = np.degrees(np.arccos(np.clip(cosines, -1.0, 1.0)))
        else:
            da_distances = np.empty(0, dtype=float)
            dha_angles = np.empty(0, dtype=float)
        for cutoff_nm, angle_cutoff in HBOND_CRITERIA:
            label = f"hbond_n_da_{cutoff_nm * 10:.1f}A_angle_{int(angle_cutoff)}"
            occupied = np.flatnonzero((da_distances <= cutoff_nm) & (dha_angles >= angle_cutoff))
            total = int(len(occupied))
            backbone_total = int(sum(bool(hbonds[index][3]) for index in occupied))
            for candidate_index in occupied:
                donor, hydrogen, acceptor, backbone, includes_sulfur = hbonds[int(candidate_index)]
                key = (
                    system,
                    generation,
                    worker_dir.name,
                    cutoff_nm,
                    angle_cutoff,
                    donor.chain,
                    donor.resid,
                    donor.resname,
                    donor.name,
                    hydrogen.name,
                    acceptor.chain,
                    acceptor.resid,
                    acceptor.resname,
                    acceptor.name,
                    backbone,
                    includes_sulfur,
                )
                hbond_counts[key] += 1
            hbond_metrics[label] = total
            hbond_metrics[f"{label}_backbone_only"] = backbone_total

        unique_contacts = set()
        domain_contacts = defaultdict(set)
        candidate_xyz = frame[candidate_mmp_indices] if len(candidate_mmp_indices) else np.empty((0, 3), dtype=float)
        for chain in CHAIN_NAMES:
            for resid in WINDOW_POSITIONS:
                collagen_atoms = [atom for atom in residue_groups[(chain, resid)] if atom.is_heavy]
                contacted_mmp_residues = set()
                if len(candidate_mmp_indices):
                    collagen_xyz = frame[[atom.index for atom in collagen_atoms]]
                    deltas = candidate_xyz[None, :, :] - collagen_xyz[:, None, :]
                    distances = np.linalg.norm(minimum_image(deltas, box, inverse_box), axis=2)
                    for candidate_index in np.flatnonzero(np.any(distances <= CONTACT_CUTOFF_NM, axis=0)):
                        mmp_atom = candidate_mmp[int(candidate_index)]
                        contacted_mmp_residues.add((mmp_atom.resid, mmp_atom.resname))
                for mmp_resid, mmp_resname in contacted_mmp_residues:
                    domain = mmp_domain(mmp_resid)
                    contact = (chain, resid, mmp_resid, mmp_resname, domain)
                    unique_contacts.add(contact)
                    domain_contacts[domain].add(contact)
                    contact_counts[(system, generation, worker_dir.name, *contact)] += 1

        for chain, resid, angle_name, indices in dihedrals:
            angle = pbc_dihedral(*(frame[index] for index in indices), box, inverse_box)
            if not np.isfinite(angle):
                continue
            key = (system, generation, worker_dir.name, chain, resid, angle_name)
            accumulator = dihedral_accumulator[key]
            radians = math.radians(angle)
            accumulator["sum_sin"] += math.sin(radians)
            accumulator["sum_cos"] += math.cos(radians)
            accumulator["count"] += 1
            if angle_name == "omega" and abs(angle) < 90.0:
                accumulator["cis_count"] += 1

        zinc_xyz = frame[catalytic_zinc.index]
        glu_zinc = min(pbc_distance(zinc_xyz, frame[atom.index], box, inverse_box) for atom in glu_oxygen)
        active_site = {
            "catalytic_zinc_resid": catalytic_zinc.resid,
            "glu200_oxygen_to_zinc_angstrom": glu_zinc * 10.0,
        }
        for histidine_resid, atom in histidines.items():
            active_site[f"his{histidine_resid}_ne2_to_zinc_angstrom"] = pbc_distance(zinc_xyz, frame[atom.index], box, inverse_box) * 10.0
        for chain in CHAIN_NAMES:
            carbon = lookup[(chain, 975, "C")]
            oxygen = lookup[(chain, 975, "O")]
            nitrogen = lookup[(chain, 976, "N")]
            active_site[f"{chain}_gly975_scissile_C_to_zinc_angstrom"] = pbc_distance(frame[carbon.index], zinc_xyz, box, inverse_box) * 10.0
            active_site[f"{chain}_gly975_scissile_O_to_zinc_angstrom"] = pbc_distance(frame[oxygen.index], zinc_xyz, box, inverse_box) * 10.0
            active_site[f"{chain}_leu976_scissile_N_to_zinc_angstrom"] = pbc_distance(frame[nitrogen.index], zinc_xyz, box, inverse_box) * 10.0
            active_site[f"{chain}_glu200_O_to_scissile_C_angstrom"] = min(
                pbc_distance(frame[atom.index], frame[carbon.index], box, inverse_box) for atom in glu_oxygen
            ) * 10.0

        row = {
            "system": system,
            "system_label": SYSTEMS[system]["label"],
            "generation_index": generation,
            "worker": worker_dir.name,
            "sampled_frame_index": frame_number,
            "step": int(step),
            "time_ps": float(time_ps),
            "trajectory_segment": str(source),
            "opening_score_angstrom": opening * 10.0,
            "chain_1_2_mean_ca_distance_angstrom": pair_means[("chain_1", "chain_2")] * 10.0,
            "chain_1_3_mean_ca_distance_angstrom": pair_means[("chain_1", "chain_3")] * 10.0,
            "chain_2_3_mean_ca_distance_angstrom": pair_means[("chain_2", "chain_3")] * 10.0,
            "mean_ca_triangle_area_angstrom2": triangle_area * 100.0,
            "mean_ca_distance_asymmetry_angstrom": asymmetry * 10.0,
            "unique_collagen_mmp_residue_contacts": len(unique_contacts),
            "unique_catalytic_domain_contacts": len(domain_contacts["catalytic"]),
            "unique_linker_contacts": len(domain_contacts["linker"]),
            "unique_hemopexin_domain_contacts": len(domain_contacts["hemopexin"]),
            **hbond_metrics,
            **active_site,
        }
        frame_rows.append(row)
        frame_writer.writerow(row)

    metrics = [
        key
        for key in frame_rows[0]
        if key not in {"system", "system_label", "generation_index", "worker", "sampled_frame_index", "step", "time_ps", "trajectory_segment", "catalytic_zinc_resid"}
    ]
    summary = {
        "system": system,
        "system_label": SYSTEMS[system]["label"],
        "generation_index": generation,
        "worker": worker_dir.name,
        "n_sampled_frames": len(frame_rows),
        "n_contact_candidate_mmp_atoms": len(candidate_mmp),
        "elapsed_seconds": time.time() - started,
    }
    for metric in metrics:
        values = np.asarray([row[metric] for row in frame_rows], dtype=float)
        summary[f"{metric}_mean"] = float(np.nanmean(values))
        summary[f"{metric}_median"] = float(np.nanmedian(values))
        summary[f"{metric}_std"] = float(np.nanstd(values, ddof=1)) if len(values) > 1 else 0.0
        summary[f"{metric}_q25"] = float(np.nanquantile(values, 0.25))
        summary[f"{metric}_q75"] = float(np.nanquantile(values, 0.75))
    return summary, len(frame_rows)


def write_csv(path: Path, rows: Iterable[dict], fieldnames: Optional[list[str]] = None):
    rows = list(rows)
    if not rows and fieldnames is None:
        path.write_text("", encoding="utf-8")
        return
    if fieldnames is None:
        fieldnames = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def output_inventory(output_dir: Path) -> list[dict]:
    rows = []
    for path in sorted(output_dir.rglob("*")):
        if path.is_file():
            rows.append({
                "path": str(path.relative_to(output_dir)),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            })
    return rows


def worker_maximum_opening(worker_dir: Path) -> float:
    """Return the largest recorded opening score for one worker."""
    score_path = worker_dir / "opening_scores.csv"
    if not score_path.exists():
        return float("nan")
    maximum = float("nan")
    with score_path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            try:
                value = float(row["opening_score_angstrom"])
            except (KeyError, TypeError, ValueError):
                continue
            if math.isfinite(value) and (not math.isfinite(maximum) or value > maximum):
                maximum = value
    return maximum


def read_gro_coordinates_angstrom(path: Path, expected_atoms: int) -> np.ndarray:
    """Read coordinates from a GRO file and return them in angstrom."""
    lines = path.read_text(errors="ignore").splitlines()
    if len(lines) < expected_atoms + 3:
        raise RuntimeError(f"Incomplete GRO file: {path}")
    atom_count = int(lines[1].strip())
    if atom_count != expected_atoms:
        raise RuntimeError(f"{path} contains {atom_count} atoms; expected {expected_atoms}")
    coordinates = np.empty((atom_count, 3), dtype=float)
    for index, line in enumerate(lines[2 : 2 + atom_count]):
        coordinates[index] = [float(line[20:28]), float(line[28:36]), float(line[36:44])]
    return coordinates * 10.0


def pdb_element(atom_name: str, resname: str) -> str:
    upper_name = atom_name.strip().upper()
    upper_resname = resname.strip().upper()
    if upper_resname == "ZN" or upper_name == "ZN":
        return "ZN"
    letters = "".join(character for character in upper_name if character.isalpha())
    return letters[0] if letters else "C"


def export_protein_pdb(gro_path: Path, topology: SystemTopology, output_path: Path) -> None:
    """Write MMP-1 and the three collagen chains from a selected GRO frame."""
    xyz = read_gro_coordinates_angstrom(gro_path, topology.total_atoms)
    chain_ids = {"mmp1": "A", "chain_1": "B", "chain_2": "C", "chain_3": "D"}
    lines = [f"REMARK Selected adaptive structure exported from {gro_path}"]
    serial = 1
    previous_chain = None
    for atom in sorted(topology.atoms, key=lambda item: item.index):
        chain_id = chain_ids[atom.chain]
        if previous_chain is not None and chain_id != previous_chain:
            lines.append("TER")
        x, y, z = xyz[atom.index]
        record = "HETATM" if atom.resname.upper() in {"ZN", "CA", "MG"} else "ATOM  "
        atom_name = atom.name[:4].rjust(4)
        resname = atom.resname[:3].rjust(3)
        element = pdb_element(atom.name, atom.resname).rjust(2)
        lines.append(
            f"{record}{serial:5d} {atom_name} {resname} {chain_id}{atom.resid:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{0.00:6.2f}          {element}"
        )
        serial += 1
        previous_chain = chain_id
    lines.extend(["TER", "END"])
    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def export_endpoint_structures(root: Path, output_dir: Path, generation: int) -> list[dict]:
    """Export worker maxima and the selected structure for the endpoint generation."""
    structure_dir = output_dir / "selected_generation_028_structures"
    structure_dir.mkdir(parents=True, exist_ok=True)
    maxima_rows = []
    selected_rows = []
    for system in SYSTEMS:
        topology = build_system_topology(root, system)
        run_dir = find_latest_run(topology.input_dir, SYSTEMS[system]["run_prefix"], [generation])
        generation_dir = run_dir / f"generation_{generation:03d}"
        candidates = []
        for worker_dir in sorted(path for path in generation_dir.glob("worker_*") if path.is_dir()):
            maximum = worker_maximum_opening(worker_dir)
            best_gro = worker_dir / "best_opening_frame.gro"
            if not math.isfinite(maximum) or not best_gro.exists():
                continue
            row = {
                "system": system,
                "system_label": SYSTEMS[system]["label"],
                "generation_index": generation,
                "worker": worker_dir.name,
                "maximum_opening_score_angstrom": maximum,
                "best_gro": str(best_gro),
            }
            maxima_rows.append(row)
            candidates.append((maximum, worker_dir.name, best_gro))
        if not candidates:
            raise RuntimeError(f"No completed worker structures found in {generation_dir}")
        maximum, worker_name, best_gro = max(candidates, key=lambda item: item[0])
        pdb_path = structure_dir / f"{system}_generation_{generation:03d}_selected.pdb"
        gro_copy = structure_dir / f"{system}_generation_{generation:03d}_selected.gro"
        export_protein_pdb(best_gro, topology, pdb_path)
        shutil.copy2(best_gro, gro_copy)
        selected_rows.append({
            "system": system,
            "system_label": SYSTEMS[system]["label"],
            "generation_index": generation,
            "selected_worker": worker_name,
            "opening_score_angstrom": maximum,
            "source_gro": str(best_gro),
            "exported_pdb": str(pdb_path.relative_to(output_dir)),
            "exported_gro": str(gro_copy.relative_to(output_dir)),
        })
    write_csv(output_dir / "generation_028_worker_maxima.csv", maxima_rows)
    write_csv(output_dir / "generation_028_selected_structures.csv", selected_rows)
    return selected_rows


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Analyse existing adaptive MD trajectories for peer-review revisions.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--root", type=Path, default=Path("."), help="Repository root on SOL.")
    parser.add_argument("--output_dir", type=Path, default=Path("peer_review_trajectory_analysis"))
    parser.add_argument("--generations", default="0,10,20,28", help="Comma-separated completed generations.")
    parser.add_argument(
        "--frame_stride",
        type=int,
        default=20,
        help="Read every Nth stored trajectory frame. The stored interval is 10 ps, so 20 gives 200 ps.",
    )
    args = parser.parse_args()
    if args.frame_stride < 1:
        raise SystemExit("--frame_stride must be at least 1")

    root = args.root.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    generations = parse_int_list(args.generations)

    script_path = Path(__file__).resolve()
    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "repository_root": str(root),
        "repository_git_commit": git_commit(root),
        "script": str(script_path),
        "script_sha256": sha256(script_path),
        "python_version": sys.version,
        "openmm_version": getattr(openmm, "__version__", "unknown"),
        "numpy_version": np.__version__,
        "new_simulations_run": False,
        "input_scope": "existing production XTC trajectories",
        "generations": generations,
        "stored_frame_stride": args.frame_stride,
        "nominal_sample_interval_ps": args.frame_stride * 10.0,
        "hbond_criteria": [
            {"donor_acceptor_cutoff_angstrom": cutoff * 10.0, "minimum_DHA_angle_deg": angle}
            for cutoff, angle in HBOND_CRITERIA
        ],
        "contact_cutoff_angstrom": CONTACT_CUTOFF_NM * 10.0,
        "systems": {},
        "warnings": [],
    }

    frame_path = output_dir / "trajectory_frame_metrics.csv"
    worker_summaries = []
    hbond_counts: dict[tuple, int] = defaultdict(int)
    contact_counts: dict[tuple, int] = defaultdict(int)
    dihedral_accumulator: dict[tuple, dict] = defaultdict(lambda: {"sum_sin": 0.0, "sum_cos": 0.0, "count": 0, "cis_count": 0})
    worker_frame_counts: dict[tuple[str, int, str], int] = {}

    # Establish the complete frame CSV schema before streaming rows.
    hbond_fields = []
    for cutoff_nm, angle_cutoff in HBOND_CRITERIA:
        label = f"hbond_n_da_{cutoff_nm * 10:.1f}A_angle_{int(angle_cutoff)}"
        hbond_fields.extend([label, f"{label}_backbone_only"])
    active_fields = ["catalytic_zinc_resid", "glu200_oxygen_to_zinc_angstrom"]
    active_fields.extend([f"his{resid}_ne2_to_zinc_angstrom" for resid in (199, 203, 209)])
    for chain in CHAIN_NAMES:
        active_fields.extend([
            f"{chain}_gly975_scissile_C_to_zinc_angstrom",
            f"{chain}_gly975_scissile_O_to_zinc_angstrom",
            f"{chain}_leu976_scissile_N_to_zinc_angstrom",
            f"{chain}_glu200_O_to_scissile_C_angstrom",
        ])
    frame_fields = [
        "system", "system_label", "generation_index", "worker", "sampled_frame_index", "step", "time_ps",
        "trajectory_segment", "opening_score_angstrom", "chain_1_2_mean_ca_distance_angstrom",
        "chain_1_3_mean_ca_distance_angstrom", "chain_2_3_mean_ca_distance_angstrom",
        "mean_ca_triangle_area_angstrom2", "mean_ca_distance_asymmetry_angstrom",
        "unique_collagen_mmp_residue_contacts", "unique_catalytic_domain_contacts", "unique_linker_contacts",
        "unique_hemopexin_domain_contacts", *hbond_fields, *active_fields,
    ]

    with frame_path.open("w", newline="", encoding="utf-8") as frame_handle:
        frame_writer = csv.DictWriter(frame_handle, fieldnames=frame_fields, extrasaction="ignore")
        frame_writer.writeheader()

        for system in SYSTEMS:
            topology = build_system_topology(root, system)
            run_dir = find_latest_run(topology.input_dir, SYSTEMS[system]["run_prefix"], generations)
            manifest["systems"][system] = {
                "input_dir": str(topology.input_dir),
                "input_gro": str(topology.input_gro),
                "run_dir": str(run_dir),
                "total_atoms": topology.total_atoms,
                "workers": [],
            }
            for generation in generations:
                generation_dir = run_dir / f"generation_{generation:03d}"
                if not generation_dir.exists():
                    manifest["warnings"].append(f"Missing {generation_dir}")
                    continue
                worker_dirs = sorted(path for path in generation_dir.glob("worker_*") if path.is_dir())
                for worker_dir in worker_dirs:
                    segments = worker_xtc_segments(worker_dir)
                    if not segments:
                        manifest["warnings"].append(f"No XTC files in {worker_dir}")
                        continue
                    print(f"[{datetime.now().isoformat(timespec='seconds')}] {system} generation {generation} {worker_dir.name}", flush=True)
                    try:
                        summary, nframes = analyse_worker(
                            system,
                            generation,
                            worker_dir,
                            topology,
                            args.frame_stride,
                            frame_writer,
                            hbond_counts,
                            contact_counts,
                            dihedral_accumulator,
                        )
                    except Exception as exc:
                        message = f"FAILED {system} generation {generation} {worker_dir.name}: {type(exc).__name__}: {exc}"
                        print(message, file=sys.stderr, flush=True)
                        manifest["warnings"].append(message)
                        continue
                    worker_summaries.append(summary)
                    worker_frame_counts[(system, generation, worker_dir.name)] = nframes
                    manifest["systems"][system]["workers"].append({
                        "generation": generation,
                        "worker": worker_dir.name,
                        "segments": [str(path) for path in segments],
                        "sampled_frames": nframes,
                    })
                    frame_handle.flush()

    write_csv(output_dir / "worker_metric_summaries.csv", worker_summaries)

    hbond_rows = []
    for key, occupied_frames in sorted(hbond_counts.items()):
        (
            system, generation, worker, cutoff_nm, angle_cutoff, donor_chain, donor_resid, donor_resname,
            donor_atom, hydrogen_atom, acceptor_chain, acceptor_resid, acceptor_resname, acceptor_atom_name,
            backbone, includes_sulfur,
        ) = key
        total_frames = worker_frame_counts.get((system, generation, worker), 0)
        hbond_rows.append({
            "system": system,
            "system_label": SYSTEMS[system]["label"],
            "generation_index": generation,
            "worker": worker,
            "donor_acceptor_cutoff_angstrom": cutoff_nm * 10.0,
            "minimum_DHA_angle_deg": angle_cutoff,
            "donor_chain": donor_chain,
            "donor_resid": donor_resid,
            "donor_resname": donor_resname,
            "donor_atom": donor_atom,
            "hydrogen_atom": hydrogen_atom,
            "acceptor_chain": acceptor_chain,
            "acceptor_resid": acceptor_resid,
            "acceptor_resname": acceptor_resname,
            "acceptor_atom": acceptor_atom_name,
            "backbone_to_backbone": backbone,
            "includes_sulfur": includes_sulfur,
            "occupied_frames": occupied_frames,
            "total_sampled_frames": total_frames,
            "occupancy_fraction": occupied_frames / total_frames if total_frames else float("nan"),
        })
    write_csv(output_dir / "interchain_hbond_occupancies.csv", hbond_rows)

    contact_rows = []
    for key, occupied_frames in sorted(contact_counts.items()):
        system, generation, worker, chain, collagen_resid, mmp_resid, mmp_resname, domain = key
        total_frames = worker_frame_counts.get((system, generation, worker), 0)
        contact_rows.append({
            "system": system,
            "system_label": SYSTEMS[system]["label"],
            "generation_index": generation,
            "worker": worker,
            "collagen_chain": chain,
            "collagen_resid": collagen_resid,
            "mmp_resid": mmp_resid,
            "mmp_resname": mmp_resname,
            "mmp_domain": domain,
            "occupied_frames": occupied_frames,
            "total_sampled_frames": total_frames,
            "occupancy_fraction": occupied_frames / total_frames if total_frames else float("nan"),
        })
    write_csv(output_dir / "collagen_mmp_residue_contact_occupancies.csv", contact_rows)

    dihedral_rows = []
    for key, values in sorted(dihedral_accumulator.items()):
        system, generation, worker, chain, resid, angle_name = key
        mean_deg, std_deg = circular_summary(values["sum_sin"], values["sum_cos"], values["count"])
        dihedral_rows.append({
            "system": system,
            "system_label": SYSTEMS[system]["label"],
            "generation_index": generation,
            "worker": worker,
            "chain": chain,
            "resid": resid,
            "dihedral": angle_name,
            "n_sampled_frames": values["count"],
            "circular_mean_deg": mean_deg,
            "circular_std_deg": std_deg,
            "cis_fraction": values["cis_count"] / values["count"] if angle_name == "omega" and values["count"] else "",
        })
    write_csv(output_dir / "residue_chain_dihedral_summaries.csv", dihedral_rows)

    manifest["completed_utc"] = datetime.now(timezone.utc).isoformat()
    manifest["n_workers_analysed"] = len(worker_summaries)
    manifest["n_sampled_frames"] = sum(worker_frame_counts.values())
    manifest["selected_endpoint_structures"] = export_endpoint_structures(root, output_dir, max(generations))
    manifest_path = output_dir / "analysis_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    manifest["outputs"] = output_inventory(output_dir)
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    archive_base = output_dir.parent / output_dir.name
    archive_path = Path(shutil.make_archive(str(archive_base), "zip", root_dir=output_dir))
    print(f"Analysis complete: {output_dir}")
    print(f"Download this ZIP through the web portal: {archive_path}")


if __name__ == "__main__":
    main()
