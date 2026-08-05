#!/usr/bin/env python3
"""
Generate a joint structural report for adaptive MMP1-collagen unwinding runs.

The report uses the best worker/frame selected for each adaptive generation.
It analyzes collagen only:
  - hydrogen-bond distance/angle/count distributions
  - collagen backbone phi/psi/omega dihedral distributions
  - relationships between those metrics and the unwinding opening score

Inputs are adaptive_run_* directories produced by adaptive_mmp1_unwinding_dual_worker.py.
"""

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


SYSTEMS = {
    "wild_type": {
        "label": "Wild type",
        "run_prefix": "adaptive_run_wild_type_",
        "relative_dir": Path("generated_mutants/salted_150mM_NaCl/wild_type"),
        "color": "#333333",
    },
    "G978S": {
        "label": "G978S",
        "run_prefix": "adaptive_run_G978S_",
        "relative_dir": Path("generated_mutants/salted_150mM_NaCl/collagen_G978S"),
        "color": "#1f77b4",
    },
    "G984C": {
        "label": "G984C",
        "run_prefix": "adaptive_run_G984C_",
        "relative_dir": Path("generated_mutants/salted_150mM_NaCl/collagen_G984C"),
        "color": "#2ca02c",
    },
    "G987R": {
        "label": "G987R",
        "run_prefix": "adaptive_run_G987R_",
        "relative_dir": Path("generated_mutants/salted_150mM_NaCl/collagen_G987R"),
        "color": "#d62728",
    },
}

BACKBONE_ATOMS = {"N", "CA", "C"}
DONOR_HEAVY_PREFIXES = ("N", "O", "S")
ACCEPTOR_PREFIXES = ("O", "S")


@dataclass
class RunSpec:
    system: str
    label: str
    run_dir: Path


@dataclass
class Atom:
    index0: int
    atom_number: int
    resid: int
    resname: str
    name: str
    xyz_a: np.ndarray
    chain: Optional[str] = None


def read_json(path: Path) -> dict:
    try:
        return json.loads(path.read_text())
    except Exception:
        return {}


def safe_mkdir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def parse_generation_index(name: str) -> Optional[int]:
    try:
        return int(str(name).split("_")[-1])
    except Exception:
        return None


def find_latest_run(system_dir: Path, run_prefix: str) -> Optional[Path]:
    runs = sorted([p for p in system_dir.glob(f"{run_prefix}*") if p.is_dir()])
    if not runs:
        runs = sorted([p for p in system_dir.glob("adaptive_run_*") if p.is_dir()])
    return runs[-1].resolve() if runs else None


def resolve_run_specs(args) -> List[RunSpec]:
    root = Path(args.root).resolve()
    explicit = {
        "wild_type": args.wild_type_run,
        "G978S": args.g978s_run,
        "G984C": args.g984c_run,
        "G987R": args.g987r_run,
    }
    specs = []
    for system, meta in SYSTEMS.items():
        if explicit[system]:
            run_dir = Path(explicit[system]).expanduser().resolve()
        else:
            run_dir = find_latest_run((root / meta["relative_dir"]).resolve(), meta["run_prefix"])
            if run_dir is None:
                continue
        if not run_dir.exists():
            raise SystemExit(f"{system} run directory does not exist: {run_dir}")
        specs.append(RunSpec(system=system, label=meta["label"], run_dir=run_dir))
    return specs


def resolve_moved_path(path_like: Optional[str], run_dir: Path) -> Optional[Path]:
    if not path_like:
        return None
    path = Path(path_like)
    if path.exists():
        return path.resolve()
    marker = "adaptive_run_"
    text = str(path_like)
    idx = text.find(marker)
    if idx >= 0:
        rel = Path(*Path(text[idx:]).parts[1:])
        repaired = run_dir / rel
        if repaired.exists():
            return repaired.resolve()
    return None


def best_frames_from_context(spec: RunSpec) -> List[dict]:
    context = read_json(spec.run_dir / "context.json")
    rows = []
    for entry in context.get("selected_frames", []):
        generation_index = entry.get("generation")
        frames = entry.get("frames", [])
        if generation_index is None or not frames:
            continue
        frame = sorted(frames, key=lambda x: x.get("rank", 999))[0]
        best_gro = resolve_moved_path(frame.get("best_gro"), spec.run_dir)
        if best_gro is None:
            continue
        rows.append(
            {
                "system": spec.system,
                "system_label": spec.label,
                "run_dir": str(spec.run_dir),
                "generation_index": int(generation_index),
                "generation": f"generation_{int(generation_index):03d}",
                "worker_dir": frame.get("worker_dir", ""),
                "best_gro": best_gro,
                "opening_score_angstrom": frame.get("best_score_angstrom"),
                "best_step": frame.get("best_step"),
            }
        )
    return sorted(rows, key=lambda x: x["generation_index"])


def best_frames_from_generation_json(spec: RunSpec) -> List[dict]:
    rows = []
    for generation_dir in sorted(spec.run_dir.glob("generation_*")):
        if not generation_dir.is_dir():
            continue
        generation_index = parse_generation_index(generation_dir.name)
        meta = read_json(generation_dir / "generation.json")
        frames = meta.get("selected_top_frames", [])
        if generation_index is None or not frames:
            continue
        frame = sorted(frames, key=lambda x: x.get("rank", 999))[0]
        best_gro = resolve_moved_path(frame.get("best_gro"), spec.run_dir)
        if best_gro is None:
            continue
        rows.append(
            {
                "system": spec.system,
                "system_label": spec.label,
                "run_dir": str(spec.run_dir),
                "generation_index": int(generation_index),
                "generation": generation_dir.name,
                "worker_dir": frame.get("worker_dir", ""),
                "best_gro": best_gro,
                "opening_score_angstrom": frame.get("best_score_angstrom", meta.get("best_score_angstrom")),
                "best_step": frame.get("best_step"),
            }
        )
    return sorted(rows, key=lambda x: x["generation_index"])


def collect_best_frame_records(specs: List[RunSpec]) -> pd.DataFrame:
    rows = []
    for spec in specs:
        frames = best_frames_from_context(spec)
        if not frames:
            frames = best_frames_from_generation_json(spec)
        rows.extend(frames)
    df = pd.DataFrame(rows)
    if not df.empty:
        df["best_gro"] = df["best_gro"].astype(str)
        df = df.sort_values(["system", "generation_index"])
    return df


def get_input_dir(run_dir: Path) -> Path:
    context = read_json(run_dir / "context.json")
    input_dir = context.get("input_dir")
    if input_dir and Path(input_dir).exists():
        return Path(input_dir).resolve()
    return run_dir.parent.resolve()


def atom_count_from_itp(path: Path) -> int:
    in_atoms = False
    count = 0
    for line in path.read_text(errors="ignore").splitlines():
        s = line.strip()
        if not s or s.startswith(";"):
            continue
        if s.startswith("["):
            in_atoms = s.lower().startswith("[ atoms")
            continue
        if in_atoms:
            parts = s.split()
            if parts and parts[0].isdigit():
                count += 1
    return count


def molecule_lines(system_top: Path) -> List[Tuple[str, int]]:
    lines = []
    in_molecules = False
    for line in system_top.read_text(errors="ignore").splitlines():
        s = line.strip()
        if not s or s.startswith(";"):
            continue
        if s.startswith("["):
            in_molecules = s.lower().startswith("[ molecules")
            continue
        if in_molecules:
            parts = s.split()
            if len(parts) >= 2 and parts[1].isdigit():
                lines.append((parts[0], int(parts[1])))
    return lines


def collagen_chain_ranges(input_dir: Path) -> Dict[str, Tuple[int, int, str]]:
    counts = {}
    for itp in input_dir.glob("*.itp"):
        if itp.name in {"forcefield.itp", "ffbonded.itp", "ffnonbonded.itp", "tip3p.itp", "ions.itp"}:
            continue
        counts[itp.stem] = atom_count_from_itp(itp)
    ranges = {}
    current = 1
    chain_i = 0
    for mol_name, mol_count in molecule_lines(input_dir / "system.top"):
        atom_count = counts.get(mol_name)
        if atom_count is None:
            continue
        for _ in range(mol_count):
            start = current
            end = current + atom_count - 1
            if mol_name.startswith("collagen"):
                chain_i += 1
                ranges[f"chain_{chain_i}"] = (start, end, mol_name)
            current = end + 1
    if len(ranges) != 3:
        raise RuntimeError(f"Expected 3 collagen chains from {input_dir / 'system.top'}, found {len(ranges)}")
    return ranges


def parse_gro(path: Path) -> List[Atom]:
    lines = path.read_text(errors="ignore").splitlines()
    atoms = []
    for i, line in enumerate(lines[2:-1]):
        try:
            resid = int(line[0:5])
            resname = line[5:10].strip()
            name = line[10:15].strip()
            atom_number = int(line[15:20])
            xyz_nm = np.array([float(line[20:28]), float(line[28:36]), float(line[36:44])], dtype=float)
        except Exception:
            continue
        atoms.append(Atom(index0=i, atom_number=atom_number, resid=resid, resname=resname, name=name, xyz_a=xyz_nm * 10.0))
    return atoms


def assign_collagen_chains(atoms: List[Atom], ranges: Dict[str, Tuple[int, int, str]]) -> List[Atom]:
    out = []
    for atom in atoms:
        for chain, (start, end, _) in ranges.items():
            if start <= atom.atom_number <= end:
                atom.chain = chain
                out.append(atom)
                break
    return out


def group_residues(atoms: Iterable[Atom]) -> Dict[str, List[List[Atom]]]:
    chains: Dict[str, List[List[Atom]]] = {}
    for chain in ["chain_1", "chain_2", "chain_3"]:
        chain_atoms = [a for a in atoms if a.chain == chain]
        residues = []
        current_key = None
        current = []
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
        chains[chain] = residues
    return chains


def atom_by_name(residue_atoms: List[Atom], name: str) -> Optional[Atom]:
    for atom in residue_atoms:
        if atom.name == name:
            return atom
    return None


def dihedral_deg(p0, p1, p2, p3) -> float:
    b0 = -(p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2
    b1 = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return float(np.degrees(np.arctan2(y, x)))


def collect_dihedrals(frame: dict, residues_by_chain: Dict[str, List[List[Atom]]]) -> pd.DataFrame:
    rows = []
    for chain, residues in residues_by_chain.items():
        for i, residue in enumerate(residues):
            n = atom_by_name(residue, "N")
            ca = atom_by_name(residue, "CA")
            c = atom_by_name(residue, "C")
            if n is None or ca is None or c is None:
                continue
            base = {
                "system": frame["system"],
                "system_label": frame["system_label"],
                "generation_index": frame["generation_index"],
                "generation": frame["generation"],
                "opening_score_angstrom": frame.get("opening_score_angstrom"),
                "chain": chain,
                "resid": residue[0].resid,
                "resname": residue[0].resname,
                "best_gro": frame["best_gro"],
            }
            if i > 0:
                c_prev = atom_by_name(residues[i - 1], "C")
                if c_prev is not None:
                    rows.append({**base, "dihedral": "phi", "angle_deg": dihedral_deg(c_prev.xyz_a, n.xyz_a, ca.xyz_a, c.xyz_a)})
            if i < len(residues) - 1:
                n_next = atom_by_name(residues[i + 1], "N")
                ca_next = atom_by_name(residues[i + 1], "CA")
                if n_next is not None:
                    rows.append({**base, "dihedral": "psi", "angle_deg": dihedral_deg(n.xyz_a, ca.xyz_a, c.xyz_a, n_next.xyz_a)})
                if n_next is not None and ca_next is not None:
                    rows.append({**base, "dihedral": "omega", "angle_deg": dihedral_deg(ca.xyz_a, c.xyz_a, n_next.xyz_a, ca_next.xyz_a)})
    return pd.DataFrame(rows)


def is_heavy_donor(atom: Atom) -> bool:
    return atom.name and atom.name[0] in DONOR_HEAVY_PREFIXES and not atom.name.startswith("H")


def is_acceptor(atom: Atom) -> bool:
    return atom.name and atom.name[0] in ACCEPTOR_PREFIXES and not atom.name.startswith("H")


def is_hydrogen(atom: Atom) -> bool:
    return atom.name.startswith("H")


def angle_dha_deg(donor: np.ndarray, hydrogen: np.ndarray, acceptor: np.ndarray) -> float:
    v1 = donor - hydrogen
    v2 = acceptor - hydrogen
    denom = np.linalg.norm(v1) * np.linalg.norm(v2)
    if denom == 0:
        return np.nan
    cosang = np.clip(np.dot(v1, v2) / denom, -1.0, 1.0)
    return float(np.degrees(np.arccos(cosang)))


def infer_donor_hydrogens(atoms: List[Atom], max_dh_a: float) -> Dict[int, List[Atom]]:
    donors = [a for a in atoms if is_heavy_donor(a)]
    hydrogens = [a for a in atoms if is_hydrogen(a)]
    mapping: Dict[int, List[Atom]] = {d.index0: [] for d in donors}
    for h in hydrogens:
        same_res_donors = [d for d in donors if d.chain == h.chain and d.resid == h.resid]
        if not same_res_donors:
            continue
        distances = [(np.linalg.norm(h.xyz_a - d.xyz_a), d) for d in same_res_donors]
        dist, donor = min(distances, key=lambda x: x[0])
        if dist <= max_dh_a:
            mapping[donor.index0].append(h)
    return mapping


def collect_hbonds(frame: dict, atoms: List[Atom], max_da_a: float, min_angle_deg: float, max_dh_a: float) -> pd.DataFrame:
    donors = [a for a in atoms if is_heavy_donor(a)]
    acceptors = [a for a in atoms if is_acceptor(a)]
    donor_h = infer_donor_hydrogens(atoms, max_dh_a)
    rows = []
    for donor in donors:
        hydrogens = donor_h.get(donor.index0, [])
        if not hydrogens:
            continue
        for acceptor in acceptors:
            if acceptor.index0 == donor.index0:
                continue
            if donor.chain == acceptor.chain and donor.resid == acceptor.resid:
                continue
            da = float(np.linalg.norm(donor.xyz_a - acceptor.xyz_a))
            if da > max_da_a:
                continue
            best = None
            for h in hydrogens:
                angle = angle_dha_deg(donor.xyz_a, h.xyz_a, acceptor.xyz_a)
                if math.isnan(angle) or angle < min_angle_deg:
                    continue
                item = (angle, h)
                if best is None or item[0] > best[0]:
                    best = item
            if best is None:
                continue
            angle, hydrogen = best
            rows.append(
                {
                    "system": frame["system"],
                    "system_label": frame["system_label"],
                    "generation_index": frame["generation_index"],
                    "generation": frame["generation"],
                    "opening_score_angstrom": frame.get("opening_score_angstrom"),
                    "donor_chain": donor.chain,
                    "donor_resid": donor.resid,
                    "donor_resname": donor.resname,
                    "donor_atom": donor.name,
                    "hydrogen_atom": hydrogen.name,
                    "acceptor_chain": acceptor.chain,
                    "acceptor_resid": acceptor.resid,
                    "acceptor_resname": acceptor.resname,
                    "acceptor_atom": acceptor.name,
                    "da_distance_angstrom": da,
                    "dha_angle_deg": angle,
                    "best_gro": frame["best_gro"],
                }
            )
    return pd.DataFrame(rows)


def summarize_frame(frame: dict, atoms: List[Atom], hbond_df: pd.DataFrame, dihedral_df: pd.DataFrame) -> dict:
    row = {
        "system": frame["system"],
        "system_label": frame["system_label"],
        "generation_index": frame["generation_index"],
        "generation": frame["generation"],
        "opening_score_angstrom": frame.get("opening_score_angstrom"),
        "best_gro": frame["best_gro"],
        "n_collagen_atoms": len(atoms),
        "n_hbonds": len(hbond_df),
        "mean_hbond_da_angstrom": np.nan,
        "median_hbond_da_angstrom": np.nan,
        "mean_hbond_angle_deg": np.nan,
    }
    if not hbond_df.empty:
        row["mean_hbond_da_angstrom"] = float(hbond_df["da_distance_angstrom"].mean())
        row["median_hbond_da_angstrom"] = float(hbond_df["da_distance_angstrom"].median())
        row["mean_hbond_angle_deg"] = float(hbond_df["dha_angle_deg"].mean())
    for dih in ["phi", "psi", "omega"]:
        subset = dihedral_df[dihedral_df["dihedral"] == dih] if not dihedral_df.empty else pd.DataFrame()
        row[f"{dih}_mean_deg"] = np.nan if subset.empty else float(subset["angle_deg"].mean())
        row[f"{dih}_circular_mean_deg"] = np.nan if subset.empty else circular_mean_deg(subset["angle_deg"].to_numpy())
        row[f"{dih}_circular_std_deg"] = np.nan if subset.empty else circular_std_deg(subset["angle_deg"].to_numpy())
    return row


def circular_mean_deg(values: np.ndarray) -> float:
    radians = np.radians(values.astype(float))
    return float(np.degrees(np.arctan2(np.nanmean(np.sin(radians)), np.nanmean(np.cos(radians)))))


def circular_std_deg(values: np.ndarray) -> float:
    radians = np.radians(values.astype(float))
    s = np.nanmean(np.sin(radians))
    c = np.nanmean(np.cos(radians))
    r = min(max(math.sqrt(s * s + c * c), 1.0e-12), 1.0)
    return float(np.degrees(math.sqrt(-2.0 * math.log(r))))


def analyze_frames(frame_df: pd.DataFrame, max_da_a: float, min_angle_deg: float, max_dh_a: float) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    hbond_tables = []
    dihedral_tables = []
    summary_rows = []
    range_cache: Dict[str, Dict[str, Tuple[int, int, str]]] = {}
    for _, frame in frame_df.iterrows():
        run_dir = Path(frame["run_dir"])
        input_dir = get_input_dir(run_dir)
        key = str(input_dir)
        if key not in range_cache:
            range_cache[key] = collagen_chain_ranges(input_dir)
        atoms = assign_collagen_chains(parse_gro(Path(frame["best_gro"])), range_cache[key])
        residues = group_residues(atoms)
        frame_dict = frame.to_dict()
        hbond_df = collect_hbonds(frame_dict, atoms, max_da_a=max_da_a, min_angle_deg=min_angle_deg, max_dh_a=max_dh_a)
        dihedral_df = collect_dihedrals(frame_dict, residues)
        if not hbond_df.empty:
            hbond_tables.append(hbond_df)
        if not dihedral_df.empty:
            dihedral_tables.append(dihedral_df)
        summary_rows.append(summarize_frame(frame_dict, atoms, hbond_df, dihedral_df))
    hbonds = pd.concat(hbond_tables, ignore_index=True) if hbond_tables else pd.DataFrame()
    dihedrals = pd.concat(dihedral_tables, ignore_index=True) if dihedral_tables else pd.DataFrame()
    summary = pd.DataFrame(summary_rows)
    return hbonds, dihedrals, summary


def add_text_page(pdf: PdfPages, title: str, lines: List[str]) -> None:
    fig = plt.figure(figsize=(8.27, 11.69))
    ax = fig.add_subplot(111)
    ax.axis("off")
    ax.text(0.06, 0.95, title, fontsize=16, weight="bold", va="top")
    ax.text(0.06, 0.89, "\n".join(lines), fontsize=10.5, va="top", family="monospace")
    pdf.savefig(fig)
    plt.close(fig)


def plot_hbond_counts(pdf: PdfPages, output_dir: Path, summary: pd.DataFrame) -> None:
    if summary.empty:
        return
    fig, ax = plt.subplots(figsize=(8.27, 5.85))
    for system, meta in SYSTEMS.items():
        sdf = summary[summary["system"] == system].sort_values("generation_index")
        if sdf.empty:
            continue
        ax.plot(sdf["generation_index"], sdf["n_hbonds"], marker="o", color=meta["color"], label=meta["label"])
    ax.set_xlabel("Adaptive generation")
    ax.set_ylabel("Collagen hydrogen bonds")
    ax.set_title("Hydrogen-bond count across best workers")
    ax.grid(True, alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_dir / "structural_hbond_count_by_generation.png", dpi=300)
    pdf.savefig(fig)
    plt.close(fig)


def plot_hbond_distance_distribution(pdf: PdfPages, output_dir: Path, hbonds: pd.DataFrame) -> None:
    if hbonds.empty:
        return
    fig, ax = plt.subplots(figsize=(8.27, 5.85))
    bins = np.linspace(2.2, 3.6, 36)
    for system, meta in SYSTEMS.items():
        vals = hbonds.loc[hbonds["system"] == system, "da_distance_angstrom"].dropna()
        if vals.empty:
            continue
        ax.hist(vals, bins=bins, histtype="step", density=True, linewidth=1.8, color=meta["color"], label=meta["label"])
    ax.set_xlabel("Donor-acceptor distance / A")
    ax.set_ylabel("Density")
    ax.set_title("Hydrogen-bond distance distribution")
    ax.grid(True, alpha=0.2)
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_dir / "structural_hbond_distance_distribution.png", dpi=300)
    pdf.savefig(fig)
    plt.close(fig)


def plot_dihedral_distribution(pdf: PdfPages, output_dir: Path, dihedrals: pd.DataFrame) -> None:
    if dihedrals.empty:
        return
    bins = np.linspace(-180, 180, 73)
    for dih in ["phi", "psi", "omega"]:
        fig, ax = plt.subplots(figsize=(8.27, 5.85))
        for system, meta in SYSTEMS.items():
            vals = dihedrals.loc[(dihedrals["system"] == system) & (dihedrals["dihedral"] == dih), "angle_deg"].dropna()
            if vals.empty:
                continue
            ax.hist(vals, bins=bins, histtype="step", density=True, linewidth=1.8, color=meta["color"], label=meta["label"])
        ax.set_xlabel(f"{dih} angle / degrees")
        ax.set_ylabel("Density")
        ax.set_title(f"Collagen backbone {dih} distribution")
        ax.grid(True, alpha=0.2)
        ax.legend()
        fig.tight_layout()
        fig.savefig(output_dir / f"structural_dihedral_{dih}_distribution.png", dpi=300)
        pdf.savefig(fig)
        plt.close(fig)


def scatter_vs_unwinding(pdf: PdfPages, output_dir: Path, summary: pd.DataFrame) -> None:
    if summary.empty:
        return
    metrics = [
        ("n_hbonds", "Hydrogen-bond count"),
        ("mean_hbond_da_angstrom", "Mean H-bond donor-acceptor distance / A"),
        ("phi_circular_std_deg", "Phi circular std / degrees"),
        ("psi_circular_std_deg", "Psi circular std / degrees"),
        ("omega_circular_std_deg", "Omega circular std / degrees"),
    ]
    for col, label in metrics:
        if col not in summary.columns:
            continue
        has_data = False
        fig, ax = plt.subplots(figsize=(8.27, 5.85))
        for system, meta in SYSTEMS.items():
            sdf = summary[(summary["system"] == system)].dropna(subset=["opening_score_angstrom", col])
            if sdf.empty:
                continue
            has_data = True
            ax.scatter(sdf["opening_score_angstrom"], sdf[col], color=meta["color"], label=meta["label"], s=45)
            if len(sdf) >= 2:
                x = sdf["opening_score_angstrom"].to_numpy(dtype=float)
                y = sdf[col].to_numpy(dtype=float)
                order = np.argsort(x)
                ax.plot(x[order], y[order], color=meta["color"], alpha=0.45)
        if not has_data:
            plt.close(fig)
            continue
        ax.set_xlabel("Best opening score / A")
        ax.set_ylabel(label)
        ax.set_title(f"{label} as a function of unwinding")
        ax.grid(True, alpha=0.25)
        ax.legend()
        fig.tight_layout()
        fig.savefig(output_dir / f"structural_vs_unwinding_{col}.png", dpi=300)
        pdf.savefig(fig)
        plt.close(fig)


def generate_pdf(output_pdf: Path, output_dir: Path, specs: List[RunSpec], frames: pd.DataFrame, hbonds: pd.DataFrame, dihedrals: pd.DataFrame, summary: pd.DataFrame, args) -> None:
    lines = [
        "Joint collagen structural report",
        "",
        "Inputs: best worker/frame selected for each adaptive generation.",
        "Scope: collagen chains only.",
        f"Hydrogen-bond criterion: donor-acceptor <= {args.max_hbond_da_a:.2f} A; D-H-A angle >= {args.min_hbond_angle_deg:.1f} deg.",
        "Backbone dihedrals: phi, psi, omega for consecutive collagen residues.",
        "",
        "Runs:",
    ]
    for spec in specs:
        lines.append(f"  {spec.label:10s} {spec.run_dir}")
    lines.extend(
        [
            "",
            f"Best frames analyzed: {len(frames)}",
            f"Hydrogen bonds detected: {len(hbonds)}",
            f"Dihedral rows: {len(dihedrals)}",
        ]
    )
    with PdfPages(output_pdf) as pdf:
        add_text_page(pdf, "Joint Collagen Structural Report", lines)
        plot_hbond_counts(pdf, output_dir, summary)
        plot_hbond_distance_distribution(pdf, output_dir, hbonds)
        plot_dihedral_distribution(pdf, output_dir, dihedrals)
        scatter_vs_unwinding(pdf, output_dir, summary)


def write_outputs(output_dir: Path, frames: pd.DataFrame, hbonds: pd.DataFrame, dihedrals: pd.DataFrame, summary: pd.DataFrame) -> None:
    safe_mkdir(output_dir)
    if hbonds.empty:
        hbonds = pd.DataFrame(
            columns=[
                "system",
                "system_label",
                "generation_index",
                "generation",
                "opening_score_angstrom",
                "donor_chain",
                "donor_resid",
                "donor_resname",
                "donor_atom",
                "hydrogen_atom",
                "acceptor_chain",
                "acceptor_resid",
                "acceptor_resname",
                "acceptor_atom",
                "da_distance_angstrom",
                "dha_angle_deg",
                "best_gro",
            ]
        )
    if dihedrals.empty:
        dihedrals = pd.DataFrame(
            columns=[
                "system",
                "system_label",
                "generation_index",
                "generation",
                "opening_score_angstrom",
                "chain",
                "resid",
                "resname",
                "best_gro",
                "dihedral",
                "angle_deg",
            ]
        )
    frames.to_csv(output_dir / "structural_best_frames.csv", index=False)
    summary.to_csv(output_dir / "structural_generation_summary.csv", index=False)
    hbonds.to_csv(output_dir / "structural_hydrogen_bonds.csv", index=False)
    dihedrals.to_csv(output_dir / "structural_backbone_dihedrals.csv", index=False)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate a joint structural PDF/CSV report from best adaptive collagen-unwinding frames.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--root", type=str, default=".", help="Repository root.")
    parser.add_argument("--wild_type_run", type=str, default=None, help="Explicit wild-type adaptive_run_* directory.")
    parser.add_argument("--g978s_run", type=str, default=None, help="Explicit G978S adaptive_run_* directory.")
    parser.add_argument("--g984c_run", type=str, default=None, help="Explicit G984C adaptive_run_* directory.")
    parser.add_argument("--g987r_run", type=str, default=None, help="Explicit G987R adaptive_run_* directory.")
    parser.add_argument("--analysis_dir", type=str, default="joint_structural_report", help="Output directory for CSV and PNG files.")
    parser.add_argument("--output", type=str, default=None, help="Output PDF path.")
    parser.add_argument("--max_hbond_da_a", type=float, default=3.5, help="Maximum donor-acceptor distance for hydrogen bonds in Angstrom.")
    parser.add_argument("--min_hbond_angle_deg", type=float, default=120.0, help="Minimum D-H-A angle for hydrogen bonds.")
    parser.add_argument("--max_dh_a", type=float, default=1.25, help="Maximum inferred donor-hydrogen bond length in Angstrom.")
    args = parser.parse_args()

    specs = resolve_run_specs(args)
    if not specs:
        raise SystemExit("No adaptive_run_* directories found. Pass explicit run directories.")
    output_dir = Path(args.analysis_dir).resolve()
    output_pdf = Path(args.output).resolve() if args.output else output_dir / "joint_collagen_structural_report.pdf"
    safe_mkdir(output_dir)

    frames = collect_best_frame_records(specs)
    if frames.empty:
        raise SystemExit("No best_opening_frame.gro records found in the supplied runs.")
    hbonds, dihedrals, summary = analyze_frames(frames, args.max_hbond_da_a, args.min_hbond_angle_deg, args.max_dh_a)
    write_outputs(output_dir, frames, hbonds, dihedrals, summary)
    generate_pdf(output_pdf, output_dir, specs, frames, hbonds, dihedrals, summary, args)

    print(f"Systems analysed: {', '.join(spec.system for spec in specs)}")
    print(f"Best frames analysed: {len(frames)}")
    print(f"Hydrogen bonds detected: {len(hbonds)}")
    print(f"Dihedral rows: {len(dihedrals)}")
    print(f"PDF written: {output_pdf}")
    print(f"Analysis directory: {output_dir}")


if __name__ == "__main__":
    main()
