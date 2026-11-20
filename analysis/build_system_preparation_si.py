#!/usr/bin/env python3
"""Analyse saved MMP-1/collagen preparation evidence and build SI figures/tables."""

from __future__ import annotations

import csv
import json
import math
import re
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SYSTEMS = ROOT / "systems" / "mmp1-bound-adaptive"
GENERATED = ROOT / "generated_mutants"
OUT = ROOT / "manuscript" / "si_figures"
DATA_OUT = ROOT / "analysis" / "si_system_preparation_metrics.json"

COLOURS = {
    "Wild type": "#111111",
    "G978S": "#2574A9",
    "G984C": "#2E8B57",
    "G987R": "#C0392B",
}


@dataclass
class GroAtom:
    index: int
    resnr: int
    resname: str
    name: str
    xyz: np.ndarray


def set_style() -> None:
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
        "font.size": 9,
        "axes.labelsize": 9,
        "axes.titlesize": 10,
        "axes.linewidth": 0.8,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 8,
        "figure.dpi": 180,
        "savefig.dpi": 300,
    })


def save(fig: plt.Figure, name: str) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / name, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(-0.14, 1.08, label, transform=ax.transAxes, fontweight="bold", fontsize=11,
            ha="left", va="top")


def parse_gro(path: Path) -> tuple[list[GroAtom], np.ndarray]:
    lines = path.read_text(encoding="utf-8").splitlines()
    n_atoms = int(lines[1].strip())
    atoms = []
    for i, line in enumerate(lines[2:2 + n_atoms]):
        atoms.append(GroAtom(
            index=i,
            resnr=int(line[:5]),
            resname=line[5:10].strip(),
            name=line[10:15].strip(),
            xyz=np.array([float(line[20:28]), float(line[28:36]), float(line[36:44])]) * 10.0,
        ))
    box = np.array([float(x) for x in lines[2 + n_atoms].split()[:3]])
    return atoms, box


def ca_groups(atoms: list[GroAtom]) -> dict[str, np.ndarray]:
    solvent_names = {"SOL", "NA", "CL"}
    protein_end = next((a.index for a in atoms if a.resname in solvent_names), len(atoms))
    all_ca = [a for a in atoms if a.index < protein_end and a.name == "CA"]
    mmp = [a for a in all_ca if a.index < 5825]
    collagen = [a for a in all_ca if a.index >= 5825]
    # OpenMM rewrites residue identifiers in exported GRO files. The analysed
    # segment is therefore selected by its invariant position within each
    # 33-residue collagen chain: residues 15-25 correspond to model 977-987.
    if len(collagen) % 33:
        raise ValueError(f"Unexpected collagen C-alpha count: {len(collagen)}")
    local = []
    for chain_start in range(0, len(collagen), 33):
        local.extend(collagen[chain_start + 14:chain_start + 25])
    return {
        "protein": np.array([a.xyz for a in all_ca]),
        "MMP-1": np.array([a.xyz for a in mmp]),
        "collagen": np.array([a.xyz for a in collagen]),
        "collagen 977-987": np.array([a.xyz for a in local]),
    }


def align_transform(reference: np.ndarray, mobile: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    ref_c = reference.mean(axis=0)
    mob_c = mobile.mean(axis=0)
    ref0 = reference - ref_c
    mob0 = mobile - mob_c
    u, _, vt = np.linalg.svd(mob0.T @ ref0)
    rot = u @ vt
    if np.linalg.det(rot) < 0:
        u[:, -1] *= -1
        rot = u @ vt
    return rot, mob_c, ref_c


def transformed(coords: np.ndarray, rot: np.ndarray, mob_c: np.ndarray, ref_c: np.ndarray) -> np.ndarray:
    return (coords - mob_c) @ rot + ref_c


def rmsd(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.sqrt(np.mean(np.sum((a - b) ** 2, axis=1))))


def unwrap_chain(coords: np.ndarray, box: np.ndarray) -> np.ndarray:
    out = np.empty_like(coords)
    out[0] = coords[0]
    for i in range(1, len(coords)):
        delta = coords[i] - coords[i - 1]
        delta -= box * np.round(delta / box)
        out[i] = out[i - 1] + delta
    return out


def aligned_chain_rmsd(reference: np.ndarray, mobile: np.ndarray, ref_box: np.ndarray,
                       mob_box: np.ndarray) -> float:
    ref = unwrap_chain(reference, ref_box * 10.0)
    mob = unwrap_chain(mobile, mob_box * 10.0)
    rot, mob_c, ref_c = align_transform(ref, mob)
    return rmsd(ref, transformed(mob, rot, mob_c, ref_c))


def mean_chunk_rmsd(reference: np.ndarray, mobile: np.ndarray, chunk: int,
                    ref_box: np.ndarray, mob_box: np.ndarray) -> float:
    values = []
    for start in range(0, len(reference), chunk):
        values.append(aligned_chain_rmsd(reference[start:start + chunk], mobile[start:start + chunk],
                                         ref_box, mob_box))
    return float(np.mean(values))


def distance_rmsd(reference: np.ndarray, mobile: np.ndarray, ref_box: np.ndarray,
                  mob_box: np.ndarray) -> float:
    def matrix(coords: np.ndarray, box_nm: np.ndarray) -> np.ndarray:
        delta = coords[:, None, :] - coords[None, :, :]
        box_a = box_nm * 10.0
        delta -= box_a * np.round(delta / box_a)
        return np.sqrt(np.sum(delta * delta, axis=-1))
    dr = matrix(reference, ref_box)
    dm = matrix(mobile, mob_box)
    tri = np.triu_indices(len(reference), k=1)
    return float(np.sqrt(np.mean((dr[tri] - dm[tri]) ** 2)))


def structural_metrics(reference_path: Path, mobile_path: Path) -> dict[str, float]:
    ref_atoms, ref_box = parse_gro(reference_path)
    mob_atoms, mob_box = parse_gro(mobile_path)
    ref = ca_groups(ref_atoms)
    mob = ca_groups(mob_atoms)
    if ref["protein"].shape != mob["protein"].shape:
        raise ValueError(f"C-alpha mismatch: {reference_path.name} vs {mobile_path.name}")
    values = {
        "MMP-1_ca_rmsd_angstrom": aligned_chain_rmsd(ref["MMP-1"], mob["MMP-1"], ref_box, mob_box),
        "collagen_chain_mean_ca_rmsd_angstrom": mean_chunk_rmsd(
            ref["collagen"], mob["collagen"], 33, ref_box, mob_box),
        "collagen_977-987_chain_mean_ca_rmsd_angstrom": mean_chunk_rmsd(
            ref["collagen 977-987"], mob["collagen 977-987"], 11, ref_box, mob_box),
        "collagen_977-987_distance_rmsd_angstrom": distance_rmsd(
            ref["collagen 977-987"], mob["collagen 977-987"], ref_box, mob_box),
    }
    values["initial_box_volume_nm3"] = float(np.prod(ref_box))
    values["final_box_volume_nm3"] = float(np.prod(mob_box))
    return values


def parse_obminimize(path: Path) -> tuple[np.ndarray, np.ndarray]:
    steps, energies = [], []
    pattern = re.compile(r"^\s*(\d+)\s+(-?\d+(?:\.\d+)?)\s+(-?\d+(?:\.\d+)?)\s*$")
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        match = pattern.match(line)
        if match:
            steps.append(int(match.group(1)))
            energies.append(float(match.group(2)))
    return np.asarray(steps), np.asarray(energies)


def parse_csv_log(path: Path) -> dict[str, np.ndarray]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
    result = {}
    for key in reader.fieldnames or []:
        clean = key.lstrip("#")
        result[clean] = np.asarray([float(r[key]) for r in rows])
    return result


def parse_state(path: Path) -> dict[str, float]:
    root = ET.parse(path).getroot()
    energies = root.find("Energies")
    vectors = root.find("PeriodicBoxVectors")
    if energies is None or vectors is None:
        raise ValueError(f"Missing state fields in {path}")
    vecs = []
    for name in ("A", "B", "C"):
        elem = vectors.find(name)
        if elem is None:
            raise ValueError(f"Missing box vector {name} in {path}")
        vecs.append([float(elem.attrib[k]) for k in ("x", "y", "z")])
    return {
        "potential_energy_kj_mol": float(energies.attrib["PotentialEnergy"]),
        "kinetic_energy_kj_mol": float(energies.attrib["KineticEnergy"]),
        "box_volume_nm3": float(abs(np.linalg.det(np.asarray(vecs)))),
        "time_ps": float(root.attrib["time"]),
        "step_count": int(root.attrib["stepCount"]),
    }


def build_workflow() -> None:
    fig, ax = plt.subplots(figsize=(8.2, 3.8))
    ax.set_xlim(0, 14)
    ax.set_ylim(0, 6)
    ax.axis("off")
    boxes = [
        (0.3, 4.3, 2.5, 1.0, "4AUO-derived\nMMP-1-collagen\ncomplex"),
        (3.8, 4.3, 2.5, 1.0, "Insert one OI-linked\nsubstitution in\ncollagen chain 1"),
        (7.3, 4.3, 2.5, 1.0, "Transfer bonded terms;\nchoose a low-clash\nside-chain rotamer"),
        (10.8, 4.3, 2.5, 1.0, "Set charge and add\n150 mM NaCl"),
        (10.8, 1.6, 2.5, 1.0, "Final 10 ps NVT +\n10 ps NPT state\nwith velocities"),
        (7.3, 1.6, 2.5, 1.0, "Generation 0:\nminimise, 50 ps NVT,\n200 ps NPT"),
        (3.8, 1.6, 2.5, 1.0, "20 ns production in\nten workers; select\nmaximum separation"),
        (0.3, 1.6, 2.5, 1.0, "Later generations:\nparent geometry, 50 ps\nNVT, 200 ps NPT"),
    ]
    for x, y, w, h, label in boxes:
        rect = plt.Rectangle((x, y), w, h, facecolor="#F3F5F7", edgecolor="#333333", lw=0.9)
        ax.add_patch(rect)
        ax.text(x + w / 2, y + h / 2, label, ha="center", va="center", fontsize=7.6)
    arrows = [
        ((2.8, 4.8), (3.8, 4.8)), ((6.3, 4.8), (7.3, 4.8)), ((9.8, 4.8), (10.8, 4.8)),
        ((12.05, 4.3), (12.05, 2.6)), ((10.8, 2.1), (9.8, 2.1)), ((7.3, 2.1), (6.3, 2.1)),
        ((3.8, 2.1), (2.8, 2.1)),
    ]
    for start, end in arrows:
        ax.annotate("", xy=end, xytext=start, arrowprops=dict(arrowstyle="->", lw=1.0, color="#333333"))
    ax.annotate("next selected parent", xy=(1.55, 2.65), xytext=(4.7, 3.35), ha="center", fontsize=7.5,
                arrowprops=dict(arrowstyle="->", lw=0.9, color="#666666", connectionstyle="arc3,rad=0.25"))
    ax.text(7.0, 5.72, "Starting-system construction and adaptive propagation", ha="center",
            va="top", fontsize=11, fontweight="bold")
    save(fig, "figure_S1_preparation_workflow.png")


def build_obminimize_plot() -> dict[str, dict[str, float]]:
    fig, ax = plt.subplots(figsize=(6.9, 3.7))
    summary = {}
    for variant in ("G978S", "G984C", "G987R"):
        steps, energy = parse_obminimize(GENERATED / f"collagen_{variant}" / "obminimize.log")
        summary[variant] = {
            "initial_recorded_energy": float(energy[0]),
            "final_recorded_energy": float(energy[-1]),
            "iterations": int(steps[-1]),
        }
        ax.plot(steps, energy, color=COLOURS[variant], lw=1.8, marker="o", ms=3.2, label=variant)
    ax.set_xlabel("Minimisation iteration")
    ax.set_ylabel("Open Babel reported energy")
    ax.set_yscale("log")
    ax.set_title("Local cluster minimisation records")
    ax.grid(axis="y", color="#D9D9D9", lw=0.6)
    ax.legend(frameon=False, ncol=3, loc="upper right")
    ax.spines[["top", "right"]].set_visible(False)
    save(fig, "figure_S2_local_cluster_minimisation.png")
    return summary


def build_g978s_thermodynamics() -> dict[str, float]:
    base = SYSTEMS / "collagen_G978S" / "openmm_validation"
    nvt = parse_csv_log(base / "nvt.log")
    npt = parse_csv_log(base / "npt.log")
    summary = json.loads((base / "summary.json").read_text(encoding="utf-8"))
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.5))
    ax = axes[0, 0]
    ax.plot(nvt["Time (ps)"], nvt["Potential Energy (kJ/mole)"], color=COLOURS["G978S"], lw=1.6)
    ax.set(xlabel="NVT time (ps)", ylabel="Potential energy (kJ mol⁻¹)")
    panel_label(ax, "A")
    ax = axes[0, 1]
    ax.plot(nvt["Time (ps)"], nvt["Temperature (K)"], color=COLOURS["G978S"], lw=1.6)
    ax.axhline(summary["temperature_k"], color="#666666", ls="--", lw=0.9, label="Target")
    ax.set(xlabel="NVT time (ps)", ylabel="Temperature (K)")
    ax.legend(frameon=False)
    panel_label(ax, "B")
    ax = axes[1, 0]
    ax.plot(npt["Time (ps)"], npt["Potential Energy (kJ/mole)"], color=COLOURS["G978S"], lw=1.6)
    ax.set(xlabel="NPT time (ps)", ylabel="Potential energy (kJ mol⁻¹)")
    panel_label(ax, "C")
    ax = axes[1, 1]
    ax.plot(npt["Time (ps)"], npt["Box Volume (nm^3)"], color=COLOURS["G978S"], lw=1.6, label="Volume")
    ax.set(xlabel="NPT time (ps)", ylabel="Box volume (nm³)")
    ax2 = ax.twinx()
    ax2.plot(npt["Time (ps)"], npt["Density (g/mL)"], color="#6B6B6B", lw=1.2, ls="--", label="Density")
    ax2.set_ylabel("Density (g mL⁻¹)")
    lines = ax.get_lines() + ax2.get_lines()
    ax.legend(lines, [line.get_label() for line in lines], frameon=False, loc="center right")
    panel_label(ax, "D")
    for a in axes.flat:
        a.grid(axis="y", color="#E0E0E0", lw=0.5)
        a.spines[["top", "right"]].set_visible(False)
    fig.suptitle("Saved G978S NVT and NPT validation record", fontsize=11, fontweight="bold", y=1.01)
    fig.tight_layout()
    save(fig, "figure_S3_G978S_validation_thermodynamics.png")
    return {
        "nvt_temperature_mean_k": float(np.mean(nvt["Temperature (K)"])),
        "nvt_temperature_final_k": float(nvt["Temperature (K)"][-1]),
        "npt_temperature_mean_k": float(np.mean(npt["Temperature (K)"])),
        "npt_volume_mean_nm3": float(np.mean(npt["Box Volume (nm^3)"])),
        "npt_density_mean_g_ml": float(np.mean(npt["Density (g/mL)"])),
        "minimized_potential_kj_mol": float(summary["minimized_potential_kj_per_mol"]),
        "post_nvt_potential_kj_mol": float(summary["post_nvt_potential_kj_per_mol"]),
        "post_npt_potential_kj_mol": float(summary["post_npt_potential_kj_per_mol"]),
    }


def build_structural_change_plot() -> tuple[dict, dict]:
    final_metrics = {}
    for variant in ("G978S", "G984C", "G987R"):
        base = SYSTEMS / f"collagen_{variant}"
        initial = base / f"NPT_eq_collagen_{variant}_150mM_NaCl.gro"
        final = base / f"NPT_eq_collagen_{variant}_150mM_NaCl_openmm_10ps_nvt_10ps_npt_with_velocities.gro"
        final_metrics[variant] = structural_metrics(initial, final)

    gbase = SYSTEMS / "collagen_G978S"
    initial = gbase / "NPT_eq_collagen_G978S_150mM_NaCl.gro"
    stage_paths = {
        "After minimisation": gbase / "openmm_validation" / "G978S_after_minimization.gro",
        "After 10 ps NVT": gbase / "openmm_validation" / "G978S_after_10ps_nvt.gro",
        "After 10 ps NPT": gbase / "openmm_validation" / "NPT_eq_collagen_G978S_150mM_NaCl_openmm_10ps_npt.gro",
    }
    stage_metrics = {label: structural_metrics(initial, path) for label, path in stage_paths.items()}

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.7))
    groups = ["MMP-1_ca_rmsd_angstrom", "collagen_chain_mean_ca_rmsd_angstrom",
              "collagen_977-987_chain_mean_ca_rmsd_angstrom", "collagen_977-987_distance_rmsd_angstrom"]
    labels = ["MMP-1", "Collagen\nchain", "Local collagen\nchain", "Local packing\ndRMSD"]
    x = np.arange(len(groups))
    width = 0.24
    for i, variant in enumerate(("G978S", "G984C", "G987R")):
        vals = [final_metrics[variant][g] for g in groups]
        axes[0].bar(x + (i - 1) * width, vals, width, color=COLOURS[variant], label=variant)
    axes[0].set_xticks(x, labels)
    axes[0].set_ylabel("Displacement from inserted structure (Å)")
    axes[0].set_title("Final saved starting states")
    axes[0].legend(frameon=False, ncol=3, loc="upper left")
    panel_label(axes[0], "A")

    stage_labels = list(stage_metrics)
    sx = np.arange(len(stage_labels))
    for group, label, colour, marker in [
        ("MMP-1_ca_rmsd_angstrom", "MMP-1", "#4D4D4D", "o"),
        ("collagen_chain_mean_ca_rmsd_angstrom", "Collagen chain mean", COLOURS["G978S"], "s"),
        ("collagen_977-987_chain_mean_ca_rmsd_angstrom", "Residues 977-987 chain mean", "#7B4AB1", "^")
    ]:
        vals = [stage_metrics[s][group] for s in stage_labels]
        axes[1].plot(sx, vals, color=colour, marker=marker, lw=1.6, ms=5, label=label)
    axes[1].set_xticks(sx, ["Minimised", "10 ps\nNVT", "10 ps\nNPT"])
    axes[1].set_ylabel("Cα RMSD from inserted G978S structure (Å)")
    axes[1].set_title("Saved G978S preparation stages")
    axes[1].legend(frameon=False, loc="upper left")
    panel_label(axes[1], "B")
    for ax in axes:
        ax.grid(axis="y", color="#E0E0E0", lw=0.5)
        ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    save(fig, "figure_S4_preparation_structural_change.png")
    return final_metrics, stage_metrics


def build_final_state_plot() -> dict[str, dict[str, float]]:
    states = {}
    for variant in ("G978S", "G984C", "G987R"):
        path = SYSTEMS / f"collagen_{variant}" / "openmm_final_10ps_10ps_with_velocities" / f"{variant}_final_state.xml"
        states[variant] = parse_state(path)
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 3.1))
    variants = list(states)
    cols = [COLOURS[v] for v in variants]
    values = [
        ([states[v]["potential_energy_kj_mol"] / 1000 for v in variants], "Potential energy\n(10³ kJ mol⁻¹)"),
        ([states[v]["kinetic_energy_kj_mol"] / 1000 for v in variants], "Kinetic energy\n(10³ kJ mol⁻¹)"),
        ([states[v]["box_volume_nm3"] for v in variants], "Box volume (nm³)"),
    ]
    for i, (ax, (vals, ylabel)) in enumerate(zip(axes, values)):
        ax.bar(np.arange(3), vals, color=cols, width=0.68)
        ax.set_xticks(np.arange(3), variants, rotation=30, ha="right")
        ax.set_ylabel(ylabel)
        ax.grid(axis="y", color="#E0E0E0", lw=0.5)
        ax.spines[["top", "right"]].set_visible(False)
        panel_label(ax, chr(65 + i))
    fig.suptitle("Quantities stored in the final 10 ps NPT states", fontsize=11, fontweight="bold", y=1.02)
    fig.tight_layout()
    save(fig, "figure_S5_final_state_quantities.png")
    return states


def main() -> None:
    set_style()
    build_workflow()
    ob = build_obminimize_plot()
    thermo = build_g978s_thermodynamics()
    final_metrics, stage_metrics = build_structural_change_plot()
    states = build_final_state_plot()
    payload = {
        "local_cluster_minimisation": ob,
        "G978S_validation": thermo,
        "initial_to_final_structural_metrics": final_metrics,
        "G978S_stage_structural_metrics": stage_metrics,
        "final_openmm_states": states,
    }
    DATA_OUT.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(DATA_OUT)
    for path in sorted(OUT.glob("figure_S*.png")):
        print(path)


if __name__ == "__main__":
    main()
