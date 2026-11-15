#!/usr/bin/env python3
"""Generation-matched manuscript analysis of the complete structural export."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


SYSTEM_ORDER = ["wild_type", "G978S", "G984C", "G987R"]
LABELS = {"wild_type": "Wild type", "G978S": "G978S", "G984C": "G984C", "G987R": "G987R"}
COLORS = {"wild_type": "#303030", "G978S": "#2574A9", "G984C": "#2E8B57", "G987R": "#C0392B"}
WINDOW = ["ALA", "GLY", "GLN", "ARG", "GLY", "ILE", "VAL", "GLY", "LEU", "HYP", "GLY"]
MUTATIONS = {"G978S": (1, "SER"), "G984C": (7, "CYS"), "G987R": (10, "ARG")}
COMMON_MAX_GENERATION = 28


def expected_window(system: str, chain: str) -> list[str]:
    names = list(WINDOW)
    if chain == "chain_1" and system in MUTATIONS:
        index, replacement = MUTATIONS[system]
        names[index] = replacement
    return names


def find_window_mapping(dihedrals: pd.DataFrame) -> dict[tuple[str, str, int], int]:
    mapping: dict[tuple[str, str, int], int] = {}
    initial = dihedrals[dihedrals["generation_index"] == 0]
    for system in SYSTEM_ORDER:
        for chain in ["chain_1", "chain_2", "chain_3"]:
            rows = (
                initial[(initial["system"] == system) & (initial["chain"] == chain)]
                .drop_duplicates(["resid"])
                .sort_values("resid")
            )
            residues = rows["resid"].astype(int).tolist()
            names = rows["resname"].tolist()
            expected = expected_window(system, chain)
            start = None
            for index in range(len(names) - len(expected) + 1):
                if names[index : index + len(expected)] == expected:
                    start = index
                    break
            if start is None:
                raise RuntimeError(f"Could not locate residues 977-987 for {system} {chain}")
            for offset, position in enumerate(range(977, 988)):
                mapping[(system, chain, residues[start + offset])] = position
    return mapping


def circular_std_deg(values: pd.Series) -> float:
    radians = np.radians(values.dropna().astype(float).to_numpy())
    if len(radians) == 0:
        return float("nan")
    s = np.mean(np.sin(radians))
    c = np.mean(np.cos(radians))
    radius = min(max(math.sqrt(s * s + c * c), 1.0e-12), 1.0)
    return float(np.degrees(math.sqrt(-2.0 * math.log(radius))))


def early_late_summary(df: pd.DataFrame, metrics: list[str]) -> pd.DataFrame:
    rows = []
    for system in SYSTEM_ORDER:
        group = df[df["system"] == system].sort_values("generation_index")
        early = group[group["generation_index"].between(0, 4)]
        late = group[group["generation_index"].between(24, 28)]
        for metric in metrics:
            early_mean = float(early[metric].mean())
            late_mean = float(late[metric].mean())
            rows.append(
                {
                    "system": system,
                    "system_label": LABELS[system],
                    "metric": metric,
                    "early_mean_generations_0_4": early_mean,
                    "late_mean_generations_24_28": late_mean,
                    "late_minus_early": late_mean - early_mean,
                    "pearson_r_with_opening_score": float(group[metric].corr(group["opening_score_angstrom"])),
                }
            )
    return pd.DataFrame(rows)


def selected_frame_chain_pair_opening(geometry: pd.DataFrame) -> pd.DataFrame:
    columns = {
        "chain_1_chain_2_ca_distance_angstrom": "chain_1-chain_2",
        "chain_1_chain_3_ca_distance_angstrom": "chain_1-chain_3",
        "chain_2_chain_3_ca_distance_angstrom": "chain_2-chain_3",
    }
    rows = []
    for (system, generation), group in geometry.groupby(["system", "generation_index"]):
        for column, pair in columns.items():
            values = group[column].astype(float)
            rows.append(
                {
                    "system": system,
                    "system_label": LABELS[system],
                    "generation_index": int(generation),
                    "chain_pair": pair,
                    "mean_ca_distance_angstrom": float(values.mean()),
                    "maximum_ca_distance_angstrom": float(values.max()),
                    "minimum_ca_distance_angstrom": float(values.min()),
                    "selected_structure_only": True,
                }
            )
    return pd.DataFrame(rows)


def selected_frame_chain_pair_hbonds(within_window: pd.DataFrame) -> pd.DataFrame:
    data = within_window.copy()
    data["chain_pair"] = ["-".join(sorted((str(row.donor_chain), str(row.acceptor_chain)))) for row in data.itertuples()]
    index = pd.MultiIndex.from_product(
        [SYSTEM_ORDER, range(COMMON_MAX_GENERATION + 1), ["chain_1-chain_2", "chain_1-chain_3", "chain_2-chain_3"]],
        names=["system", "generation_index", "chain_pair"],
    )
    counts = (
        data.groupby(["system", "generation_index", "chain_pair"])
        .size()
        .reindex(index, fill_value=0)
        .rename("n_window_interchain_hbonds")
        .reset_index()
    )
    counts["system_label"] = counts["system"].map(LABELS)
    counts["selected_structure_only"] = True
    return counts[["system", "system_label", "generation_index", "chain_pair", "n_window_interchain_hbonds", "selected_structure_only"]]


def endpoint_summary(frames: pd.DataFrame, hb_summary: pd.DataFrame, window_hbonds: pd.DataFrame, window_dihedrals: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for system in SYSTEM_ORDER:
        for generation in [0, COMMON_MAX_GENERATION]:
            f = frames[(frames["system"] == system) & (frames["generation_index"] == generation)].iloc[0]
            h = hb_summary[(hb_summary["system"] == system) & (hb_summary["generation_index"] == generation)].iloc[0]
            wh = window_hbonds[(window_hbonds["system"] == system) & (window_hbonds["generation_index"] == generation)]
            wd = window_dihedrals[(window_dihedrals["system"] == system) & (window_dihedrals["generation_index"] == generation)]
            rows.append(
                {
                    "system": system,
                    "system_label": LABELS[system],
                    "generation_index": generation,
                    "opening_score_angstrom": f["opening_score_angstrom"],
                    "mean_triangle_area_angstrom2": f["mean_window_triangle_area_angstrom2"],
                    "window_radius_of_gyration_angstrom": f["window_ca_radius_of_gyration_angstrom"],
                    "window_rmsd_from_generation_0_angstrom": f["window_ca_rmsd_from_system_generation_0_angstrom"],
                    "mmp1_contact_pairs": f["window_mmp1_heavy_atom_contact_pairs"],
                    "minimum_mmp1_distance_angstrom": f["minimum_window_mmp1_distance_angstrom"],
                    "all_collagen_hbonds": h["n_hbonds"],
                    "all_collagen_interchain_hbonds": h["n_interchain_hbonds"],
                    "all_collagen_intrachain_hbonds": h["n_intrachain_hbonds"],
                    "window_interchain_hbonds": int(wh["n_window_interchain_hbonds"].iloc[0]),
                    "window_phi_circular_std_deg": float(wd.loc[wd["dihedral"] == "phi", "circular_std_deg"].iloc[0]),
                    "window_psi_circular_std_deg": float(wd.loc[wd["dihedral"] == "psi", "circular_std_deg"].iloc[0]),
                    "window_omega_circular_std_deg": float(wd.loc[wd["dihedral"] == "omega", "circular_std_deg"].iloc[0]),
                }
            )
    return pd.DataFrame(rows)


def style_axes(axis, ylabel: str) -> None:
    axis.set_ylabel(ylabel)
    axis.set_xlabel("Adaptive generation")
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.grid(axis="y", color="#D9D9D9", linewidth=0.6, alpha=0.7)


def plot_generation_metrics(combined: pd.DataFrame, output: Path) -> None:
    panels = [
        ("opening_score_angstrom", "Opening score (Å)", "A)"),
        ("n_window_interchain_hbonds", "Interchain hydrogen bonds\nwithin residues 977-987", "B)"),
        ("mean_window_triangle_area_angstrom2", "Mean Cα triangle area (Å²)", "C)"),
        ("window_mmp1_heavy_atom_contact_pairs", "MMP-1 heavy-atom contact pairs", "D)"),
    ]
    # Match Figure 5's effective typography when reproduced at manuscript width.
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10.0})
    fig, axes = plt.subplots(2, 2, figsize=(8.6, 6.45), constrained_layout=True)
    for axis, (metric, ylabel, letter) in zip(axes.flat, panels):
        for system in SYSTEM_ORDER:
            group = combined[combined["system"] == system].sort_values("generation_index")
            axis.plot(
                group["generation_index"],
                group[metric],
                color=COLORS[system],
                marker="o",
                markersize=3.2,
                linewidth=1.3,
                label=LABELS[system],
            )
        style_axes(axis, ylabel)
        axis.text(-0.12, 1.04, letter, transform=axis.transAxes, fontweight="bold", fontsize=12)
    axes[0, 0].legend(frameon=False, ncol=2, fontsize=10)
    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_residue_changes(geometry: pd.DataFrame, contacts: pd.DataFrame, output: Path) -> None:
    geom_ep = geometry[geometry["generation_index"].isin([0, COMMON_MAX_GENERATION])]
    geom = geom_ep.pivot_table(index=["system", "residue_position"], columns="generation_index", values="mean_ca_distance_angstrom").reset_index()
    geom["change"] = geom[COMMON_MAX_GENERATION] - geom[0]

    contact_sum = contacts.groupby(["system", "generation_index", "residue_position"], as_index=False)["heavy_atom_contact_pairs"].sum()
    contact_ep = contact_sum[contact_sum["generation_index"].isin([0, COMMON_MAX_GENERATION])]
    con = contact_ep.pivot_table(index=["system", "residue_position"], columns="generation_index", values="heavy_atom_contact_pairs").reset_index()
    con["change"] = con[COMMON_MAX_GENERATION] - con[0]

    # Use the same physical plot proportions and base type size as Figure 5.
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8.5})
    fig, axes = plt.subplots(2, 1, figsize=(7.2, 5.7), sharex=True, constrained_layout=True)
    for system in SYSTEM_ORDER:
        g = geom[geom["system"] == system]
        c = con[con["system"] == system]
        axes[0].plot(g["residue_position"], g["change"], marker="o", linewidth=1.4, color=COLORS[system], label=LABELS[system])
        axes[1].plot(c["residue_position"], c["change"], marker="o", linewidth=1.4, color=COLORS[system], label=LABELS[system])
    axes[0].axhline(0, color="#777777", linewidth=0.8)
    axes[1].axhline(0, color="#777777", linewidth=0.8)
    axes[0].set_ylabel("Change in mean Cα separation (Å)")
    axes[1].set_ylabel("Change in MMP-1 contact pairs")
    axes[1].set_xlabel("Collagen residue")
    axes[0].spines[["top", "right"]].set_visible(False)
    axes[1].spines[["top", "right"]].set_visible(False)
    axes[0].grid(axis="y", color="#D9D9D9", linewidth=0.6, alpha=0.7)
    axes[1].grid(axis="y", color="#D9D9D9", linewidth=0.6, alpha=0.7)
    axes[0].legend(frameon=False, ncol=4, fontsize=8.5)
    axes[0].text(-0.09, 1.04, "A)", transform=axes[0].transAxes, fontweight="bold", fontsize=11)
    axes[1].text(-0.09, 1.04, "B)", transform=axes[1].transAxes, fontweight="bold", fontsize=11)
    axes[1].set_xticks(range(977, 988))
    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    data = args.data_dir.resolve()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)

    hb_dir = data / "hydrogen_bonds_and_backbone"
    context_dir = data / "opening_and_mmp1_context"
    hb_summary = pd.read_csv(hb_dir / "structural_generation_summary.csv")
    hbonds = pd.read_csv(hb_dir / "structural_hydrogen_bonds.csv")
    dihedrals = pd.read_csv(hb_dir / "structural_backbone_dihedrals.csv")
    frames = pd.read_csv(context_dir / "structural_frame_summary.csv")
    geometry = pd.read_csv(context_dir / "residue_resolved_opening_geometry.csv")
    contacts = pd.read_csv(context_dir / "residue_resolved_mmp1_contacts.csv")

    tables = [hb_summary, hbonds, dihedrals, frames, geometry, contacts]
    for table in tables:
        table.drop(table[table["generation_index"] > COMMON_MAX_GENERATION].index, inplace=True)

    mapping = find_window_mapping(dihedrals)
    dihedrals["residue_position"] = [mapping.get((row.system, row.chain, int(row.resid))) for row in dihedrals.itertuples()]
    window_dih = dihedrals[dihedrals["residue_position"].notna()].copy()
    window_dihedral_summary = (
        window_dih.groupby(["system", "generation_index", "dihedral"])["angle_deg"]
        .apply(circular_std_deg)
        .rename("circular_std_deg")
        .reset_index()
    )

    hbonds["donor_position"] = [mapping.get((row.system, row.donor_chain, int(row.donor_resid))) for row in hbonds.itertuples()]
    hbonds["acceptor_position"] = [mapping.get((row.system, row.acceptor_chain, int(row.acceptor_resid))) for row in hbonds.itertuples()]
    within_window = hbonds[
        (hbonds["bond_context"] == "interchain")
        & hbonds["donor_position"].notna()
        & hbonds["acceptor_position"].notna()
    ]
    index = pd.MultiIndex.from_product(
        [SYSTEM_ORDER, range(COMMON_MAX_GENERATION + 1)], names=["system", "generation_index"]
    )
    window_hbond_counts = (
        within_window.groupby(["system", "generation_index"]).size().reindex(index, fill_value=0).rename("n_window_interchain_hbonds").reset_index()
    )

    combined = (
        frames.merge(
            hb_summary[
                [
                    "system",
                    "generation_index",
                    "n_hbonds",
                    "n_interchain_hbonds",
                    "n_intrachain_hbonds",
                    "mean_hbond_da_angstrom",
                    "phi_circular_std_deg",
                    "psi_circular_std_deg",
                    "omega_circular_std_deg",
                ]
            ],
            on=["system", "generation_index"],
            how="left",
        )
        .merge(window_hbond_counts, on=["system", "generation_index"], how="left")
    )

    endpoint = endpoint_summary(frames, hb_summary, window_hbond_counts, window_dihedral_summary)
    metrics = [
        "mean_window_triangle_area_angstrom2",
        "window_ca_radius_of_gyration_angstrom",
        "window_ca_rmsd_from_system_generation_0_angstrom",
        "window_mmp1_heavy_atom_contact_pairs",
        "minimum_window_mmp1_distance_angstrom",
        "n_hbonds",
        "n_interchain_hbonds",
        "n_intrachain_hbonds",
        "n_window_interchain_hbonds",
    ]
    early_late = early_late_summary(combined, metrics)

    geometry_endpoints = geometry[geometry["generation_index"].isin([0, COMMON_MAX_GENERATION])].copy()
    contact_endpoints = contacts[contacts["generation_index"].isin([0, COMMON_MAX_GENERATION])].copy()
    chain_pair_opening = selected_frame_chain_pair_opening(geometry)
    chain_pair_hbonds = selected_frame_chain_pair_hbonds(within_window)
    chain_pair_opening_endpoints = chain_pair_opening[
        chain_pair_opening["generation_index"].isin([0, COMMON_MAX_GENERATION])
    ].copy()

    chain_pair_hbond_early_late_rows = []
    for (system, chain_pair), group in chain_pair_hbonds.groupby(["system", "chain_pair"]):
        early = group[group["generation_index"].between(0, 4)]["n_window_interchain_hbonds"].mean()
        late = group[group["generation_index"].between(24, 28)]["n_window_interchain_hbonds"].mean()
        chain_pair_hbond_early_late_rows.append(
            {
                "system": system,
                "system_label": LABELS[system],
                "chain_pair": chain_pair,
                "early_mean_generations_0_4": float(early),
                "late_mean_generations_24_28": float(late),
                "late_minus_early": float(late - early),
                "selected_structure_only": True,
            }
        )
    chain_pair_hbond_early_late = pd.DataFrame(chain_pair_hbond_early_late_rows)

    endpoint.to_csv(output / "endpoint_structural_summary.csv", index=False)
    early_late.to_csv(output / "early_late_structural_summary.csv", index=False)
    window_hbond_counts.to_csv(output / "window_interchain_hbond_counts.csv", index=False)
    window_dihedral_summary.to_csv(output / "window_dihedral_dispersion.csv", index=False)
    geometry_endpoints.to_csv(output / "residue_geometry_endpoints.csv", index=False)
    contact_endpoints.to_csv(output / "residue_contact_endpoints.csv", index=False)
    chain_pair_opening.to_csv(output / "chain_pair_opening_selected_frames.csv", index=False)
    chain_pair_opening_endpoints.to_csv(output / "chain_pair_opening_endpoints.csv", index=False)
    chain_pair_hbonds.to_csv(output / "chain_pair_hbonds_selected_frames.csv", index=False)
    chain_pair_hbond_early_late.to_csv(output / "chain_pair_hbond_early_late.csv", index=False)

    plot_generation_metrics(combined, output / "figure_structural_metrics_by_generation.png")
    plot_residue_changes(geometry, contacts, output / "figure_residue_resolved_structural_changes.png")

    metadata = {
        "common_generation_range": [0, COMMON_MAX_GENERATION],
        "systems": SYSTEM_ORDER,
        "records": {
            "generation_summary": len(hb_summary),
            "hydrogen_bonds": len(hbonds),
            "backbone_dihedrals": len(dihedrals),
            "structural_frames": len(frames),
            "opening_geometry": len(geometry),
            "mmp1_contacts": len(contacts),
        },
        "window_mapping": {f"{system}:{chain}:{resid}": position for (system, chain, resid), position in mapping.items()},
    }
    (output / "analysis_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    print(endpoint.to_string(index=False))
    print(f"\nOutputs written to {output}")


if __name__ == "__main__":
    main()
