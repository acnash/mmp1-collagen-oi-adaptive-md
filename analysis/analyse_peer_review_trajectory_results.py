#!/usr/bin/env python3
"""Summarise the trajectory-wide peer-review analysis for the manuscript.

The input tables are produced by run_peer_review_trajectory_analysis_standalone.py.
Each worker is treated as the descriptive unit.  No inferential tests are
performed because workers and generations belong to adaptively selected
lineages.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


SYSTEM_ORDER = ["wild_type", "G978S", "G984C", "G987R"]
SYSTEM_LABELS = {
    "wild_type": "Wild type",
    "G978S": "G978S",
    "G984C": "G984C",
    "G987R": "G987R",
}
COLOURS = {
    "wild_type": "#4D4D4D",
    "G978S": "#2A6FBB",
    "G984C": "#C95A28",
    "G987R": "#7A4EA3",
}
GENERATIONS = [0, 10, 20, 28]
KEYS = ["system", "generation_index", "worker"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--prior-input-dir", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def describe_workers(table: pd.DataFrame, metrics: list[str]) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for (system, generation), group in table.groupby(["system", "generation_index"]):
        row: dict[str, object] = {
            "system": system,
            "system_label": SYSTEM_LABELS[system],
            "generation_index": generation,
            "n_workers": group["worker"].nunique(),
        }
        for metric in metrics:
            values = group[metric].dropna()
            row[f"{metric}_median"] = values.median()
            row[f"{metric}_q25"] = values.quantile(0.25)
            row[f"{metric}_q75"] = values.quantile(0.75)
        rows.append(row)
    return pd.DataFrame(rows).sort_values(["system", "generation_index"])


def validate(frame: pd.DataFrame, worker: pd.DataFrame) -> None:
    expected = pd.MultiIndex.from_product(
        [SYSTEM_ORDER, GENERATIONS, [f"worker_{index:03d}" for index in range(10)]],
        names=KEYS,
    )
    observed = pd.MultiIndex.from_frame(frame[KEYS].drop_duplicates())
    missing = expected.difference(observed)
    extra = observed.difference(expected)
    if len(missing) or len(extra):
        raise ValueError(f"Unexpected worker coverage; missing={list(missing)}, extra={list(extra)}")
    counts = frame.groupby(KEYS).size()
    if not (counts == 100).all():
        raise ValueError("Every worker must contain exactly 100 sampled frames")
    if frame.duplicated(KEYS + ["sampled_frame_index"]).any():
        raise ValueError("Duplicate sampled frames were found")
    if frame.isna().any().any():
        raise ValueError("Missing frame-level metric values were found")
    if len(worker) != 160 or not (worker["n_sampled_frames"] == 100).all():
        raise ValueError("Worker summary coverage is incomplete")


def panel_boxes(ax, table: pd.DataFrame, metric: str, ylabel: str, letter: str) -> None:
    offsets = dict(zip(SYSTEM_ORDER, [-0.27, -0.09, 0.09, 0.27]))
    width = 0.15
    for system in SYSTEM_ORDER:
        series = []
        positions = []
        for index, generation in enumerate(GENERATIONS):
            values = table.loc[
                (table["system"] == system) & (table["generation_index"] == generation), metric
            ].to_numpy()
            series.append(values)
            positions.append(index + offsets[system])
        box = ax.boxplot(
            series,
            positions=positions,
            widths=width,
            patch_artist=True,
            showfliers=False,
            medianprops={"color": "white", "linewidth": 1.2},
            whiskerprops={"color": COLOURS[system], "linewidth": 0.8},
            capprops={"color": COLOURS[system], "linewidth": 0.8},
            boxprops={"edgecolor": COLOURS[system], "linewidth": 0.8},
        )
        for patch in box["boxes"]:
            patch.set_facecolor(COLOURS[system])
            patch.set_alpha(0.86)
        ax.plot([], [], color=COLOURS[system], linewidth=6, label=SYSTEM_LABELS[system])
    ax.set_xticks(range(len(GENERATIONS)), [str(value) for value in GENERATIONS])
    ax.set_xlabel("Adaptive generation")
    ax.set_ylabel(ylabel)
    ax.text(-0.13, 1.04, letter, transform=ax.transAxes, fontweight="bold", fontsize=11)
    ax.grid(axis="y", color="#D9D9D9", linewidth=0.6)


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    input_dirs = [args.input_dir]
    if args.prior_input_dir is not None:
        input_dirs.insert(0, args.prior_input_dir)
    frame = pd.concat(
        [pd.read_csv(path / "trajectory_frame_metrics.csv") for path in input_dirs], ignore_index=True
    )
    worker = pd.concat(
        [pd.read_csv(path / "worker_metric_summaries.csv") for path in input_dirs], ignore_index=True
    )
    dihedral = pd.concat(
        [pd.read_csv(path / "residue_chain_dihedral_summaries.csv") for path in input_dirs], ignore_index=True
    )
    frame = frame[frame["generation_index"].isin(GENERATIONS)].copy()
    worker = worker[worker["generation_index"].isin(GENERATIONS)].copy()
    dihedral = dihedral[dihedral["generation_index"].isin(GENERATIONS)].copy()
    validate(frame, worker)

    frame_metrics = [
        "opening_score_angstrom",
        "mean_ca_triangle_area_angstrom2",
        "mean_ca_distance_asymmetry_angstrom",
        "hbond_n_da_3.5A_angle_135",
        "hbond_n_da_3.5A_angle_135_backbone_only",
        "unique_collagen_mmp_residue_contacts",
        "unique_catalytic_domain_contacts",
        "unique_linker_contacts",
        "unique_hemopexin_domain_contacts",
        "glu200_oxygen_to_zinc_angstrom",
        "his199_ne2_to_zinc_angstrom",
        "his203_ne2_to_zinc_angstrom",
        "his209_ne2_to_zinc_angstrom",
    ]
    worker_means = frame.groupby(KEYS, as_index=False)[frame_metrics].mean()

    for atom in ("C", "O"):
        columns = [f"chain_{chain}_gly975_scissile_{atom}_to_zinc_angstrom" for chain in (1, 2, 3)]
        frame[f"closest_scissile_{atom}_to_zinc_angstrom"] = frame[columns].min(axis=1)
    columns = [f"chain_{chain}_leu976_scissile_N_to_zinc_angstrom" for chain in (1, 2, 3)]
    frame["closest_scissile_N_to_zinc_angstrom"] = frame[columns].min(axis=1)
    columns = [f"chain_{chain}_glu200_O_to_scissile_C_angstrom" for chain in (1, 2, 3)]
    frame["closest_glu200_O_to_scissile_C_angstrom"] = frame[columns].min(axis=1)
    catalytic_metrics = [
        "closest_scissile_C_to_zinc_angstrom",
        "closest_scissile_O_to_zinc_angstrom",
        "closest_scissile_N_to_zinc_angstrom",
        "closest_glu200_O_to_scissile_C_angstrom",
    ]
    catalytic_worker = frame.groupby(KEYS, as_index=False)[catalytic_metrics].mean()
    worker_means = worker_means.merge(catalytic_worker, on=KEYS, validate="one_to_one")

    backbone = (
        dihedral.loc[dihedral["dihedral"].isin(["phi", "psi"])]
        .groupby(KEYS + ["dihedral"], as_index=False)["circular_std_deg"]
        .mean()
        .pivot(index=KEYS, columns="dihedral", values="circular_std_deg")
        .rename(columns={"phi": "phi_circular_std_deg", "psi": "psi_circular_std_deg"})
        .reset_index()
    )
    worker_means = worker_means.merge(backbone, on=KEYS, validate="one_to_one")

    summary_metrics = [column for column in worker_means.columns if column not in KEYS]
    summary = describe_workers(worker_means, summary_metrics)
    summary.to_csv(args.output_dir / "trajectory_worker_summary_by_system_generation.csv", index=False)

    endpoints = summary.loc[summary["generation_index"].isin([0, 28])].copy()
    endpoints.to_csv(args.output_dir / "trajectory_endpoint_summary.csv", index=False)

    changes = []
    for system in SYSTEM_ORDER:
        subset = worker_means.loc[worker_means["system"] == system].set_index(["worker", "generation_index"])
        for metric in summary_metrics:
            delta = subset[metric].unstack()[28] - subset[metric].unstack()[0]
            changes.append(
                {
                    "system": system,
                    "system_label": SYSTEM_LABELS[system],
                    "metric": metric,
                    "median_change": delta.median(),
                    "q25_change": delta.quantile(0.25),
                    "q75_change": delta.quantile(0.75),
                }
            )
    pd.DataFrame(changes).to_csv(args.output_dir / "trajectory_endpoint_changes.csv", index=False)

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8.5})
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.4), constrained_layout=True)
    panel_boxes(axes[0, 0], worker_means, "opening_score_angstrom", "Mean opening score (Å)", "A)")
    panel_boxes(
        axes[0, 1],
        worker_means,
        "hbond_n_da_3.5A_angle_135",
        "Interchain hydrogen bonds",
        "B)",
    )
    panel_boxes(axes[1, 0], worker_means, "phi_circular_std_deg", "Mean φ circular SD (°)", "C)")
    panel_boxes(
        axes[1, 1],
        worker_means,
        "unique_catalytic_domain_contacts",
        "Catalytic-domain residue contacts",
        "D)",
    )
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=4, frameon=False, bbox_to_anchor=(0.5, -0.035))
    for ax in axes.flat:
        legend = ax.get_legend()
        if legend is not None:
            legend.remove()
    fig.savefig(args.output_dir / "figure_trajectory_wide_validation.png", dpi=350, bbox_inches="tight")
    plt.close(fig)

    print(f"Validated {len(worker_means)} workers and {len(frame)} sampled frames")
    print(f"Wrote trajectory-wide summaries to {args.output_dir}")


if __name__ == "__main__":
    main()
