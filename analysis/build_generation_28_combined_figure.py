#!/usr/bin/env python3
"""Combine adaptive progression and endpoint worker variability in one figure."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


SYSTEM_ORDER = ["wild_type", "G978S", "G984C", "G987R"]
LABELS = {"wild_type": "Wild type", "G978S": "G978S", "G984C": "G984C", "G987R": "G987R"}
COLORS = {"wild_type": "#4D4D4D", "G978S": "#2878B5", "G984C": "#D05A2A", "G987R": "#7A4EA3"}
ENDPOINT = 28


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--selected-frames", type=Path, required=True)
    parser.add_argument("--worker-maxima", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    selected = pd.read_csv(args.selected_frames)
    selected = selected.rename(
        columns={"generation": "generation_index", "best_score_angstrom": "opening_score_angstrom"}
    )
    selected = selected[selected["generation_index"].between(0, ENDPOINT)]
    workers = pd.read_csv(args.worker_maxima)
    workers = workers.rename(columns={"max_opening_score_angstrom": "maximum_opening_score_angstrom"})
    workers = workers[workers["generation_index"] == ENDPOINT]

    # Figure 5 establishes the manuscript's effective printed type size.  Figure
    # 1 is wider, so its source type must be larger to remain comparable after
    # both figures are scaled to the same manuscript width.
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 12.0})
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.35), gridspec_kw={"width_ratios": [1.55, 1.0]})
    for system in SYSTEM_ORDER:
        group = selected[selected["system"] == system].sort_values("generation_index")
        axes[0].plot(
            group["generation_index"],
            group["opening_score_angstrom"],
            color=COLORS[system],
            marker="o",
            markersize=3.3,
            linewidth=1.7,
            label=LABELS[system],
        )
    axes[0].set_xlim(-0.5, ENDPOINT + 0.5)
    axes[0].set_xlabel("Adaptive generation")
    axes[0].set_ylabel("Selected opening score (Å)")
    axes[0].legend(frameon=False, ncol=2, loc="upper left", fontsize=12.0)

    values = [
        workers.loc[workers["system"] == system, "maximum_opening_score_angstrom"].dropna().to_numpy()
        for system in SYSTEM_ORDER
    ]
    boxes = axes[1].boxplot(
        values,
        tick_labels=[LABELS[system] for system in SYSTEM_ORDER],
        patch_artist=True,
        showfliers=True,
        widths=0.58,
        medianprops={"color": "white", "linewidth": 1.4},
    )
    for patch, system in zip(boxes["boxes"], SYSTEM_ORDER):
        patch.set_facecolor(COLORS[system])
        patch.set_edgecolor(COLORS[system])
        patch.set_alpha(0.88)
    axes[1].set_ylabel("Worker maximum opening score (Å)")
    axes[1].tick_params(axis="x", rotation=18)

    for letter, axis in zip(("A)", "B)"), axes):
        axis.text(-0.12, 1.04, letter, transform=axis.transAxes, fontweight="bold", fontsize=15)
        axis.grid(axis="y", color="#D8D8D8", linewidth=0.6)
        axis.spines[["top", "right"]].set_visible(False)

    fig.tight_layout()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=350, bbox_inches="tight", facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    main()
