#!/usr/bin/env python3
"""Analyse the compact and frame-level outputs from the joint unwinding report.

The script treats adaptive generations as a directed sampling campaign. It therefore
reports descriptive changes and worker variability without assigning equilibrium
probabilities or treating generations as statistically independent replicates.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


SYSTEM_ORDER = ["wild_type", "G978S", "G984C", "G987R"]
LABELS = {
    "wild_type": "Wild type",
    "G978S": "G978S",
    "G984C": "G984C",
    "G987R": "G987R",
}
COLORS = {
    "wild_type": "#222222",
    "G978S": "#2166ac",
    "G984C": "#1b9e77",
    "G987R": "#b2182b",
}


def safe_sd(values: pd.Series) -> float:
    return float(values.std(ddof=1)) if len(values) > 1 else float("nan")


def write_compact_analyses(input_dir: Path, output_dir: Path) -> dict:
    best = pd.read_csv(input_dir / "best_worker_unwinding_by_generation.csv")
    generations = pd.read_csv(input_dir / "joint_generation_unwinding_summary.csv")
    workers = pd.read_csv(input_dir / "joint_worker_unwinding_summary.csv")

    complete_generations = generations[generations["generation_status"] == "completed"].copy()
    last_completed = complete_generations.groupby("system")["generation_index"].max()
    last_common = int(last_completed.reindex(SYSTEM_ORDER).min())

    best_long = best.melt(
        id_vars="generation",
        value_vars=SYSTEM_ORDER,
        var_name="system",
        value_name="best_score_angstrom",
    ).dropna()
    best_long["system_label"] = best_long["system"].map(LABELS)
    wt = best[["generation", "wild_type"]].rename(columns={"wild_type": "wild_type_score_angstrom"})
    best_long = best_long.merge(wt, on="generation", how="left")
    best_long["difference_from_wild_type_angstrom"] = (
        best_long["best_score_angstrom"] - best_long["wild_type_score_angstrom"]
    )
    best_long.to_csv(output_dir / "matched_generation_best_scores.csv", index=False)

    completed_workers = workers[workers["status"].isin(["completed", "completed_approximate"])].copy()
    worker_rows = []
    for (system, generation), frame in completed_workers.groupby(["system", "generation_index"], sort=True):
        values = frame["max_opening_score_angstrom"].dropna()
        worker_rows.append(
            {
                "system": system,
                "system_label": LABELS.get(system, system),
                "generation_index": int(generation),
                "n_workers": int(len(values)),
                "mean_worker_max_angstrom": float(values.mean()),
                "sd_worker_max_angstrom": safe_sd(values),
                "median_worker_max_angstrom": float(values.median()),
                "q25_worker_max_angstrom": float(values.quantile(0.25)),
                "q75_worker_max_angstrom": float(values.quantile(0.75)),
                "maximum_worker_max_angstrom": float(values.max()),
            }
        )
    worker_stats = pd.DataFrame(worker_rows)
    worker_stats.to_csv(output_dir / "worker_maximum_descriptive_statistics.csv", index=False)

    progression_rows = []
    threshold_rows = []
    for system in SYSTEM_ORDER:
        subset = best.loc[best["generation"] <= last_common, ["generation", system]].dropna()
        initial = float(subset.loc[subset["generation"] == 0, system].iloc[0])
        endpoint = float(subset.loc[subset["generation"] == last_common, system].iloc[0])
        slope = float(np.polyfit(subset["generation"], subset[system], 1)[0])
        endpoint_workers = completed_workers[
            (completed_workers["system"] == system)
            & (completed_workers["generation_index"] == last_common)
        ]["max_opening_score_angstrom"].dropna()
        progression_rows.append(
            {
                "system": system,
                "system_label": LABELS[system],
                "last_common_completed_generation": last_common,
                "initial_best_score_angstrom": initial,
                "endpoint_best_score_angstrom": endpoint,
                "absolute_change_angstrom": endpoint - initial,
                "relative_change_percent": 100.0 * (endpoint / initial - 1.0),
                "linear_descriptive_slope_angstrom_per_generation": slope,
                "endpoint_n_completed_workers": int(len(endpoint_workers)),
                "endpoint_worker_median_max_angstrom": float(endpoint_workers.median()),
                "endpoint_worker_q25_max_angstrom": float(endpoint_workers.quantile(0.25)),
                "endpoint_worker_q75_max_angstrom": float(endpoint_workers.quantile(0.75)),
                "endpoint_worker_mean_max_angstrom": float(endpoint_workers.mean()),
                "endpoint_worker_sd_max_angstrom": safe_sd(endpoint_workers),
            }
        )
        available = best[["generation", system]].dropna()
        for threshold in [7.0, 8.0, 9.0, 10.0]:
            crossed = available.loc[available[system] >= threshold, "generation"]
            threshold_rows.append(
                {
                    "system": system,
                    "system_label": LABELS[system],
                    "threshold_angstrom": threshold,
                    "first_generation_at_or_above_threshold": (
                        int(crossed.iloc[0]) if not crossed.empty else pd.NA
                    ),
                }
            )

    progression = pd.DataFrame(progression_rows)
    thresholds = pd.DataFrame(threshold_rows)
    progression.to_csv(output_dir / "campaign_progression_summary.csv", index=False)
    thresholds.to_csv(output_dir / "opening_threshold_crossings.csv", index=False)

    overview = (
        generations.groupby("system", as_index=False)
        .agg(
            first_generation=("generation_index", "min"),
            latest_generation=("generation_index", "max"),
            completed_generations=("generation_status", lambda x: int((x == "completed").sum())),
            running_generations=("generation_status", lambda x: int((x == "running").sum())),
        )
    )
    overview["system_label"] = overview["system"].map(LABELS)
    overview.to_csv(output_dir / "campaign_overview.csv", index=False)

    return {
        "best": best,
        "workers": completed_workers,
        "progression": progression,
        "last_common": last_common,
    }


def stream_frame_statistics(input_dir: Path, output_dir: Path, chunksize: int) -> pd.DataFrame:
    source = input_dir / "joint_opening_score_frames.csv"
    usecols = ["system", "generation_index", "opening_score_angstrom", "status"]
    partials = []
    for chunk in pd.read_csv(source, usecols=usecols, chunksize=chunksize):
        chunk = chunk[chunk["status"].isin(["completed", "completed_approximate"])]
        grouped = chunk.groupby(["system", "generation_index"])["opening_score_angstrom"]
        part = grouped.agg(["count", "sum", "min", "max"])
        part["sum_sq"] = grouped.apply(lambda x: float(np.square(x.to_numpy(dtype=float)).sum()))
        partials.append(part.reset_index())
    merged = pd.concat(partials, ignore_index=True)
    rows = []
    for (system, generation), frame in merged.groupby(["system", "generation_index"], sort=True):
        count = int(frame["count"].sum())
        total = float(frame["sum"].sum())
        total_sq = float(frame["sum_sq"].sum())
        mean = total / count
        variance = max(0.0, (total_sq - count * mean * mean) / (count - 1)) if count > 1 else np.nan
        rows.append(
            {
                "system": system,
                "system_label": LABELS.get(system, system),
                "generation_index": int(generation),
                "n_frames": count,
                "mean_opening_score_angstrom": mean,
                "sd_opening_score_angstrom": float(np.sqrt(variance)),
                "minimum_opening_score_angstrom": float(frame["min"].min()),
                "maximum_opening_score_angstrom": float(frame["max"].max()),
            }
        )
    result = pd.DataFrame(rows)
    result.to_csv(output_dir / "frame_level_generation_statistics.csv", index=False)
    return result


def make_figures(results: dict, output_dir: Path) -> None:
    best = results["best"]
    last_common = results["last_common"]
    progression = results["progression"]
    workers = results["workers"]

    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    for system in SYSTEM_ORDER:
        ax.plot(
            best["generation"],
            best[system],
            marker="o",
            markersize=3.2,
            linewidth=1.6,
            color=COLORS[system],
            label=LABELS[system],
        )
    ax.axvline(last_common, color="#777777", linestyle="--", linewidth=1.0)
    ax.text(last_common - 0.4, ax.get_ylim()[0], " last common completed generation", rotation=90,
            va="bottom", ha="right", fontsize=8, color="#555555")
    ax.set_xlabel("Adaptive generation")
    ax.set_ylabel("Best opening score (Å)")
    ax.legend(frameon=False, ncol=2)
    ax.grid(alpha=0.18)
    fig.tight_layout()
    fig.savefig(output_dir / "figure_adaptive_unwinding_progression.png", dpi=300)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.2, 4.4))
    values = [
        workers[(workers["system"] == system) & (workers["generation_index"] == last_common)][
            "max_opening_score_angstrom"
        ].dropna().to_numpy()
        for system in SYSTEM_ORDER
    ]
    bp = ax.boxplot(
        values,
        tick_labels=[LABELS[s] for s in SYSTEM_ORDER],
        patch_artist=True,
        showfliers=True,
    )
    for patch, system in zip(bp["boxes"], SYSTEM_ORDER):
        patch.set_facecolor(COLORS[system])
        patch.set_alpha(0.75)
    ax.set_ylabel("Worker maximum opening score (Å)")
    ax.set_title(f"Generation {last_common}")
    ax.grid(axis="y", alpha=0.18)
    fig.tight_layout()
    fig.savefig(output_dir / "figure_worker_variability_common_endpoint.png", dpi=300)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.2, 4.4))
    bars = ax.barh(
        progression["system_label"],
        progression["absolute_change_angstrom"],
        color=[COLORS[s] for s in progression["system"]],
    )
    ax.axvline(0, color="#333333", linewidth=0.8)
    ax.set_xlabel(f"Change through generation {last_common} (Å)")
    ax.bar_label(bars, fmt="%.2f", padding=3, fontsize=9)
    ax.grid(axis="x", alpha=0.18)
    fig.tight_layout()
    fig.savefig(output_dir / "figure_change_from_generation_zero.png", dpi=300)
    plt.close(fig)


def write_text_summary(results: dict, output_dir: Path) -> None:
    last_common = results["last_common"]
    progression = results["progression"].set_index("system")
    lines = [
        "Adaptive MMP-1 collagen unwinding analysis",
        "",
        f"The last generation completed by all four systems was generation {last_common}.",
        "Generation values are adaptively selected and are not independent equilibrium samples.",
        "",
    ]
    for system in SYSTEM_ORDER:
        row = progression.loc[system]
        lines.append(
            f"{LABELS[system]}: {row.initial_best_score_angstrom:.3f} Å at generation 0; "
            f"{row.endpoint_best_score_angstrom:.3f} Å at generation {last_common}; "
            f"change {row.absolute_change_angstrom:+.3f} Å "
            f"({row.relative_change_percent:+.1f}%)."
        )
    (output_dir / "analysis_summary.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_dir", type=Path, help="Extracted joint_unwinding_report directory")
    parser.add_argument("--output_dir", type=Path, default=Path("analysis_results/unwinding_report"))
    parser.add_argument("--chunksize", type=int, default=250_000)
    parser.add_argument("--skip_frame_table", action="store_true")
    args = parser.parse_args()

    input_dir = args.input_dir.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    results = write_compact_analyses(input_dir, output_dir)
    if not args.skip_frame_table:
        stream_frame_statistics(input_dir, output_dir, args.chunksize)
    make_figures(results, output_dir)
    write_text_summary(results, output_dir)
    print(f"Analysis written to {output_dir}")


if __name__ == "__main__":
    main()
