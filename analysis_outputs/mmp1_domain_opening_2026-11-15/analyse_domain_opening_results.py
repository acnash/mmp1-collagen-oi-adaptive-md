#!/usr/bin/env python3
"""Analyse completed MMP-1 domain-opening outputs without altering source data."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

try:
    from scipy import stats
except ImportError:  # Keep the portable plotting analysis usable without SciPy.
    from types import SimpleNamespace

    class _StatsFallback:
        @staticmethod
        def linregress(x: np.ndarray, y: np.ndarray) -> SimpleNamespace:
            slope, intercept = np.polyfit(np.asarray(x, dtype=float), np.asarray(y, dtype=float), 1)
            return SimpleNamespace(slope=float(slope), intercept=float(intercept))

        @staticmethod
        def pearsonr(x: np.ndarray, y: np.ndarray) -> SimpleNamespace:
            coefficient = np.corrcoef(np.asarray(x, dtype=float), np.asarray(y, dtype=float))[0, 1]
            return SimpleNamespace(statistic=float(coefficient))

        @staticmethod
        def spearmanr(x: np.ndarray, y: np.ndarray) -> SimpleNamespace:
            x_rank = pd.Series(np.asarray(x, dtype=float)).rank(method="average").to_numpy()
            y_rank = pd.Series(np.asarray(y, dtype=float)).rank(method="average").to_numpy()
            coefficient = np.corrcoef(x_rank, y_rank)[0, 1]
            return SimpleNamespace(statistic=float(coefficient))

    stats = _StatsFallback()


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "incoming_data" / "mmp1_domain_opening_2026-11-15" / "extracted"
OUT = Path(__file__).resolve().parent

SYSTEMS = ["wild_type", "G978S", "G984C", "G987R"]
MUTANTS = SYSTEMS[1:]
LABELS = {"wild_type": "Wild type", "G978S": "G978S", "G984C": "G984C", "G987R": "G987R"}
# Manuscript-wide system palette, established in Figure 3.
COLOURS = {"wild_type": "#303030", "G978S": "#2574A9", "G984C": "#2E8B57", "G987R": "#C0392B"}
METRICS = {
    "catalytic_hpx_com_nm": "Catalytic–HPX centre-of-mass distance",
    "ser123_ser347_ca_nm": "SER123–SER347 Cα distance",
    "catalytic_opening_n171_t230_ca_nm": "N171–T230 catalytic-opening distance",
}


def moving_block_bootstrap_ci(values: np.ndarray, block: int = 5, draws: int = 20000) -> tuple[float, float]:
    """Return a circular moving-block bootstrap CI for the sequence mean."""
    values = np.asarray(values, dtype=float)
    n = values.size
    rng = np.random.default_rng(20261115)
    starts = np.arange(n)
    means = np.empty(draws)
    blocks_needed = int(np.ceil(n / block))
    offsets = np.arange(block)
    for i in range(draws):
        chosen = rng.choice(starts, size=blocks_needed, replace=True)
        indices = ((chosen[:, None] + offsets[None, :]) % n).ravel()[:n]
        means[i] = values[indices].mean()
    return tuple(np.quantile(means, [0.025, 0.975]))


def load_data() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict]:
    generation = pd.read_csv(SOURCE / "system_generation_summary.csv")
    overall = pd.read_csv(SOURCE / "overall_system_summary.csv")
    workers = pd.read_csv(SOURCE / "worker_summary.csv")
    manifest = json.loads((SOURCE / "analysis_manifest.json").read_text(encoding="utf-8"))
    return generation, overall, workers, manifest


def validate(generation: pd.DataFrame, overall: pd.DataFrame, workers: pd.DataFrame, manifest: dict) -> list[str]:
    notes: list[str] = []
    assert set(generation.system) == set(SYSTEMS)
    assert set(overall.system) == set(SYSTEMS)
    assert set(workers.system) == set(SYSTEMS)
    assert sorted(generation.generation_index.unique()) == list(range(30))
    assert len(generation) == 120
    assert len(workers) == 1200
    for metric in METRICS:
        ncol = metric + "_n"
        bad = generation[generation[ncol] != 1000][["system", "generation_index", ncol]]
        for row in bad.itertuples(index=False):
            notes.append(f"{metric}: {row.system} generation {row.generation_index} has {getattr(row, ncol)} frames, expected 1000")
    notes.append(f"Manifest created: {manifest.get('created_utc')}")
    notes.append(f"Repository commit: {manifest.get('repository_git_commit')}")
    notes.append(f"Frame stride: {manifest.get('frame_stride')}")
    return notes


def calculate_effects(generation: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    delta_rows = []
    effect_rows = []
    for metric, metric_label in METRICS.items():
        mean_col = metric + "_mean"
        wide = generation.pivot(index="generation_index", columns="system", values=mean_col).sort_index()
        for mutant in MUTANTS:
            delta_nm = wide[mutant] - wide["wild_type"]
            delta_a = 10.0 * delta_nm
            wild_mean = float(wide["wild_type"].mean())
            ci_low, ci_high = moving_block_bootstrap_ci(delta_a.to_numpy())
            slope = stats.linregress(delta_a.index.to_numpy(float), delta_a.to_numpy(float))
            rho = stats.spearmanr(delta_a.index.to_numpy(float), delta_a.to_numpy(float))
            for generation_index, value in delta_a.items():
                delta_rows.append({
                    "generation_index": int(generation_index),
                    "mutant": mutant,
                    "metric": metric,
                    "metric_label": metric_label,
                    "mutant_minus_wild_type_angstrom": float(value),
                    "rolling_mean_5_generations_angstrom": np.nan,
                })
            effect_rows.append({
                "metric": metric,
                "metric_label": metric_label,
                "mutant": mutant,
                "generations": int(delta_a.size),
                "mean_delta_angstrom": float(delta_a.mean()),
                "median_delta_angstrom": float(delta_a.median()),
                "sd_across_generations_angstrom": float(delta_a.std(ddof=1)),
                "moving_block_bootstrap_95ci_low_angstrom": float(ci_low),
                "moving_block_bootstrap_95ci_high_angstrom": float(ci_high),
                "positive_generations_fraction": float((delta_a > 0).mean()),
                "first_5_generations_mean_delta_angstrom": float(delta_a.iloc[:5].mean()),
                "last_5_generations_mean_delta_angstrom": float(delta_a.iloc[-5:].mean()),
                "linear_delta_slope_angstrom_per_generation": float(slope.slope),
                "spearman_rho_delta_vs_generation": float(rho.statistic),
                "percent_difference_from_generation_matched_wild_type": float(100.0 * delta_nm.mean() / wild_mean),
            })

    deltas = pd.DataFrame(delta_rows)
    deltas["rolling_mean_5_generations_angstrom"] = (
        deltas.groupby(["metric", "mutant"], sort=False)["mutant_minus_wild_type_angstrom"]
        .transform(lambda s: s.rolling(5, center=True, min_periods=3).mean())
    )
    effects = pd.DataFrame(effect_rows)
    return deltas, effects


def calculate_delta_correlations(deltas: pd.DataFrame) -> pd.DataFrame:
    rows = []
    metric_pairs = [
        ("catalytic_hpx_com_nm", "ser123_ser347_ca_nm"),
        ("catalytic_hpx_com_nm", "catalytic_opening_n171_t230_ca_nm"),
        ("ser123_ser347_ca_nm", "catalytic_opening_n171_t230_ca_nm"),
    ]
    for mutant in MUTANTS:
        wide = (
            deltas[deltas.mutant == mutant]
            .pivot(index="generation_index", columns="metric", values="mutant_minus_wild_type_angstrom")
            .sort_index()
        )
        for metric_a, metric_b in metric_pairs:
            pearson = stats.pearsonr(wide[metric_a], wide[metric_b])
            spearman = stats.spearmanr(wide[metric_a], wide[metric_b])
            rows.append({
                "mutant": mutant,
                "metric_a": metric_a,
                "metric_b": metric_b,
                "pearson_r": float(pearson.statistic),
                "spearman_rho": float(spearman.statistic),
            })
    return pd.DataFrame(rows)


def plot_generation_profiles(generation: pd.DataFrame) -> None:
    fig, axes = plt.subplots(3, 1, figsize=(12.5, 12), sharex=True)
    for ax, (metric, label) in zip(axes, METRICS.items()):
        mean_col = metric + "_mean"
        for system in SYSTEMS:
            d = generation[generation.system == system].sort_values("generation_index")
            ax.plot(
                d.generation_index,
                10.0 * d[mean_col],
                color=COLOURS[system],
                lw=1.8,
                marker="o",
                ms=3.2,
                label=LABELS[system],
            )
        ax.set_ylabel("Distance (Å)")
        ax.set_title(label, loc="left", fontsize=12, weight="bold")
        ax.grid(axis="y", color="#D9D9D9", lw=0.7)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].legend(frameon=False, ncol=4, loc="upper right")
    axes[-1].set_xlabel("Adaptive generation")
    axes[-1].set_xticks(np.arange(0, 30, 2))
    fig.suptitle("MMP-1 distances across adaptive generations", fontsize=16, x=0.08, ha="left")
    fig.text(0.08, 0.015, "Lines show generation means across sampled worker frames.", fontsize=9, color="#444444")
    fig.tight_layout(rect=(0.04, 0.035, 1, 0.96))
    fig.savefig(OUT / "generation_distance_profiles.png", dpi=300, facecolor="white", bbox_inches="tight")
    plt.close(fig)


def plot_differences(deltas: pd.DataFrame) -> None:
    fig, axes = plt.subplots(3, 1, figsize=(12.5, 12), sharex=True)
    for ax, (metric, label) in zip(axes, METRICS.items()):
        dmetric = deltas[deltas.metric == metric]
        ax.axhline(0, color="#666666", lw=1)
        for mutant in MUTANTS:
            d = dmetric[dmetric.mutant == mutant].sort_values("generation_index")
            ax.plot(d.generation_index, d.mutant_minus_wild_type_angstrom, color=COLOURS[mutant], alpha=0.35, lw=1)
            ax.plot(d.generation_index, d.rolling_mean_5_generations_angstrom, color=COLOURS[mutant], lw=2.4, label=mutant)
        ax.set_ylabel("Mutant − wild type (Å)")
        ax.set_title(label, loc="left", fontsize=12, weight="bold")
        ax.grid(axis="y", color="#D9D9D9", lw=0.7)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].legend(frameon=False, ncol=3, loc="upper right")
    axes[-1].set_xlabel("Adaptive generation")
    axes[-1].set_xticks(np.arange(0, 30, 2))
    fig.suptitle("Generation-matched mutation effects", fontsize=16, x=0.08, ha="left")
    fig.text(0.08, 0.015, "Faint lines show individual generation contrasts; bold lines show centred five-generation rolling means.", fontsize=9, color="#444444")
    fig.tight_layout(rect=(0.04, 0.035, 1, 0.96))
    fig.savefig(OUT / "mutant_minus_wild_type_by_generation.png", dpi=300, facecolor="white", bbox_inches="tight")
    plt.close(fig)


def plot_effect_summary(effects: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(15, 5.6), sharey=True)
    y = np.arange(len(MUTANTS))[::-1]
    for ax, (metric, label) in zip(axes, METRICS.items()):
        d = effects[effects.metric == metric].set_index("mutant").loc[MUTANTS]
        means = d.mean_delta_angstrom.to_numpy()
        lows = d.moving_block_bootstrap_95ci_low_angstrom.to_numpy()
        highs = d.moving_block_bootstrap_95ci_high_angstrom.to_numpy()
        ax.axvline(0, color="#666666", lw=1)
        for i, mutant in enumerate(MUTANTS):
            yy = y[i]
            ax.errorbar(
                means[i], yy,
                xerr=[[means[i] - lows[i]], [highs[i] - means[i]]],
                fmt="o", color=COLOURS[mutant], ms=7, capsize=4, lw=1.7,
            )
            ax.text(means[i], yy + 0.19, f"{means[i]:+.2f} Å", ha="center", va="bottom", fontsize=9)
        ax.set_title(label, fontsize=11, weight="bold")
        ax.set_xlabel("Mean mutant − wild type (Å)")
        ax.grid(axis="x", color="#D9D9D9", lw=0.7)
        ax.spines[["top", "right", "left"]].set_visible(False)
        ax.set_yticks(y, [LABELS[m] for m in MUTANTS])
    fig.suptitle("Mutation effects across 30 generation-matched contrasts", fontsize=16, x=0.04, ha="left")
    fig.text(0.04, 0.015, "Intervals are 95% circular moving-block bootstrap intervals for the mean contrast (block length 5 generations).", fontsize=9, color="#444444")
    fig.tight_layout(rect=(0.02, 0.05, 1, 0.93))
    fig.savefig(OUT / "mutation_effect_summary.png", dpi=300, facecolor="white", bbox_inches="tight")
    plt.close(fig)


def write_findings(effects: pd.DataFrame, notes: list[str]) -> None:
    ordered = effects.sort_values(["metric", "mean_delta_angstrom"])
    lines = ["Quality-control notes", *[f"- {n}" for n in notes], "", "Generation-level effects"]
    for row in ordered.itertuples(index=False):
        lines.append(
            f"- {row.metric_label}; {row.mutant}: mean {row.mean_delta_angstrom:+.3f} Å, "
            f"block-bootstrap 95% CI [{row.moving_block_bootstrap_95ci_low_angstrom:+.3f}, "
            f"{row.moving_block_bootstrap_95ci_high_angstrom:+.3f}], "
            f"positive in {100*row.positive_generations_fraction:.1f}% of generations, "
            f"first five {row.first_5_generations_mean_delta_angstrom:+.3f} Å, "
            f"last five {row.last_5_generations_mean_delta_angstrom:+.3f} Å."
        )
    (OUT / "analysis_findings.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    generation, overall, workers, manifest = load_data()
    notes = validate(generation, overall, workers, manifest)
    deltas, effects = calculate_effects(generation)
    correlations = calculate_delta_correlations(deltas)
    deltas.to_csv(OUT / "generation_matched_mutant_differences.csv", index=False)
    effects.to_csv(OUT / "generation_level_effect_summary.csv", index=False)
    correlations.to_csv(OUT / "metric_delta_correlations.csv", index=False)
    plot_generation_profiles(generation)
    plot_differences(deltas)
    plot_effect_summary(effects)
    write_findings(effects, notes)
    print(effects.to_string(index=False))
    print("\nMetric-delta correlations\n" + correlations.to_string(index=False))
    print("\n" + "\n".join(notes))


if __name__ == "__main__":
    main()
