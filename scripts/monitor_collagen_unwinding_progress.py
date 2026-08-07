#!/usr/bin/env python3
"""
Generate a joint collagen-unwinding progress report for adaptive MMP1 runs.

This monitor reads adaptive_run_* output directories and summarizes only the
PBC-corrected collagen opening score reported by the MD runner. It does not
perform structural clustering, contact analysis, hydrogen-bond analysis, or
trajectory reanalysis.
"""

import argparse
import json
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


@dataclass
class RunSpec:
    system: str
    label: str
    run_dir: Path


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
            system_dir = (root / meta["relative_dir"]).resolve()
            run_dir = find_latest_run(system_dir, meta["run_prefix"])
            if run_dir is None:
                continue
        if not run_dir.exists():
            raise SystemExit(f"{system} run directory does not exist: {run_dir}")
        specs.append(RunSpec(system=system, label=meta["label"], run_dir=run_dir))
    return specs


def find_generation_dirs(run_dir: Path) -> List[Path]:
    return sorted([p for p in run_dir.glob("generation_*") if p.is_dir()])


def find_worker_dirs(generation_dir: Path) -> List[Path]:
    return sorted([p for p in generation_dir.glob("worker_*") if p.is_dir()])


def worker_status_counts(worker_dirs: Iterable[Path]) -> Dict[str, int]:
    counts: Dict[str, int] = {}
    for worker_dir in worker_dirs:
        meta = read_json(worker_dir / "worker.json")
        status = str(meta.get("status", "unknown"))
        counts[status] = counts.get(status, 0) + 1
    return counts


def read_worker_score_table(spec: RunSpec, generation_dir: Path, worker_dir: Path) -> Tuple[Optional[pd.DataFrame], dict]:
    meta = read_json(worker_dir / "worker.json")
    generation_index = parse_generation_index(generation_dir.name)
    score_csv = worker_dir / "opening_scores.csv"
    row = {
        "system": spec.system,
        "system_label": spec.label,
        "run_dir": str(spec.run_dir),
        "generation": generation_dir.name,
        "generation_index": generation_index,
        "worker": worker_dir.name,
        "worker_dir": str(worker_dir),
        "status": meta.get("status", "unknown"),
        "is_replacement": bool(meta.get("is_replacement", False)),
        "score_file_found": score_csv.exists(),
        "n_score_frames": 0,
        "max_opening_score_angstrom": np.nan,
        "mean_opening_score_angstrom": np.nan,
        "final_opening_score_angstrom": np.nan,
        "best_time_ps": np.nan,
        "best_step": np.nan,
    }
    if not score_csv.exists():
        return None, row
    try:
        df = pd.read_csv(score_csv)
    except Exception:
        return None, row
    if df.empty or "opening_score_angstrom" not in df.columns:
        return None, row

    df = df.copy()
    if "time_ps" not in df.columns:
        df["time_ps"] = np.nan
    if "step" not in df.columns:
        df["step"] = np.arange(len(df))
    df["system"] = spec.system
    df["system_label"] = spec.label
    df["run_dir"] = str(spec.run_dir)
    df["generation"] = generation_dir.name
    df["generation_index"] = generation_index
    df["worker"] = worker_dir.name
    df["worker_dir"] = str(worker_dir)
    df["status"] = row["status"]

    best_idx = df["opening_score_angstrom"].idxmax()
    row.update(
        {
            "n_score_frames": int(len(df)),
            "max_opening_score_angstrom": float(df["opening_score_angstrom"].max()),
            "mean_opening_score_angstrom": float(df["opening_score_angstrom"].mean()),
            "final_opening_score_angstrom": float(df["opening_score_angstrom"].iloc[-1]),
            "best_time_ps": float(df.loc[best_idx, "time_ps"]),
            "best_step": int(df.loc[best_idx, "step"]),
        }
    )
    return df, row


def collect_run_tables(specs: List[RunSpec]) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    frame_tables = []
    worker_rows = []
    generation_rows = []

    for spec in specs:
        context = read_json(spec.run_dir / "context.json")
        completed = {
            int(entry.get("generation")): entry
            for entry in context.get("completed_generations", [])
            if entry.get("generation") is not None
        }
        for generation_dir in find_generation_dirs(spec.run_dir):
            generation_index = parse_generation_index(generation_dir.name)
            worker_dirs = find_worker_dirs(generation_dir)
            status_counts = worker_status_counts(worker_dirs)
            generation_worker_rows = []

            for worker_dir in worker_dirs:
                frame_df, worker_row = read_worker_score_table(spec, generation_dir, worker_dir)
                worker_rows.append(worker_row)
                generation_worker_rows.append(worker_row)
                if frame_df is not None:
                    frame_tables.append(frame_df)

            generation_meta = read_json(generation_dir / "generation.json")
            worker_df = pd.DataFrame(generation_worker_rows)
            scored = worker_df.dropna(subset=["max_opening_score_angstrom"]) if not worker_df.empty else pd.DataFrame()
            completed_entry = completed.get(generation_index, {}) if generation_index is not None else {}
            generation_rows.append(
                {
                    "system": spec.system,
                    "system_label": spec.label,
                    "run_dir": str(spec.run_dir),
                    "generation": generation_dir.name,
                    "generation_index": generation_index,
                    "generation_status": generation_meta.get("status", ""),
                    "n_workers": len(worker_dirs),
                    "n_completed": int(status_counts.get("completed", 0) + status_counts.get("completed_approximate", 0)),
                    "n_running": int(status_counts.get("running", 0)),
                    "n_queued": int(status_counts.get("queued", 0)),
                    "n_failed_archived": int(status_counts.get("failed_archived", 0)),
                    "n_workers_with_scores": int(len(scored)),
                    "best_worker_score_angstrom": np.nan if scored.empty else float(scored["max_opening_score_angstrom"].max()),
                    "mean_worker_max_score_angstrom": np.nan if scored.empty else float(scored["max_opening_score_angstrom"].mean()),
                    "selected_best_score_angstrom": generation_meta.get("best_score_angstrom", completed_entry.get("best_score_angstrom", np.nan)),
                }
            )

    frame_df = pd.concat(frame_tables, ignore_index=True) if frame_tables else pd.DataFrame()
    worker_df = pd.DataFrame(worker_rows)
    generation_df = pd.DataFrame(generation_rows)
    if not generation_df.empty:
        generation_df = generation_df.sort_values(["system", "generation_index"])
    return frame_df, worker_df, generation_df


def latest_per_system(generation_df: pd.DataFrame) -> pd.DataFrame:
    if generation_df.empty:
        return generation_df
    rows = []
    for _, sdf in generation_df.dropna(subset=["generation_index"]).groupby("system"):
        rows.append(sdf.sort_values("generation_index").iloc[-1])
    return pd.DataFrame(rows)


def best_worker_unwinding_matrix(generation_df: pd.DataFrame) -> pd.DataFrame:
    columns = ["generation", "wild_type", "G978S", "G984C", "G987R"]
    if generation_df.empty:
        return pd.DataFrame(columns=columns)
    required = {"generation_index", "system", "best_worker_score_angstrom"}
    if not required.issubset(generation_df.columns):
        return pd.DataFrame(columns=columns)

    df = generation_df.dropna(subset=["generation_index"]).copy()
    if df.empty:
        return pd.DataFrame(columns=columns)

    pivot = df.pivot_table(
        index="generation_index",
        columns="system",
        values="best_worker_score_angstrom",
        aggfunc="max",
    )
    pivot = pivot.reindex(columns=["wild_type", "G978S", "G984C", "G987R"])
    pivot = pivot.reset_index().rename(columns={"generation_index": "generation"})
    pivot["generation"] = pivot["generation"].astype(int)
    return pivot[columns]


def write_csv_outputs(output_dir: Path, frame_df: pd.DataFrame, worker_df: pd.DataFrame, generation_df: pd.DataFrame) -> None:
    safe_mkdir(output_dir)
    frame_df.to_csv(output_dir / "joint_opening_score_frames.csv", index=False)
    worker_df.to_csv(output_dir / "joint_worker_unwinding_summary.csv", index=False)
    generation_df.to_csv(output_dir / "joint_generation_unwinding_summary.csv", index=False)
    latest_per_system(generation_df).to_csv(output_dir / "latest_generation_by_system.csv", index=False)
    best_worker_unwinding_matrix(generation_df).to_csv(output_dir / "best_worker_unwinding_by_generation.csv", index=False)


def add_text_page(pdf: PdfPages, title: str, lines: List[str]) -> None:
    fig = plt.figure(figsize=(8.27, 11.69))
    ax = fig.add_subplot(111)
    ax.axis("off")
    ax.text(0.06, 0.95, title, fontsize=16, weight="bold", va="top")
    ax.text(0.06, 0.89, "\n".join(lines), fontsize=10.5, va="top", family="monospace")
    pdf.savefig(fig)
    plt.close(fig)


def plot_generation_progress(pdf: PdfPages, output_dir: Path, generation_df: pd.DataFrame) -> None:
    if generation_df.empty:
        return
    fig = plt.figure(figsize=(8.27, 5.85))
    ax = fig.add_subplot(111)
    for system, meta in SYSTEMS.items():
        sdf = generation_df[generation_df["system"] == system].dropna(subset=["generation_index"]).copy()
        if sdf.empty:
            continue
        ax.plot(
            sdf["generation_index"],
            sdf["best_worker_score_angstrom"],
            marker="o",
            color=meta["color"],
            label=f"{meta['label']} best worker",
        )
        ax.plot(
            sdf["generation_index"],
            sdf["mean_worker_max_score_angstrom"],
            marker=".",
            linestyle="--",
            color=meta["color"],
            alpha=0.55,
            label=f"{meta['label']} mean worker max",
        )
    ax.set_xlabel("Adaptive generation")
    ax.set_ylabel("PBC-corrected opening score / A")
    ax.set_title("Collagen unwinding progress by generation")
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=7, ncol=2)
    fig.tight_layout()
    fig.savefig(output_dir / "joint_generation_unwinding_progress.png", dpi=300)
    pdf.savefig(fig)
    plt.close(fig)


def plot_worker_completion(pdf: PdfPages, output_dir: Path, generation_df: pd.DataFrame) -> None:
    if generation_df.empty:
        return
    fig = plt.figure(figsize=(8.27, 5.85))
    ax = fig.add_subplot(111)
    for system, meta in SYSTEMS.items():
        sdf = generation_df[generation_df["system"] == system].dropna(subset=["generation_index"]).copy()
        if sdf.empty:
            continue
        ax.plot(sdf["generation_index"], sdf["n_completed"], marker="o", color=meta["color"], label=meta["label"])
    ax.set_xlabel("Adaptive generation")
    ax.set_ylabel("Completed workers")
    ax.set_title("Completed workers available per generation")
    ax.grid(True, alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_dir / "joint_completed_workers_by_generation.png", dpi=300)
    pdf.savefig(fig)
    plt.close(fig)


def plot_frame_score_envelopes(pdf: PdfPages, output_dir: Path, frame_df: pd.DataFrame) -> None:
    if frame_df.empty or "opening_score_angstrom" not in frame_df.columns:
        return
    fig = plt.figure(figsize=(8.27, 5.85))
    ax = fig.add_subplot(111)
    for system, meta in SYSTEMS.items():
        sdf = frame_df[frame_df["system"] == system].dropna(subset=["generation_index", "opening_score_angstrom"]).copy()
        if sdf.empty:
            continue
        grouped = sdf.groupby("generation_index")["opening_score_angstrom"]
        q = grouped.quantile([0.1, 0.5, 0.9]).unstack()
        ax.plot(q.index, q[0.5], marker="o", color=meta["color"], label=f"{meta['label']} median frame")
        ax.fill_between(q.index.to_numpy(dtype=float), q[0.1].to_numpy(dtype=float), q[0.9].to_numpy(dtype=float), color=meta["color"], alpha=0.12)
    ax.set_xlabel("Adaptive generation")
    ax.set_ylabel("Frame opening score / A")
    ax.set_title("Frame-score envelope per generation")
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(output_dir / "joint_frame_score_envelopes.png", dpi=300)
    pdf.savefig(fig)
    plt.close(fig)


def build_title_lines(specs: List[RunSpec], frame_df: pd.DataFrame, worker_df: pd.DataFrame, generation_df: pd.DataFrame) -> List[str]:
    lines = [
        "Joint adaptive MMP1-collagen unwinding monitor",
        "",
        "Metric: PBC-corrected C-alpha opening score reported by the MD runner.",
        "Scope: monitoring output only; no structural/contact/H-bond/PCA reanalysis is performed.",
        "",
        "Runs:",
    ]
    for spec in specs:
        lines.append(f"  {spec.label:10s} {spec.run_dir}")
    lines.extend(
        [
            "",
            f"Systems detected: {len(specs)}",
            f"Generations detected: {0 if generation_df.empty else len(generation_df)}",
            f"Workers detected: {0 if worker_df.empty else len(worker_df)}",
            f"Score frames detected: {0 if frame_df.empty else len(frame_df)}",
            "",
            "CSV outputs:",
            "  joint_generation_unwinding_summary.csv",
            "  joint_worker_unwinding_summary.csv",
            "  joint_opening_score_frames.csv",
            "  latest_generation_by_system.csv",
            "  best_worker_unwinding_by_generation.csv",
        ]
    )
    return lines


def add_latest_summary_page(pdf: PdfPages, generation_df: pd.DataFrame) -> None:
    latest = latest_per_system(generation_df)
    if latest.empty:
        return
    lines = [
        "Latest generation summary",
        "",
        "system       gen  completed  scored  best_A   mean_worker_max_A  failed",
    ]
    for _, row in latest.sort_values("system").iterrows():
        best = row.get("best_worker_score_angstrom", np.nan)
        mean = row.get("mean_worker_max_score_angstrom", np.nan)
        lines.append(
            f"{row['system_label']:<10s} "
            f"{int(row['generation_index']):>3d}  "
            f"{int(row['n_completed']):>9d}  "
            f"{int(row['n_workers_with_scores']):>6d}  "
            f"{best:>7.3f}  "
            f"{mean:>17.3f}  "
            f"{int(row['n_failed_archived']):>6d}"
        )
    add_text_page(pdf, "Joint Unwinding Report: Latest Status", lines)


def generate_pdf(output_pdf: Path, output_dir: Path, specs: List[RunSpec], frame_df: pd.DataFrame, worker_df: pd.DataFrame, generation_df: pd.DataFrame) -> None:
    with PdfPages(output_pdf) as pdf:
        add_text_page(pdf, "Joint Adaptive Collagen Unwinding Progress", build_title_lines(specs, frame_df, worker_df, generation_df))
        add_latest_summary_page(pdf, generation_df)
        plot_generation_progress(pdf, output_dir, generation_df)
        plot_frame_score_envelopes(pdf, output_dir, frame_df)
        plot_worker_completion(pdf, output_dir, generation_df)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate a joint PDF/CSV monitor report for WT and OI-mutant adaptive collagen unwinding runs.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--root", type=str, default=".", help="Repository root.")
    parser.add_argument("--wild_type_run", type=str, default=None, help="Explicit wild-type adaptive_run_* directory.")
    parser.add_argument("--g978s_run", type=str, default=None, help="Explicit G978S adaptive_run_* directory.")
    parser.add_argument("--g984c_run", type=str, default=None, help="Explicit G984C adaptive_run_* directory.")
    parser.add_argument("--g987r_run", type=str, default=None, help="Explicit G987R adaptive_run_* directory.")
    parser.add_argument("--output", type=str, default=None, help="Output PDF path.")
    parser.add_argument("--analysis_dir", type=str, default="joint_unwinding_report", help="Output directory for CSV and PNG files.")
    args = parser.parse_args()

    specs = resolve_run_specs(args)
    if not specs:
        raise SystemExit("No adaptive_run_* directories found. Pass explicit --wild_type_run/--g978s_run/--g984c_run/--g987r_run paths.")

    output_dir = Path(args.analysis_dir).resolve()
    safe_mkdir(output_dir)
    output_pdf = Path(args.output).resolve() if args.output else output_dir / "joint_collagen_unwinding_progress_report.pdf"

    frame_df, worker_df, generation_df = collect_run_tables(specs)
    write_csv_outputs(output_dir, frame_df, worker_df, generation_df)
    generate_pdf(output_pdf, output_dir, specs, frame_df, worker_df, generation_df)

    print(f"Systems analysed: {', '.join(spec.system for spec in specs)}")
    print(f"PDF written: {output_pdf}")
    print(f"Analysis directory: {output_dir}")
    print(f"Generation rows: {len(generation_df)}")
    print(f"Worker rows: {len(worker_df)}")
    print(f"Score frames: {len(frame_df)}")


if __name__ == "__main__":
    main()
