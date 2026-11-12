#!/usr/bin/env python3
"""Run the complete selected-frame structural analysis without reading trajectories.

This launcher combines the collagen hydrogen-bond/backbone analysis with the
residue-resolved opening and MMP-1 contact-context analysis. Each component reads
only the best_opening_frame.gro selected for a completed adaptive generation.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


RUN_ARGUMENTS = {
    "wild_type_run": "--wild_type_run",
    "g978s_run": "--g978s_run",
    "g984c_run": "--g984c_run",
    "g987r_run": "--g987r_run",
}


def git_commit(root: Path) -> str | None:
    result = subprocess.run(
        ["git", "-C", str(root), "rev-parse", "HEAD"],
        check=False,
        capture_output=True,
        text=True,
    )
    return result.stdout.strip() if result.returncode == 0 else None


def append_run_arguments(command: list[str], args: argparse.Namespace) -> None:
    for attribute, option in RUN_ARGUMENTS.items():
        value = getattr(args, attribute)
        if value is not None:
            command.extend([option, str(value.expanduser().resolve())])


def output_inventory(output_dir: Path) -> list[dict]:
    rows = []
    for path in sorted(output_dir.rglob("*")):
        if path.is_file():
            rows.append(
                {
                    "path": str(path.relative_to(output_dir)),
                    "bytes": path.stat().st_size,
                }
            )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run hydrogen-bond, backbone, opening-geometry and MMP-1 contact analyses on selected GRO frames.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--root", type=Path, default=Path("."), help="Repository root.")
    parser.add_argument("--output_dir", type=Path, default=Path("complete_structural_analysis"))
    parser.add_argument("--wild_type_run", type=Path)
    parser.add_argument("--g978s_run", type=Path)
    parser.add_argument("--g984c_run", type=Path)
    parser.add_argument("--g987r_run", type=Path)
    parser.add_argument("--max_hbond_da_angstrom", type=float, default=3.5)
    parser.add_argument("--min_hbond_angle_deg", type=float, default=120.0)
    parser.add_argument("--max_donor_hydrogen_angstrom", type=float, default=1.25)
    parser.add_argument("--mmp1_contact_cutoff_angstrom", type=float, default=4.5)
    args = parser.parse_args()

    root = args.root.expanduser().resolve()
    scripts = root / "scripts"
    output_dir = args.output_dir.expanduser().resolve()
    hbond_dir = output_dir / "hydrogen_bonds_and_backbone"
    context_dir = output_dir / "opening_and_mmp1_context"
    output_dir.mkdir(parents=True, exist_ok=True)

    hbond_command = [
        sys.executable,
        str(scripts / "report_collagen_structural_metrics.py"),
        "--root",
        str(root),
        "--analysis_dir",
        str(hbond_dir),
        "--max_hbond_da_a",
        str(args.max_hbond_da_angstrom),
        "--min_hbond_angle_deg",
        str(args.min_hbond_angle_deg),
        "--max_dh_a",
        str(args.max_donor_hydrogen_angstrom),
    ]
    append_run_arguments(hbond_command, args)

    context_command = [
        sys.executable,
        str(scripts / "analyse_collagen_unwinding_structures.py"),
        "--root",
        str(root),
        "--analysis_dir",
        str(context_dir),
        "--contact_cutoff_angstrom",
        str(args.mmp1_contact_cutoff_angstrom),
    ]
    append_run_arguments(context_command, args)

    print("Running periodic-boundary-corrected collagen hydrogen-bond and backbone analysis")
    subprocess.run(hbond_command, check=True)
    print("Running residue-resolved opening and MMP-1 contact-context analysis")
    subprocess.run(context_command, check=True)

    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "repository_root": str(root),
        "git_commit": git_commit(root),
        "trajectory_files_read": False,
        "input_scope": "best_opening_frame.gro selected for each completed adaptive generation",
        "criteria": {
            "maximum_hbond_donor_acceptor_angstrom": args.max_hbond_da_angstrom,
            "minimum_hbond_dha_angle_deg": args.min_hbond_angle_deg,
            "maximum_donor_hydrogen_angstrom": args.max_donor_hydrogen_angstrom,
            "mmp1_heavy_atom_contact_cutoff_angstrom": args.mmp1_contact_cutoff_angstrom,
        },
        "commands": {
            "hydrogen_bonds_and_backbone": hbond_command,
            "opening_and_mmp1_context": context_command,
        },
    }
    manifest_path = output_dir / "analysis_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    manifest["outputs"] = output_inventory(output_dir)
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    print(f"Complete structural analysis written to: {output_dir}")
    print(f"Reproducibility manifest: {manifest_path}")


if __name__ == "__main__":
    main()
