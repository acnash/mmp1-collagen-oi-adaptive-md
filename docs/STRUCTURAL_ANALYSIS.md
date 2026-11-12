# Complete structural analysis of the adaptive unwinding campaigns

## Purpose

The complete workflow analyses one selected structure from every completed adaptive generation. It calculates:

- periodic-boundary-corrected collagen hydrogen bonds, separated into interchain and intrachain contacts;
- collagen backbone phi, psi and omega dihedral angles;
- residue-resolved three-chain C-alpha separations and triangular opening area for residues 977-987;
- local radius of gyration and aligned structural change from generation 0;
- residue-resolved MMP-1 heavy-atom contacts and minimum MMP-1 distances.

The workflow does not read, copy or package trajectory files. Its structural inputs are the compact `best_opening_frame.gro` files already selected by the adaptive simulations.

## Installation

The repository environment contains the required packages:

```bash
conda env create -f environment.yml
conda activate openmm
```

An existing Python environment may be used if it provides Python 3.9 or later, NumPy, pandas and Matplotlib.

## One-command server run

From the repository root on the simulation server:

```bash
cd /scratch/anash19/mmp1-collagen-oi-adaptive-md
git pull origin main

python scripts/run_complete_structural_analysis.py \
  --root . \
  --output_dir complete_structural_analysis
```

The launcher automatically selects the newest adaptive run found for wild type, G978S, G984C and G987R.

## Explicit run directories

If the latest directory is not the desired campaign, specify all four runs:

```bash
python scripts/run_complete_structural_analysis.py \
  --root . \
  --wild_type_run generated_mutants/salted_150mM_NaCl/wild_type/adaptive_run_wild_type_YYYYMMDD_HHMMSS \
  --g978s_run generated_mutants/salted_150mM_NaCl/collagen_G978S/adaptive_run_G978S_YYYYMMDD_HHMMSS \
  --g984c_run generated_mutants/salted_150mM_NaCl/collagen_G984C/adaptive_run_G984C_YYYYMMDD_HHMMSS \
  --g987r_run generated_mutants/salted_150mM_NaCl/collagen_G987R/adaptive_run_G987R_YYYYMMDD_HHMMSS \
  --output_dir complete_structural_analysis
```

## Outputs to download

The output directory contains only compact tables, plots, a PDF report and a JSON reproducibility manifest. The principal CSV files are:

```text
complete_structural_analysis/hydrogen_bonds_and_backbone/structural_best_frames.csv
complete_structural_analysis/hydrogen_bonds_and_backbone/structural_generation_summary.csv
complete_structural_analysis/hydrogen_bonds_and_backbone/structural_hydrogen_bonds.csv
complete_structural_analysis/hydrogen_bonds_and_backbone/structural_backbone_dihedrals.csv

complete_structural_analysis/opening_and_mmp1_context/residue_resolved_opening_geometry.csv
complete_structural_analysis/opening_and_mmp1_context/residue_resolved_mmp1_contacts.csv
complete_structural_analysis/opening_and_mmp1_context/structural_frame_summary.csv
complete_structural_analysis/opening_and_mmp1_context/structural_changes_from_generation_zero.csv

complete_structural_analysis/analysis_manifest.json
```

For a compact transfer archive containing no trajectories:

```bash
tar -czf complete_structural_analysis_csvs.tar.gz \
  complete_structural_analysis/analysis_manifest.json \
  complete_structural_analysis/hydrogen_bonds_and_backbone/*.csv \
  complete_structural_analysis/opening_and_mmp1_context/*.csv
```

## Default structural criteria

- Hydrogen-bond donor-acceptor distance: at most 3.5 angstrom.
- Hydrogen-bond D-H-A angle: at least 120 degrees.
- Donor-hydrogen assignment distance: at most 1.25 angstrom.
- MMP-1 heavy-atom contact distance: at most 4.5 angstrom.

All distances and angles used for hydrogen bonds, backbone dihedrals and MMP-1 contacts are corrected for periodic boundary conditions. The JSON manifest records the thresholds, commands and Git commit used for each analysis.
