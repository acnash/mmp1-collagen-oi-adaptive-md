# Adaptive MMP-1/collagen workflow reference

This repository contains prepared GROMACS/OpenMM systems and scripts for adaptive molecular dynamics of human MMP-1 bound to the synthetic homotrimeric model collagen peptide used in PDB 4AUO. The scientific aim is to compare the unmodified model peptide with three single-chain glycine substitutions near its MMP-1 collagenase cleavage region. The peptide is treated as a general model collagen substrate rather than a particular native collagen type. Residue numbers follow the 4AUO model-peptide coordinate register and are not COL1A1 or COL1A2 HGVS positions.

## Environment

Create the documented `openmm` Conda environment from the repository root:

```bash
conda env create -f environment.yml
conda activate openmm
```

The environment includes OpenMM and the NumPy, pandas, and Matplotlib dependencies used by the simulation, monitoring, and structural-report scripts. GPU execution additionally requires a compatible NVIDIA driver on the host system.

## Systems

The modeled collagen register follows the numbering assigned to the 4AUO model collagen peptide in the prepared topology. The peptide contains a collagenase-recognition region and terminal Gly-Pro-Hyp repeats that stabilize its triple-helical structure. The scissile Gly-Leu bond is Gly975-Leu976 in this topology. The analysis window, residues 977-987, corresponds to the P2′-P12′ region; G978S, G984C and G987R are therefore model substitutions at P3′, P9′ and P12′, respectively. They are model-peptide coordinates, not established COL1A1 or COL1A2 clinical-variant coordinates. The prepared systems are:

| System | Mutated chain | Collagen residue change | Input folder |
|---|---:|---|---|
| `wild_type` | none | none | `systems/mmp1-bound-adaptive/wild_type/` |
| `G978S` | chain 1 | `GLY978 -> SER978` | `systems/mmp1-bound-adaptive/collagen_G978S/` |
| `G984C` | chain 1 | `GLY984 -> CYS984` | `systems/mmp1-bound-adaptive/collagen_G984C/` |
| `G987R` | chain 1 | `GLY987 -> ARG987` | `systems/mmp1-bound-adaptive/collagen_G987R/` |

Each system directory contains a complete local GROMACS include set:

```text
system.top
topol_MMP1_active.itp
collagen.itp
collagen_G*.itp       # mutant systems only
forcefield.itp
ffbonded.itp
ffnonbonded.itp
tip3p.itp
ions.itp
*.gro
```

Position-restraint include calls have been removed from the committed `system.top` files. Restraints should be added programmatically in OpenMM workflows when needed.

## Salt And Charge

The systems were salted to a target of 150 mM NaCl by replacing water molecules with ions. For this box volume, that corresponds to 57 NaCl pairs. Existing neutralizing chloride ions were preserved.

Final molecule counts:

| System | SOL | NA | CL | Net charge |
|---|---:|---:|---:|---:|
| `wild_type` | 18745 | 57 | 67 | ~0 |
| `G978S` | 18745 | 57 | 67 | ~0 |
| `G984C` | 18745 | 57 | 67 | ~0 |
| `G987R` | 18744 | 57 | 68 | ~0 |

## OpenMM Validation Runs

Short validation runs were attempted locally in the `kcc2` conda environment using OpenMM CPU platform with 8 CPU threads. The protocol was:

```text
energy minimization
10 ps NVT
10 ps NPT
```

The original 2 fs startup was unstable for `G978S`, so validation was repeated with a 0.5 fs timestep. That protocol completed successfully for `G978S`.

Validated output:

```text
systems/mmp1-bound-adaptive/collagen_G978S/
  NPT_eq_collagen_G978S_150mM_NaCl_openmm_10ps_npt.gro
  openmm_validation/
    G978S_after_minimization.gro
    G978S_after_10ps_nvt.gro
    NPT_eq_collagen_G978S_150mM_NaCl_openmm_10ps_npt.gro
    nvt.log
    npt.log
    summary.json
```

Initial `G984C` and `G987R` structures did not pass the same validation. Both developed NaN coordinates during NVT. Additional diagnostics showed:

```text
G984C: failed during NVT at 0.5 fs; also failed at 10 K / 0.1 fs.
G987R: failed during NVT at 0.5 fs; also failed at 10 K / 0.1 fs.
```

Those systems were later rescued by a more aggressive OpenMM-only relaxation campaign followed by 10 ps NVT and 10 ps NPT at 312.5 K with a 2 fs timestep. Final velocity-bearing structures are committed for the three mutant systems:

```text
systems/mmp1-bound-adaptive/collagen_G978S/
  NPT_eq_collagen_G978S_150mM_NaCl_openmm_10ps_nvt_10ps_npt_with_velocities.gro

systems/mmp1-bound-adaptive/collagen_G984C/
  NPT_eq_collagen_G984C_150mM_NaCl_openmm_10ps_nvt_10ps_npt_with_velocities.gro

systems/mmp1-bound-adaptive/collagen_G987R/
  NPT_eq_collagen_G987R_150mM_NaCl_openmm_10ps_nvt_10ps_npt_with_velocities.gro
```

Validation scripts:

```text
scripts/openmm_validate_mutants.py
scripts/openmm_diagnose_mutant.py
```

Example validation command:

```bash
conda run -n kcc2 python scripts/openmm_validate_mutants.py \
  --threads 8 \
  --timestep_ps 0.0005 \
  --variants G978S G984C G987R
```

## Adaptive Unwinding Runner

The adaptive runner is:

```text
scripts/adaptive_mmp1_unwinding_dual_worker.py
```

It launches generations of independent OpenMM workers. Generation 0 performs minimization, NVT, NPT, and production MD. Later generations use single-parent conformational seeding: the highest-opening parent frame seeds all workers in the next generation, using kinetic-energy-matched stochastic velocity reinitialisation. During production, a PBC-corrected collagen-opening metric is reported.

The opening metric uses residues 977-987 from all three collagen chains. For each residue position, it measures all three pairwise C-alpha distances between collagen chains using minimum-image periodic distances. The score is the mean of 33 distances:

```text
11 residue positions x 3 inter-chain distances = 33 distances
```

The runner is now system-aware:

```bash
--system_variant wild_type|G978S|G984C|G987R
--input_gro <coordinate.gro>
--input_top system.top
```

Fresh run directories include the variant label:

```text
adaptive_run_wild_type_YYYYMMDD_HHMMSS
adaptive_run_G978S_YYYYMMDD_HHMMSS
adaptive_run_G984C_YYYYMMDD_HHMMSS
adaptive_run_G987R_YYYYMMDD_HHMMSS
```

The runner validates that the opening-score window contains 33 CA atoms across the three collagen chains and that the expected mutation is present on chain 1 only.

## SLURM Launchers

Use the system-specific launchers:

```bash
cd /scratch/anash19/mmp1-collagen-oi-adaptive-md
sbatch scripts/run-unwinding-wild_type.sh
sbatch scripts/run-unwinding-G978S.sh
sbatch scripts/run-unwinding-G984C.sh
sbatch scripts/run-unwinding-G987R.sh
```

Each launcher is deliberately explicit: it loads the `mamba` module, activates the `openmm` environment, changes into `/scratch/anash19/mmp1-collagen-oi-adaptive-md`, then runs one `python scripts/adaptive_mmp1_unwinding_dual_worker.py ...` command.

The generic template is also available:

```bash
cd /scratch/anash19/mmp1-collagen-oi-adaptive-md
sbatch scripts/run-unwinding.sh
```

For the generic template, edit `SYSTEM_VARIANT`, `INPUT_DIR`, and `INPUT_GRO` together near the top of the file:

```bash
SYSTEM_VARIANT="wild_type"
INPUT_DIR="systems/mmp1-bound-adaptive/wild_type"
INPUT_GRO="NPT_eq_wild_type_150mM_NaCl.gro"
```

Each bespoke launcher has unique SLURM job, stdout, and stderr names. The mutant launchers use the final 10 ps NVT + 10 ps NPT velocity-bearing `.gro` files as their fresh-run inputs.

## Joint Unwinding Monitor

The joint unwinding monitor is:

```text
scripts/monitor_collagen_unwinding_progress.py
```

It reads only the adaptive runner's unwinding outputs: per-worker `opening_scores.csv`, `worker.json`, and per-generation `generation.json`. It does not perform structural clustering, contact analysis, hydrogen-bond analysis, PCA, or trajectory reanalysis.

From the repository root, it can auto-detect the newest adaptive run for wild type and the three OI mutants:

```bash
python scripts/monitor_collagen_unwinding_progress.py \
  --root . \
  --analysis_dir joint_unwinding_report
```

For explicit resumed campaigns, pass the four run directories:

```bash
python scripts/monitor_collagen_unwinding_progress.py \
  --wild_type_run systems/mmp1-bound-adaptive/wild_type/adaptive_run_wild_type_YYYYMMDD_HHMMSS \
  --g978s_run systems/mmp1-bound-adaptive/collagen_G978S/adaptive_run_G978S_YYYYMMDD_HHMMSS \
  --g984c_run systems/mmp1-bound-adaptive/collagen_G984C/adaptive_run_G984C_YYYYMMDD_HHMMSS \
  --g987r_run systems/mmp1-bound-adaptive/collagen_G987R/adaptive_run_G987R_YYYYMMDD_HHMMSS \
  --analysis_dir joint_unwinding_report
```

The default PDF is `joint_unwinding_report/joint_collagen_unwinding_progress_report.pdf`. The monitor also writes joint frame, worker, generation, and latest-generation CSV summaries. For independent plotting, it writes `best_worker_unwinding_by_generation.csv` with five columns: generation, wild type, G978S, G984C, and G987R best-worker opening scores.

## Joint Structural Report

The structural report script is:

```text
scripts/report_collagen_structural_metrics.py
```

It uses the best worker/frame selected for each adaptive generation and analyzes collagen only. It reports hydrogen-bond distributions, backbone `phi`/`psi`/`omega` dihedral distributions, and plots structural metrics against the PBC-corrected unwinding score.

Run it from the repository root with auto-detected newest runs:

```bash
python scripts/report_collagen_structural_metrics.py \
  --root . \
  --analysis_dir joint_structural_report
```

For explicit resumed campaigns:

```bash
python scripts/report_collagen_structural_metrics.py \
  --wild_type_run systems/mmp1-bound-adaptive/wild_type/adaptive_run_wild_type_YYYYMMDD_HHMMSS \
  --g978s_run systems/mmp1-bound-adaptive/collagen_G978S/adaptive_run_G978S_YYYYMMDD_HHMMSS \
  --g984c_run systems/mmp1-bound-adaptive/collagen_G984C/adaptive_run_G984C_YYYYMMDD_HHMMSS \
  --g987r_run systems/mmp1-bound-adaptive/collagen_G987R/adaptive_run_G987R_YYYYMMDD_HHMMSS \
  --analysis_dir joint_structural_report
```

The default PDF is `joint_structural_report/joint_collagen_structural_report.pdf`. CSV outputs include `structural_best_frames.csv`, `structural_generation_summary.csv`, `structural_hydrogen_bonds.csv`, and `structural_backbone_dihedrals.csv`.

## Compact Structural Analysis

The compact structural analysis reads only the selected `best_opening_frame.gro` from each completed generation. It does not read or transfer trajectory files. It reports residue-resolved three-chain opening geometry, local collagen radius of gyration, aligned change from generation 0, and MMP-1 heavy-atom contacts:

```bash
python scripts/analyse_collagen_unwinding_structures.py \
  --root . \
  --analysis_dir compact_structural_analysis
```

The four CSV outputs are `residue_resolved_opening_geometry.csv`, `residue_resolved_mmp1_contacts.csv`, `structural_frame_summary.csv`, and `structural_changes_from_generation_zero.csv`.

## Complete Structural Analysis

The complete workflow calculates periodic-boundary-corrected collagen hydrogen bonds, interchain and intrachain bond counts, backbone dihedrals, residue-resolved opening geometry and MMP-1 contact context. It reads selected `best_opening_frame.gro` structures only and does not read or transfer trajectories.

On the simulation server, update the repository and run the analysis with:

```bash
cd /scratch/anash19/mmp1-collagen-oi-adaptive-md
git pull origin main

python scripts/run_complete_structural_analysis.py \
  --root . \
  --output_dir complete_structural_analysis
```

The principal downloadable tables are:

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

## Peer-review trajectory analysis

`run_peer_review_trajectory_analysis_standalone.py` extends the structural analysis across existing production trajectories. It does not run new molecular dynamics. By default it samples every 200 ps from all workers at generations 0, 10, 20 and 28 and reports chain-resolved opening, hydrogen-bond sensitivity and occupancy, residue-specific dihedral statistics, unique collagen-MMP-1 residue contacts, MMP-domain contacts, catalytic-zinc coordination and scissile-region distances.

When the server is not connected to GitHub, upload this one Python file into the repository root with the SOL web portal and run this single command from the repository root:

```bash
python run_peer_review_trajectory_analysis_standalone.py --root . --output_dir peer_review_trajectory_analysis
```

The script creates `peer_review_trajectory_analysis.zip` for download through the web portal. Its default 200 ps sampling interval is intended to make the full four-system analysis practical while retaining worker-level replication. Use `--frame_stride 10` for 100 ps sampling if additional temporal resolution is required.

Create a compact archive containing the CSV data and reproducibility manifest, but no trajectories:

```bash
tar -czf complete_structural_analysis_csvs.tar.gz \
  complete_structural_analysis/analysis_manifest.json \
  complete_structural_analysis/hydrogen_bonds_and_backbone/*.csv \
  complete_structural_analysis/opening_and_mmp1_context/*.csv
```

To create a fresh checkout on a different machine:

```bash
git clone https://github.com/acnash/mmp1-collagen-oi-adaptive-md.git
cd mmp1-collagen-oi-adaptive-md
```

The manuscript is then available under `manuscript/`. Full structural-analysis instructions, including explicit adaptive-run directory options and analysis thresholds, are also retained in [`docs/STRUCTURAL_ANALYSIS.md`](docs/STRUCTURAL_ANALYSIS.md).

## PyMOL open-state scenes

`analysis/pymol_open_state_scripts/all_open_states_standalone.py` is a self-contained PyMOL script containing the wild-type, G978S, G984C and G987R generation-28 selected structures. Run it through PyMOL's **File > Run Script** command. It creates four named scenes with light-grey MMP-1, purple collagen and conventionally coloured residue sticks. The structures are continuous and no periodic or crystallographic copies are generated.

## Manuscript draft

The current manuscript draft is available in Markdown and PDF under [`manuscript/`](manuscript/). The draft reports the generation-matched opening, structural, hydrogen-bond, backbone-dihedral and MMP-1-contact analyses through generation 28.
