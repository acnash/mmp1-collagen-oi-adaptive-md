# MMP-1 / Osteogenesis Imperfecta Adaptive Collagen Unwinding

This repository contains prepared GROMACS/OpenMM systems and scripts for adaptive molecular dynamics of human MMP-1 bound to a type-I-collagen-like triple-helical collagen peptide. The scientific aim is to compare wild-type collagen with clinically motivated osteogenesis imperfecta glycine substitutions near the MMP-1 collagenase cleavage region.

## Systems

The modeled collagen register follows the 4AUO collagen peptide numbering. The MMP-1 scissile-region glycine is residue 981 in the 4AUO coordinate register. The prepared systems are:

| System | Mutated chain | Collagen residue change | Input folder |
|---|---:|---|---|
| `wild_type` | none | none | `generated_mutants/salted_150mM_NaCl/wild_type/` |
| `G978S` | chain 1 | `GLY978 -> SER978` | `generated_mutants/salted_150mM_NaCl/collagen_G978S/` |
| `G984C` | chain 1 | `GLY984 -> CYS984` | `generated_mutants/salted_150mM_NaCl/collagen_G984C/` |
| `G987R` | chain 1 | `GLY987 -> ARG987` | `generated_mutants/salted_150mM_NaCl/collagen_G987R/` |

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
generated_mutants/salted_150mM_NaCl/collagen_G978S/
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
generated_mutants/salted_150mM_NaCl/collagen_G978S/
  NPT_eq_collagen_G978S_150mM_NaCl_openmm_10ps_nvt_10ps_npt_with_velocities.gro

generated_mutants/salted_150mM_NaCl/collagen_G984C/
  NPT_eq_collagen_G984C_150mM_NaCl_openmm_10ps_nvt_10ps_npt_with_velocities.gro

generated_mutants/salted_150mM_NaCl/collagen_G987R/
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
INPUT_DIR="generated_mutants/salted_150mM_NaCl/wild_type"
INPUT_GRO="NPT_eq_wild_type_150mM_NaCl.gro"
```

Each bespoke launcher has unique SLURM job, stdout, and stderr names. The mutant launchers use the final 10 ps NVT + 10 ps NPT velocity-bearing `.gro` files as their fresh-run inputs.
