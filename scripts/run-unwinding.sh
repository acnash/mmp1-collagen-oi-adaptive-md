#!/bin/bash
#SBATCH -p public
#SBATCH -q public
#SBATCH -t 7-00:00:00
#SBATCH -c 16
#SBATCH -N 1
#SBATCH --mem=66G
#SBATCH -G a100:4
#SBATCH -J unwind_wild_type
#SBATCH -e run_unwinding.err
#SBATCH -o run_unwinding.out
#SBATCH --export=NONE

set -euo pipefail

# Template launcher. Edit SYSTEM_VARIANT, INPUT_DIR, and INPUT_GRO together.
SYSTEM_VARIANT="wild_type"
INPUT_DIR="generated_mutants/salted_150mM_NaCl/wild_type"
INPUT_GRO="NPT_eq_wild_type_150mM_NaCl.gro"

# Other system options:
# SYSTEM_VARIANT="G978S"
# INPUT_DIR="generated_mutants/salted_150mM_NaCl/collagen_G978S"
# INPUT_GRO="NPT_eq_collagen_G978S_150mM_NaCl_openmm_10ps_nvt_10ps_npt_with_velocities.gro"
#
# SYSTEM_VARIANT="G984C"
# INPUT_DIR="generated_mutants/salted_150mM_NaCl/collagen_G984C"
# INPUT_GRO="NPT_eq_collagen_G984C_150mM_NaCl_openmm_10ps_nvt_10ps_npt_with_velocities.gro"
#
# SYSTEM_VARIANT="G987R"
# INPUT_DIR="generated_mutants/salted_150mM_NaCl/collagen_G987R"
# INPUT_GRO="NPT_eq_collagen_G987R_150mM_NaCl_openmm_10ps_nvt_10ps_npt_with_velocities.gro"

module load mamba
source activate openmm

cd /scratch/anash19/mmp1-collagen-oi-adaptive-md

python scripts/adaptive_mmp1_unwinding_dual_worker.py \
    --input_dir "${INPUT_DIR}" \
    --system_variant "${SYSTEM_VARIANT}" \
    --input_gro "${INPUT_GRO}" \
    --input_top system.top \
    --platform CUDA \
    --concurrent_workers 4 \
    --mode auto \
    --target_generation 50
