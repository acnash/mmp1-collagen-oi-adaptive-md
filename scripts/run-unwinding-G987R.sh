#!/bin/bash
#SBATCH -p public
#SBATCH -q public
#SBATCH -t 7-00:00:00
#SBATCH -c 16
#SBATCH -N 1
#SBATCH --mem=66G
#SBATCH -G a100:4
#SBATCH -J unwind_G987R
#SBATCH -e run_unwinding_G987R.err
#SBATCH -o run_unwinding_G987R.out
#SBATCH --export=NONE

set -euo pipefail

module load mamba
source activate openmm

cd /scratch/anash19/mmp1-collagen-oi-adaptive-md

python scripts/adaptive_mmp1_unwinding_dual_worker.py \
    --input_dir generated_mutants/salted_150mM_NaCl/collagen_G987R \
    --system_variant G987R \
    --input_gro NPT_eq_collagen_G987R_150mM_NaCl_openmm_10ps_nvt_10ps_npt_with_velocities.gro \
    --input_top system.top \
    --platform CUDA \
    --concurrent_workers 4 \
    --mode auto \
    --target_generation 50
