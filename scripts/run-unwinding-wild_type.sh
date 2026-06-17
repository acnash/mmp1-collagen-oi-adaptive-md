#!/bin/bash
#SBATCH -p public
#SBATCH -q public
#SBATCH -t 7-00:00:00
#SBATCH -c 16
#SBATCH -N 1
#SBATCH --mem=66G
#SBATCH -G a100:4
#SBATCH -J unwind_wild_type
#SBATCH -e run_unwinding_wild_type.err
#SBATCH -o run_unwinding_wild_type.out
#SBATCH --export=NONE

set -euo pipefail

module load mamba
source activate openmm

cd /scratch/anash19/mmp1-collagen-oi-adaptive-md

python scripts/adaptive_mmp1_unwinding_dual_worker.py \
    --input_dir generated_mutants/salted_150mM_NaCl/wild_type \
    --system_variant wild_type \
    --input_gro NPT_eq_wild_type_150mM_NaCl.gro \
    --input_top system.top \
    --platform CUDA \
    --concurrent_workers 4 \
    --mode auto \
    --target_generation 50
