#!/bin/bash
#SBATCH -p public
#SBATCH -q public
#SBATCH -t 7-00:00:00
#SBATCH -c 16
#SBATCH -N 1
#SBATCH --mem=66G
#SBATCH -G a100:4
#SBATCH -J unwind_G984C
#SBATCH -e run_unwinding_G984C.err
#SBATCH -o run_unwinding_G984C.out
#SBATCH --export=NONE

set -euo pipefail

# Run from a clone of the repository on SLURM scratch space:
#   cd /scratch/$USER/mmp1-collagen-oi-adaptive-md
#   sbatch scripts/run-unwinding-G984C.sh
#
# Resume example:
#   RUN_DIR="${REPO_ROOT}/generated_mutants/salted_150mM_NaCl/collagen_G984C/adaptive_run_G984C_YYYYMMDD_HHMMSS"
#   Leave RUN_DIR="" for a fresh run.

SYSTEM_VARIANT="G984C"
RUN_DIR=""

TARGET_GENERATION=50
CONCURRENT_WORKERS=4
PLATFORM="CUDA"
MODE="auto"

module load mamba
source activate openmm

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
PYTHON_SCRIPT="${REPO_ROOT}/scripts/adaptive_mmp1_unwinding_dual_worker.py"

INPUT_DIR="${REPO_ROOT}/generated_mutants/salted_150mM_NaCl/collagen_G984C"
INPUT_GRO="NPT_eq_collagen_G984C_150mM_NaCl_openmm_10ps_nvt_10ps_npt_with_velocities.gro"

echo "Repository: ${REPO_ROOT}"
echo "System variant: ${SYSTEM_VARIANT}"
echo "Input directory: ${INPUT_DIR}"
echo "Input GRO: ${INPUT_GRO}"
echo "Target generation: ${TARGET_GENERATION}"
echo "Concurrent workers: ${CONCURRENT_WORKERS}"
echo "Platform: ${PLATFORM}"

if [[ -n "${RUN_DIR}" ]]; then
    conda run -n openmm python "${PYTHON_SCRIPT}" \
        --resume "${RUN_DIR}" \
        --platform "${PLATFORM}" \
        --concurrent_workers "${CONCURRENT_WORKERS}" \
        --mode "${MODE}" \
        --target_generation "${TARGET_GENERATION}"
else
    conda run -n openmm python "${PYTHON_SCRIPT}" \
        --input_dir "${INPUT_DIR}" \
        --system_variant "${SYSTEM_VARIANT}" \
        --input_gro "${INPUT_GRO}" \
        --input_top system.top \
        --platform "${PLATFORM}" \
        --concurrent_workers "${CONCURRENT_WORKERS}" \
        --mode "${MODE}" \
        --target_generation "${TARGET_GENERATION}"
fi
