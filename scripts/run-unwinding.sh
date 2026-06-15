#!/bin/bash
#SBATCH -p public
#SBATCH -q public
#SBATCH -t 7-00:00:00
#SBATCH -c 16
#SBATCH -N 1
#SBATCH --mem=66G
#SBATCH -G a100:4
#SBATCH -e run_unwinding.err
#SBATCH -o run_unwinding.out
#SBATCH --export=NONE

set -euo pipefail

# Run this script from a clone of the repository on SLURM scratch space:
#   cd /scratch/$USER/mmp1-collagen-oi-adaptive-md
#   sbatch scripts/run-unwinding.sh
#
# Available SYSTEM_VARIANT values:
#   wild_type
#   G978S
#   G984C
#   G987R
#
# Fresh-run examples. Change SYSTEM_VARIANT below to one of:
#   SYSTEM_VARIANT="wild_type"
#   SYSTEM_VARIANT="G978S"
#   SYSTEM_VARIANT="G984C"
#   SYSTEM_VARIANT="G987R"
#
# Resume example:
#   Set RUN_DIR to an existing adaptive run directory, for example:
#   RUN_DIR="${REPO_ROOT}/generated_mutants/salted_150mM_NaCl/collagen_G978S/adaptive_run_G978S_20260615_120000"
#   Leave RUN_DIR="" for a fresh run.

SYSTEM_VARIANT="wild_type"
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

case "${SYSTEM_VARIANT}" in
    wild_type)
        INPUT_DIR="${REPO_ROOT}/generated_mutants/salted_150mM_NaCl/wild_type"
        INPUT_GRO="NPT_eq_wild_type_150mM_NaCl.gro"
        ;;
    G978S)
        INPUT_DIR="${REPO_ROOT}/generated_mutants/salted_150mM_NaCl/collagen_G978S"
        INPUT_GRO="NPT_eq_collagen_G978S_150mM_NaCl.gro"
        ;;
    G984C)
        INPUT_DIR="${REPO_ROOT}/generated_mutants/salted_150mM_NaCl/collagen_G984C"
        INPUT_GRO="NPT_eq_collagen_G984C_150mM_NaCl.gro"
        ;;
    G987R)
        INPUT_DIR="${REPO_ROOT}/generated_mutants/salted_150mM_NaCl/collagen_G987R"
        INPUT_GRO="NPT_eq_collagen_G987R_150mM_NaCl.gro"
        ;;
    *)
        echo "Unknown SYSTEM_VARIANT: ${SYSTEM_VARIANT}" >&2
        echo "Choose one of: wild_type, G978S, G984C, G987R" >&2
        exit 2
        ;;
esac

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
