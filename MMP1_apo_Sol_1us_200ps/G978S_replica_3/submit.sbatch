#!/bin/bash
#SBATCH -J m1-g978s-r3
#SBATCH -N 1
#SBATCH -c 4
#SBATCH -p general
#SBATCH -q public
#SBATCH -G a100:1
#SBATCH -C a100_80
#SBATCH -t 7-00:00:00
#SBATCH --requeue
#SBATCH --signal=USR1@600
#SBATCH --open-mode=append
#SBATCH -o /scratch/anash19/MMP1_apo_Sol_1us_200ps/G978S_replica_3/slurm-%j.out
#SBATCH -e /scratch/anash19/MMP1_apo_Sol_1us_200ps/G978S_replica_3/slurm-%j.err
#SBATCH --mail-type=FAIL

module load mamba/latest
source activate openmm
python -u /scratch/anash19/MMP1_apo_Sol_1us_200ps/G978S_replica_3/run.py
