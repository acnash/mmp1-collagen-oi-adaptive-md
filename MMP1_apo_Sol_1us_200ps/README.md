# MMP-1 apo collagen: nine 1 microsecond trajectories

Copy this complete directory from the repository to
`/scratch/anash19/MMP1_apo_Sol_1us_200ps/`. Before submission, verify the
portable inputs from inside that directory with:

```bash
sha256sum -c SHA256SUMS.txt
```

Each job continues one portable 5 ns endpoint for another 995 ns, giving 1,000 ns
of unrestrained production per replica. There are no positional restraints and
no adaptive sampling. The ensemble is NPT at 310.15 K and 1 bar with a 2 fs
timestep, using the same AMBER99SB-ILDN/TIP3P topology as the local runs.

Submit the nine independent GPU jobs with:

```bash
sbatch /scratch/anash19/MMP1_apo_Sol_1us_200ps/G978S_replica_1/submit.sbatch
sbatch /scratch/anash19/MMP1_apo_Sol_1us_200ps/G978S_replica_2/submit.sbatch
sbatch /scratch/anash19/MMP1_apo_Sol_1us_200ps/G978S_replica_3/submit.sbatch
sbatch /scratch/anash19/MMP1_apo_Sol_1us_200ps/G984C_replica_1/submit.sbatch
sbatch /scratch/anash19/MMP1_apo_Sol_1us_200ps/G984C_replica_2/submit.sbatch
sbatch /scratch/anash19/MMP1_apo_Sol_1us_200ps/G984C_replica_3/submit.sbatch
sbatch /scratch/anash19/MMP1_apo_Sol_1us_200ps/G987R_replica_1/submit.sbatch
sbatch /scratch/anash19/MMP1_apo_Sol_1us_200ps/G987R_replica_2/submit.sbatch
sbatch /scratch/anash19/MMP1_apo_Sol_1us_200ps/G987R_replica_3/submit.sbatch
```

Each batch file requests one A100 80 GB GPU for the seven-day maximum. Its
executable body contains only the Mamba module load, activation of the existing
`openmm` environment and one absolute Python command.

Every 200 ps, the job writes:

- a compressed NumPy snapshot in `output/snapshots_200ps/` containing every
  atom's coordinates and velocities, the periodic box, energies, temperature,
  production step and physical time;
- one row in `output/thermodynamics_200ps.csv`; and
- an atomic rolling binary checkpoint in `output/restart.chk`.

Open a snapshot with `numpy.load`. Arrays are named `positions_nm`,
`velocities_nm_per_ps` and `box_vectors_nm`; scalar metadata are stored in the
same file. The estimated snapshot storage is approximately 6.5 GB per replica,
or about 58 GB across all nine jobs.

If a job approaches the seven-day wall time, Slurm sends a warning. The Python
process finishes the current 200 ps block, commits the checkpoint and asks
Slurm to requeue the same job. A manually resubmitted job also resumes from
`restart.chk`. The workstation binary checkpoints are deliberately excluded;
the first Sol allocation starts from the portable XML state, and subsequent
allocations use only the Sol-generated checkpoint on the same requested GPU
model and OpenMM environment.

Progress is in each job's `output/progress.json`. Completion is recorded in
`output/complete.json` and `output/final_1us_state.xml`.
