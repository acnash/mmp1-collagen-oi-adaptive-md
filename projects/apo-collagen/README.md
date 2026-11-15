# MMP-1-free apo-collagen simulations

## Question

Can substituted model collagen populate locally opened conformations without
MMP-1 and without adaptive selection?

## Design

MMP-1 was removed from compact, pre-adaptive collagen coordinates. Each system
was freshly solvated with TIP3P water, neutralised and brought to approximately
150 mM NaCl while retaining the existing periodic collagen topology and force
field. Three independent replicas were prepared for each substituted system:

- G978S
- G984C
- G987R

Replicas share their starting coordinates within a substitution but use
different initial velocities and stochastic-integrator seeds. After energy
minimisation and staged restraint release, production dynamics proceed without
positional restraints or adaptive reseeding.

## Current production programme

The nine substituted-collagen trajectories have completed local validation and
short unrestrained continuation. Long-timescale GPU simulations are being used
to extend each trajectory towards one microsecond, with complete coordinates,
velocities and thermodynamic observables recorded at regular intervals.

## Interpretation boundary

These simulations can identify spontaneous structural opening within the model
systems. They cannot establish a mutation-specific effect until a matched
wild-type apo ensemble is simulated and analysed under the same protocol.
Neither spontaneous opening nor altered remodelling, by itself, establishes a
causal mechanism for OI pathology.

## Reproducibility resources

- [Prepare the MMP-1-free systems](../../scripts/prepare_mmp1_free_collagen.py)
- [Run staged equilibration](../../scripts/run_mmp1_free_collagen_equilibration.py)
- [Run local unrestrained continuation](../../scripts/run_mmp1_free_5ns_production.py)
- [Audit the equilibrated replicas](../../scripts/audit_mmp1_free_equilibration.py)
- [Build the explicit Sol production bundle](../../scripts/build_mmp1_apo_sol_1us_package.py)

Large trajectories and server-output snapshots are excluded from the public Git
repository. Portable preparation scripts, topology records, checksums and
analysis code should accompany a formal data release.
