# MMP-1-bound adaptive collagen sampling

## Question

How does replacing glycine within the model collagen triple helix alter the
local opening pathways accessible when collagen is bound to MMP-1?

## Systems

The programme compares the unmodified 4AUO model collagen peptide with three
single-chain substitutions: G978S, G984C and G987R. These labels refer to the
model-peptide coordinate register, not to `COL1A1` or `COL1A2` clinical variant
coordinates.

## Approach

The systems contain MMP-1, periodically connected model collagen, TIP3P water
and approximately physiological NaCl. Adaptive OpenMM simulations select parent
frames using a periodic-boundary-corrected collagen-opening score over model
residues 977-987. Structural analyses examine interchain separation,
hydrogen-bonding, backbone geometry and MMP-1 contact context.

## Interpretation boundary

Adaptive selection is intended to explore and compare rare opening pathways.
It does not yield an unbiased equilibrium ensemble and should not be used alone
to infer spontaneous opening rates, cleavage kinetics or clinical causality.

## Resources

- [Detailed adaptive workflow](../../docs/ADAPTIVE_WORKFLOW_REFERENCE.md)
- [Structural-analysis methods](../../docs/STRUCTURAL_ANALYSIS.md)
- [Working manuscript](../../manuscript/MMP1_collagen_OI_adaptive_unwinding_manuscript.md)
- [Simulation runner](../../scripts/adaptive_mmp1_unwinding_dual_worker.py)
