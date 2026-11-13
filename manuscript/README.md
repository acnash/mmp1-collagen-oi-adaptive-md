# Manuscript draft

This directory contains the portable manuscript draft generated from the joint adaptive-unwinding report:

- `MMP1_collagen_OI_adaptive_unwinding_manuscript.md`: editable Markdown source;
- `MMP1_collagen_OI_adaptive_unwinding_manuscript.pdf`: rendered journal-style draft;
- `figures/`: manuscript figures referenced by the Markdown file.

The draft compares generations 0-28 across the unmodified synthetic model collagen peptide, G978S, G984C and G987R. Adaptive maxima are interpreted as evidence of pathway accessibility, not as equilibrium probabilities or free energies.

The draft includes generation-matched selected-structure analysis through generation 28 and worker-wide analysis at generations 0, 10, 20 and 28. It reports opening, hydrogen bonding, backbone conformational dispersion, collagen-MMP-1 contacts and catalytic-site geometry. It identifies the 4AUO substrate as a synthetic homotrimeric model collagen peptide without assigning it to a particular native collagen type, distinguishes aggregate worker time from continuous sampling and records the current model and preparation limitations. Regenerate the PDF from the Markdown source with `python analysis/build_manuscript_pdf.py` from the repository root.
