# Prepared simulation systems

This directory separates the two production programmes by physical system and
sampling method.

| Directory | Enzyme state | Systems | Intended use |
|---|---|---|---|
| [`mmp1-bound-adaptive/`](mmp1-bound-adaptive/) | MMP-1 present | wild type, G978S, G984C and G987R | adaptive sampling of local collagen opening |
| [`mmp1-free-apo/`](mmp1-free-apo/) | MMP-1 absent | G978S, G984C and G987R | conventional, non-adaptive apo-collagen dynamics |

The separately named [`MMP1_apo_Sol_1us_200ps/`](../MMP1_apo_Sol_1us_200ps/)
directory is the immutable deployment bundle for the nine Sol jobs. It retains
its original paths so that this repository reorganisation cannot affect jobs
already copied to, or running on, Sol.
