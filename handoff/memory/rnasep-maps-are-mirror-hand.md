---
name: rnasep-maps-are-mirror-hand
description: FACT (2026-10-14): the RNase P consensus (rnasep_states/24_refine3D_auto/recvol_state01.mrc) and therefore every reconstruction in rnasep_coupled is the MIRROR hand of the deposited Degenhardt et al. maps (EMD-74131 docks at cc 0.99 into the z-flipped consensus vs 0.86 unflipped, dock_volpair lp 8); flip their maps/models before any comparison, and flip our maps before depositing or fitting models
metadata:
  type: project
---

Found while testing the closed TL-TLR state: the deposited open conformer docked into our consensus at 0.86 at 15 A
but the core correlation collapsed below 30 A (0.28 at 30-20 A, 0.08 at 20-12 A); a magnification scan did not help;
docking into the z-flipped consensus gave cc 0.9895 with the same in-plane rotation and no shift. The ab initio start
chose the mirror hand; heterogeneity results are unaffected, but all comparisons with deposited maps or atomic models
made before this check (section 18 of `coupled_2d3d_math_formulation_2026_10_13.md`: sub-class closed scores, site
correlations, the first two-template E-step) used mirror-image templates and were redone with mirrored deposited maps
(`emd_*_m160.mrc` in `rnasep_coupled/degenhardt`).

**How to apply:** before docking a deposited map or fitting a model to any RNase P map from this project, test both
hands (dock_volpair against the map and its z-flip) and use the hand with the higher cc; expect the flip. Related:
[[rnasep-24-refine3d-auto-quarter-unrefined]].
