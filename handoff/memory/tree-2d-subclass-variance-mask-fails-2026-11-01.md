---
name: tree-2d-subclass-variance-mask-fails-2026-11-01
description: "VERDICT 2026-11-01: cls_expansion 2D sub-classes as the children of the node's cross-half variance map (mask + discover_sites) FAIL on r131 and the light arms: the bin-centred covariance is noise peaking at the molecule centre (peak/null 3-7 vs cheap children 19, which peak INSIDE the closed arm / arm 4); recall of the arms 0.00 at every matched mask size; 7-20 core-h vs 1; not wired in"
metadata:
  type: project
---

2026-11-01, side session mystifying-mcclintock (branch claude/tree-2d-subclass-mask off claude/tree-dist 9005a8c77), note
doc/implementation_notes/tree_2d_subclass_mask_test_2026_11.md; data hetero_tree_20261028/subcls_mask/ on batch2; scripts
scripts/cls_expansion_proto/tree_subclass_mask/ (Python port of half_variance / falloff_mask / discover_sites: Dice 0.98-1.00 and
variance corr 1.0000 against the tree's own mask.mrc / sites.txt, so the port is the Fortran).
- Fixtures: r131 (guardfix fixture run, 44k rows, maps kept) and a fresh arms root (routes/toy55/simarms, 160k non-junk rows,
  maxdepth 0, dist build12, 26 core-h). Both baselines are the audit's failure: the cheap-children mask is 93 % (r131) / 100 %
  (arms) of the molecule under the 0.2-of-peak core.
- Sub-classes (pose-mean cls_expansion, 8 per direction bin, bins = n/1300 as cheap_children): 272 (r131) / 984 (arms) single-view
  children, sizes median 36 / 108, min 1. Their covariance must be bin-centred (vs the node mean the missing wedge is the deviation).
- Result: the bin-centred sub-class variance peaks at the molecule CENTRE (13 A / 5 A from the centre), peak/null 3.1 / 7.0, and
  ranks the closed arm / the four arms at chance or below at every matched size (top 2-40 % of the support: recall 0.00-0.26 r131,
  0.00 arms); discovery fragments (48 sites at r131, best closed-arm coverage 0.22; 0-2 off-arm sites on arms). A 100-particle size
  floor (86 / 525 children, 89 % of the particles) changes nothing. The cheap-children variance peaks inside the closed arm (65 A,
  peak/null 19) / inside arm 4 (50 A) with 4x precision enrichment at the top 2 %: its only defect is the 20 % core rule.
- Why: v = sum_j w_j dE_j dO_j with w_j ~ n_j gives every child the same noise power whatever its size, so noise grows as sqrt(K);
  270-984 children of ~40-160 particles at snr 0.0012 (12 kDa arms) are far below the per-child floor at 15 A; the within-view
  contrast that survives sits on the densest region (scale/defocus/registration), not on peripheral motion. The 2026-10-17 small-domain
  result was about sub-class LABELS, not sub-class HALF MAPS.
- Cost: 7.1 (r131) / 20.3 (arms) core-h for the sub-class route vs 1.06 for the tree's cheap step.

**Why:** closes the last untested use of cls_expansion in the tree (embedding failed 2026-10-21, variance map fails here).
**How to apply:** do not wire sub-class children into cheap_children; the mask's fix is the null-calibrated core (audit 6.1, math-reworks
session), applied to the cheap children's variance, which is already localized. Related: [[subclass-node-embedding-2026-10-21]],
[[pooled-2d-route-sees-small-domains-synthetic-2026-10-17]], [[hetero-tree-algorithmic-audit-2026-11-01]].
