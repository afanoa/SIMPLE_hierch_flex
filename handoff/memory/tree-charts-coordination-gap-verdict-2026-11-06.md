---
name: tree-charts-coordination-gap-verdict-2026-11-06
description: "VERDICT 2026-11-06: hetero_tree as coordinated local PPCA charts (child->parent affine maps on fflex outputs) gives a per-particle global coordinate with a modest gain over the root chart (qtoyflex s1 ceiling .59 vs .53, closed/holo AUROC .82/.92 vs .78/.87; qtoy2 closed d' 4.7 vs 3.0); the XD density-gap split criterion has no false gaps but is blind below chart d' ~2.5 (misses closed/holo on qtoyflex and the r131 closed split): annotation only, no Fortran modelling layer"
metadata:
  node_type: memory
  type: project
  originSessionId: 7da7e826-c848-444f-a66e-ad78379bb739
  modified: 2026-11-07T00:49:46.894Z
---

Task (paper coordinator, 2026-11-06): the tree as a hierarchical mixture of local PPCA charts with LLC-style global coordination, and an
LLE-style density-gap split criterion. Python only: scripts/hetero_tree_mest/mest_charts.py (worktree mest-formulations; cluster copy
math_reworks/mest, runs/charts_<tree>/charts.txt). Full tables in the side-session note mest_formulations_status.md.

Facts about the chart files: flex_pca_coordinates.txt is the delivered MAP z after the gauge rotation; flex_pca_embedding.bin holds the
raw MAP z, per-particle posterior precisions (k x k), and a trailing block with flex's calibrated noise scale and deconvolved z (layout in
simple_flex_pca_embedding_io.f90). The LS map raw -> delivered is a rotation only approximately (|R^T R - I| .7-.9); the variance budget
that balances is the un-shrunk z (MAP x P^-1 (P + Pi) with flex_pca_map_prior.txt) with covariance ns * P^-1.

- Coordination: every node already is a chart; composing precision-weighted child -> parent affine maps gives one global frame. On the
  30-conformer continuum (qtoyflex) the composed coordinate tracks s1 at .57-.59 vs .52-.53 for the root chart and orders conformer means
  at .78; within nodes no chart exceeds .35 (noise swamps adjacent conformers). On qtoy2 the local charts resolve discrete sub-states the
  root chart blurs (within-node s1 ceiling .50/.58 vs .40/.42; closed d' 4.7 vs 3.0, holo 6.6 vs 3.7) and the composition keeps it.
- Gap criterion (1-D XD along the split direction, D between the children's modes, uniform-continuum null): qtoy2 root 8/10 pairs
  (D .5-.86, p .01, chart d' 2.8-5.3), r4 closed-vs-22 borderline (D .42, p .33, d' 3.5), misses at d' <= 2.6; qtoyflex all D 0 (arm cuts
  correctly, closed/holo wrongly: chart d' 1.3/2.3); nulls and r131 D 0. The uniform null's own D is .3-.5 at K* >= 4 (mixture ripple).

**How to apply:** do not build a chart-mixture modelling layer in Fortran; if the paper wants a particle-level landscape, port only the
output (composed global coordinates per particle, plus per-split D/p/d' as annotations of discrete vs continuum cuts). The gap test cannot
replace the cross-half split test. Related: [[mest-formulations-verdict-2026-11-06]], [[discrete-splits-capture-the-span-states-that-matter]],
[[hetero-tree-paper-plan-2026-11-06]].
