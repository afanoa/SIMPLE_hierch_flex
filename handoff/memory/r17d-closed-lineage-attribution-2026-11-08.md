---
name: r17d-closed-lineage-attribution-2026-11-08
description: "VERDICT: the r17d 692k tree reaches its closed nodes by selection along ONE dimension of each node's own flex latent (r13->r131 z1 .82, r23->r232 z1 .18, 93 % of the parent's closed rows); noise / ice / micrograph / defocus flat at every split; the corner r1311_s2 (1.00 own) comes from a split the reference rows barely track (4.2->6.8 %)"
metadata:
  type: project
---

rnasep-misfit lineage attribution, 2026-11-08 (`handoff/rnasep-misfit/lineage/LINEAGE.md`; `math_reworks/rnasep_misfit/lineage/`,
`lin_attr.py`). Tree `integration/f692_r17d_s1`; reference = the 8,720 closed rows (1.26 % of the root).

- Closed share by depth: r 1.26 % -> r1 1.94 -> r13 4.11 -> r131 6.94 -> r1311 5.57 (r1312 8.34) -> r1311_s2 6.79; r2 1.31 -> r23
  2.92 -> r232 5.71 -> r2322 9.37.
- The concentrating steps are binary splits along a single node latent dimension (LDA .84-.86): r13->r131, r23->r232, r232->r2322.
  K-way levels add "good particles" (SIMPLE corr / amplitude .59-.60, poses that moved less in the parent's stage B).
- Flat at every split (|AUROC - .5| <= .05): the empirical-Bayes noise level (6 bands), defocus, micrograph, FoilHole. Views /
  2D classes .53-.58. Inside nodes the closed rows follow amplitude (.62-.72), not noise.
- Refitting the latent on a subset adds .02-.04 closed AUROC (root latent on the same particles vs own), nothing at depth 4, and
  never beats the root latent's .786 over all 692k: the closed direction is already in the root latent.
- Corner puzzle: r131->r1311 went AGAINST the closed rows; r1311->s2 is a strong split along r1311's z1 with mild enrichment and the
  weakest passing re-test (z 4.87 vs 3.71); nothing nuisance-like follows it -> the 1.00 own readout is the imprint pattern; honest
  fission reads (extremes-repro) decide.
- Node means (pmean.mrc) and the parents' recB / stageB maps are PRUNED in this tree; node projects, rows, flex coordinates, split
  logs survive. Noise-model sparse sites per node need one reconstruction per node first (proposed, not run).

**Why:** the user asked why the closed states appear at the nodes they do; this names the mechanism.
**How to apply:** treat r17d's closed leaves as products of latent-direction selection, not of nuisance separation; judge their
closed content by honest (fission) reads, not own-data maps. Related: [[integration-closed-two-leaves-2026-10-29]],
[[extremes-fission-verdict-2026-11-07]], [[rnasep-misfit-noise-spectrum-verdict-2026-11-08]].
- 2026-11-08 coordinator: approved the per-node rebuild (r13, r1311, r1312, r23, r232, r2322; ~20 core-h). The question: does a
  shared-basis site carry the closed direction where the tree cuts, and could the shared model have made these carves itself?
