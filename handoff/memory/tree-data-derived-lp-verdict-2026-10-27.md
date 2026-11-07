---
name: tree-data-derived-lp-verdict-2026-10-27
description: "hetero_tree node lp from the node mean's FSC=0.5 (test switch SIMPLE_TREE_LP=fsc) tested 2026-10-27: no root gain, r1211 closed recall .94->.78, deep nodes coarsen into the split re-test band; keep fixed 15/10"
metadata:
  node_type: memory
  type: project
  originSessionId: 1e61b3a2-c9de-41d3-b849-157ff77baec6
  modified: 2026-10-28T03:43:42.571Z
---

Test (worktree jovial-jennings-821b3d, uncommitted; cluster tree_data_lp_20261027): stage B lp = node-mean FSC=0.5, A/flex 1.5x, box by autoscale.
- Roots: picks 7.0-9.8 Å. qtoy2 unchanged. r1211 closed stage-B recall .937 -> .776: the 2.4 % flex-lp change re-seeded it (seeds ARI -0.02).
  qtoy null: same verdict as off under both split tests.
- Saved m33 692k nodes (offline, recB halves): picks coarsen with size (5-10k: 13.8 Å; <5k: 14.6, up to 18 Å; flex/A 18-27 Å).
  103 of 122 nodes fall inside or coarser than the 14.7-10.2 Å re-test band.
- r2121212 (lost real split): on arm made a split unrelated to m33's (ARI 0.00) that passes even the rho floor.
- TRAP: ~/SIMPLE_smr4/build/bin/simple_exec predates its SIMPLE_SPLIT_RHOMIN source (no RELATIVE EFFECT line); the agreed split
  setting (contrast=maps, rhomin 0.03) needs ~/SIMPLE_smr5o (the qtoynull_fixm build).

**Why:** FSC-driven limits coarsen exactly where real RNase P differences live (<= 15 Å at small nodes).
**How to apply:** do not propose data-derived per-node lp of this form again without addressing node-size coarsening.
Related: [[rnasep-closed-state-must-be-retained]], [[tree-design-feature-freeze-2026-10-25]], [[beegfs-transient-write-failure-2026-10-27]].
- REOPENED 2026-10-29 by the user: "certain molecules have heterogeneity appear at different resolutions, so you would think a data
  derived schedule would make more sense". Next form: drive the CLASSIFICATION band by the between-state difference-SNR spectrum
  (half maps; SNR_d > 1, capped by consensus FSC), keep alignment on the consensus limit - replay first (data-derived-resolution
  session), not the consensus-FSC form above.
- 2026-10-29 DIFFERENCE-SNR REPLAY (332 saved nodes, beegfs hetero_tree_20261028/dsnr_band/): heterogeneity lives at 20-40 A in
  every dataset (stage-B I50 26-34 A, I90 17-25 A; RNase P closed incl. the user's one-helix pair reach 90 % info by ~16.7 A);
  datasets differ only ~15-25 A. Rule (finest shell with pair SNR_d >= 1, consensus-capped) keeps >= 90 % info for every named
  real split; still coarsens with node size. KEY: the split test's re-test band 14.7-10.2 A holds 32 % of NULL difference info but
  only 3-19 % of real splits' -> re-test looks in the wrong band (null false positives + r2121212 power). Coarser re-test replay
  -> adaptive; stage-B band fixture replays (~150 core-h, recipe17b) -> data-derived-resolution.
