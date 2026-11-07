---
name: split-test-single-draw-variance-2026-11-02
description: KEY - marginal split verdicts (r131 r1 corner node) are set by the one random equal-N subsample (old T .074-.166 over 4 seeds, 1 PASS in 4); FIXED by the draw-free pair statistic (full-N maps, analytic noise at the common n), which is seed-identical and 3x tighter
metadata:
  type: project
---
2026-11-02, math reworks (hetero_tree_20261028/math_reworks, note hetero_tree_reworks_2026_11.md sec. 8): standalone state_split
replays at SIMPLE_SPLIT_SEED 1-4 (instrumented build_df) show the production pair statistic spreading sd .02-.04 at marginal nodes:
r131 r1 (seed-2 tree, the closed corner) .074 .137 .166 .097 = 1 PASS in 4 (the tree's own draw .154 passed by luck); variant-D r131
r1 seed 1 0/4; qtoy2 r3 2/4, r5 1/4. The draw-free statistic (per-shell cross-half product of the full-N difference maps over that
plus the groups' half-difference noise scaled to the common n; null = same at full N from within-strata shuffles) is identical at
every seed, sits at the old T's mean (.117 vs .119) and has a null a third as wide, so it passes every one of those nodes at z 4.6-11.
**Why:** a cross-half FSC of ONE seeded random subsample has draw variance as large as the margins at ~20k-particle nodes.
**How to apply:** gate split-test changes with the draw-free statistic (state_mergeup test_once, commit 9ee5d5c72), never with one
draw; its extra power adds ~1 marginal split per tree (qtoy2 within combo 28, NMI .764 -> .742) - the user decides the default.
Related: [[math-reworks-verdicts-2026-11-02]], [[rnasep-closed-state-must-be-retained]].
