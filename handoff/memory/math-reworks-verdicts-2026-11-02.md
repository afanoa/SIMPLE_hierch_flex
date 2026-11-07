---
name: math-reworks-verdicts-2026-11-02
description: verdicts of the hetero_tree math reworks - subspace site test GO (26c350abb, variant D, upstream); single-pass split test NOT default (draw-free statistic gated - corner 2/2, nulls clean, qtoy2 NMI .742 vs .764; user's call); mask z-core NO-GO; SVD/angle partial
metadata:
  type: project
---
2026-11-02, branch claude/tree-math-reworks-30ff30 on claude/tree-dist (cleanup e12be4e8d + variant D 26c350abb); note
doc/implementation_notes/hetero_tree_reworks_2026_11.md; batch2 math_reworks/ (builds build_sub4 = D, build_df = instrumented draw-free):
- Rework 1 (audit 6.3) GO as variant D (upstream 26c350abb): one subspace test on the endpoints' cross-half deviation Grams; pole
  strength = weakest endpoint hold-out z; floor of 6 on every site_rank state. Nulls 1 leaf; qtoy2 NMI .764 vs .712; r131 corner 2/2.
- Rework 2 (6.7) single pass: 9c88430cd (not default), 9ee5d5c72 draw-free statistic, 680bda2d5 caller back on split%test. Gates on
  the draw-free pass: nulls 1 leaf, r131 corner retained both seeds (z 7-11; the production test's seed-2 retention is 1 draw in 4),
  qtoy2 9 leaves NMI .742 (D .764, base .712: one within-combo split at r3), r131 one extra 0.93-closed-arm split per seed, OLE = D.
  Default change = user's decision (one line in exec_state_split: test -> test_once with stats).
- 6.1 z-core NO-GO; 6.2 SVD partial; 6.4 angle partial; 6.5 design only. Item 13: p_empty = 0 on all fixtures; flag in split_stats.txt.
**How to apply:** build on 26c350abb for the site test; gate any split-test change with the draw-free statistic; campaign ~1500 core-h.
Related: [[hetero-tree-math-reworks-setup-2026-11-01]], [[split-test-single-draw-variance-2026-11-02]],
[[site-strength-ranking-decides-r131-corner-2026-11-01]].
