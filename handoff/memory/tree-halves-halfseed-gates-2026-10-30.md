---
name: tree-halves-halfseed-gates-2026-10-30
description: "SIMPLE_TREE_HALVES + SIMPLE_TREE_HALFSEED (diff 4f7724ed on 21fac34db): all nulls pass, qtoy2 better, but r131 loses the closed corner (children z ~1) and OLE narrows; not a default"
metadata:
  type: project
---

2026-10-30 gold-standard node fix (adaptive-split-merge session), beegfs hetero_tree_20261028/route9/halves (v2 build, diff
halves_halfseed_r17d_v2.diff md5 4f7724ed, off = 21fac34db). Both switches on, re-align (no fixposes in this base):
- nulls one leaf: micnull 0.47 deg z -0.91 (off 61.7 guard-merged; halves alone -0.32, seed alone 52.3), qtoynull, olenull.
- qtoy2 md2 9 leaves NMI .712 vs .686, c18 F1 .89. r131 md1: root split kept, both children z ~1.1 vs 3.71 (base 6.1 / 5.3) ->
  closed corner .88/.52 vs .97/.54: HARD GATE FAIL ([[rnasep-closed-state-must-be-retained]]). OLE md2: 6 leaves, hs covers the r14
  landscape only .11.
- Mechanism: with independent half maps the split test's beyond-EM pair statistic (shells past the EM band) loses the
  shared-reference correlation that carried depth-1 splits; the micnull confound (re-alignment scatters poses to ~7.5 deg) means
  good-pose seeding stays untested ([[tree-good-pose-null-fails-halves-and-seeding-2026-10-30]]).

**Why:** the user asked for true e/o halves + per-half seeding gated on nulls, qtoy2, r131 corner and OLE.
**How to apply:** do not default either switch; if pursued, first attribute the r131 loss (halves-only vs seed-only arms, ~25 core-h
each) and recalibrate the split test band/null under halves before re-gating.
