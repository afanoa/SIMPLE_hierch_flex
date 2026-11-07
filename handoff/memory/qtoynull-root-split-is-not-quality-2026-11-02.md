---
name: qtoynull-root-split-is-not-quality-2026-11-02
description: CORRECTION - the qtoynull root's accepted (pre-guard) split does not separate snrq or B-factor (d 0.03 / 0.01 by the truth); it is a non-quality false positive the locality guard merges; the stage-B corr IS a good snrq proxy (r 0.66-0.70)
metadata:
  type: project
---
2026-11-02 (math reworks, flag2/ standalone split replays on base_qtoynull root, build 9005a8c77 line): pre-guard groups 31139 / 28861,
truth snrq 1.074 vs 1.050 (within sd 0.38), bfac 100.4 vs 98.9, junk none -> no quality separation; guard L/L_ref 0.225 merges.
Per-particle proxies in state_mergeup: stage-B corr (qv) correlates 0.66-0.70 with snrq on qtoynull / qtoy2; fit_contrast's contrast
0.22-0.55; the whitening noise variance 0.00; p_empty is 0.0000 for every particle on every fixture (qtoy junk has contrast).
**How to apply:** do not expect a quality flag to fire on qtoynull; to test a quality-split detector, build a fixture whose split
follows snrq (or check a node's groups against the truth snrq first). Supersedes the "quality split" reading in
[[split-test-accepts-quality-split-on-null-2026-10-27]] for this fixture. Related: [[hetero-tree-math-reworks-setup-2026-11-01]].
