---
name: tree-good-pose-null-fails-halves-and-seeding-2026-10-30
description: "KEY open issue: at GOOD fixed poses (0.47 deg / 0.31 px) hetero_tree splits the micnull null (z 10, merged-ref stage B); SIMPLE_TREE_HALVES (own e/o half maps) cuts it to z 4.8; the residual comes from seeds = flex_pca poles estimated on BOTH halves (open design issue, whole tree)"
metadata:
  node_type: memory
  type: project
  originSessionId: c82ee840-ba2c-4b73-a246-e79b0fb2c726
  modified: 2026-10-30T04:25:36.480Z
---

2026-10-30, micelle session (fixture trees in hetero_tree_20261028/fixtures/pfcrt_micelle/trees/*_a1s025_*, *_in_halves).
- fixposes + order roots at 0.47 deg / 0.31 px: micnull K*=2 z 10.10 guard 0.854 (null FAIL, groups track nothing); mic3 z 11.35, NMI
  with state 0.003. At 1.41 deg / 0.92 px the same roots passed (K*=1): sharp references at good poses make the merged-map stage B
  ([[refine3d-lp-set-scores-merged-map-trap]]) lock noise both halves share, and the split test certifies it.
- SIMPLE_TREE_HALVES=1 (worktree tree-micelle-synth, read_mask_filter_refvols: own e/o half maps when lp is set, fixposes only, logs
  '>>> TREE HALVES'; build recipe17d_micsup_fph, gate 15/15): micnull 0.47 deg z 4.83 (still splits, groups track nothing), micnull
  1.41 deg passes, mic3 0.47 deg z 18 now partly by state (NMI 0.034, Cramer V 0.24, leaf angles 3.9 / 8.8 deg vs truth 4.5 / 7.4).
- Residual mechanism (OPEN, whole tree, no fix yet): the root's seeds are flex_pca poles estimated on all particles, i.e. one
  direction from both halves; splitting each half along a shared direction makes group-difference maps agree between halves even on
  noise, and the split test's label-shuffle null cannot see it. A fix would seed each half separately (or match per-half seeds).
- The (sparse) null at mixB poses (1/3 random) splits aligned vs misaligned (z 19-24); a misaligned leaf map can fit a fake rocking
  angle (15.4 deg at protein corr 0.34): report protein corr and pose agreement before any angle.

**Why:** the user asked whether PfCRT states can be pulled out at good poses (EMPIAR-10330 deposited orientations).
**How to apply:** at good poses read a tree split only if it clearly beats the residual null (z ~5 here) and survives nuisance checks;
SIMPLE_TREE_HALVES needs the full gate suite (r131 closed corner, both nulls, qtoy2) before any default. Related:
[[micelle-fixture-pose-scatter-verdict-2026-10-29]], [[soft-multistate-fixed-pose-verdict-2026-10-30]].
