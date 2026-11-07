---
name: soft-multistate-fixed-pose-verdict-2026-10-30
description: "Soft multi-state reconstruction at fixed poses (SIMPLE_SOFT_STATES, scratch builds _soft/_soft2/_soft3) fails on the micelle fixture: beta 1 = effectively hard, noise-locks, nulls fail; true cross-half + difference shrink only delays the runaway (4 iterations)"
metadata:
  node_type: memory
  type: project
  originSessionId: c82ee840-ba2c-4b73-a246-e79b0fb2c726
  modified: 2026-10-30T03:42:24.257Z
---

2026-10-29/30, micelle session (W = hetero_tree_20261028/fixtures/pfcrt_micelle/soft; scratch builds recipe17d_micsup_soft*,
never on hetero-tree). Soft EM inside refine3D refine=cont fixposes=yes: per-state log-evidence -L P/4 (Cartesian T_i),
w = pi exp(beta llk), pi = previous mean weights, sigma2 from the weight-averaged residual, all states gridded with wtab.
- beta 1, merged-map scoring (lp set, [[refine3d-lp-set-scores-merged-map-trap]]): identical to hard beta 100; mic3 NMI 0.147 -> 0.010
  over 8 its; 87 % of weights > 0.99 at accuracy 0.45; maps fake states (fitted 3/7.6/12 deg, members random); micnull grows confident
  splits from every start; random starts never collapse. The unchanged split test ACCEPTED those groups (micnull z 22, guard 1.025).
- True cross-half (soft3, start halves from reconstruct3D): mic3 peaks NMI 0.11 at it 2 then decays; micnull still noise-locks.
  The split test REJECTS these groups (K* 1). SIMPLE_SOFT_SHRINK (state difference x max(0, its e/o FSC)) holds the null 4
  iterations (gains <= 0.04) then the cross-half loop runs away (gains 0.4-0.6 by it 4, maxw 0.94 by it 8).
- Root cause: point-estimate state maps from ~9k particles carry noise the likelihood treats as signal (gaps 8-14 nats vs ~1 ideal);
  even cross-half, sorting one half by the other's map noise makes that noise reproducible.

**Why:** the coordinator asked (2026-10-29) for a soft multi-state prototype to find the 8/16 deg helix states at held poses.
**How to apply:** do not revive soft/posterior state weighting at fixed poses without a map-noise model that survives the cross-half
loop; any half-map statistic computed after classification that used the merged map is not independent. Related:
[[micelle-fixture-pose-scatter-verdict-2026-10-29]], [[pfcrt-real-poses-do-not-reproduce-2026-10-30]].
