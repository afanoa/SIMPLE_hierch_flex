---
name: igm-choose-vs-em-replays-2026-10-28
description: IgM toy (armsigm1) replays of 11 m69M split nodes: choose delivers stage-B at 8/11, sharper single-arm motion at 5/7 (+.01-.17 eta^2), worse at r121; junk node +28 % junk; 0/20 arm extremes at any single split
metadata:
  type: project
---

2026-10-28 replays (igm_r14_vs_r12b/replay, igm_replay.sh, recipe14 state_split on the old m69M tree's saved stage-B, maps + rho 0.03,
132 tests, 20-120 s per arm; truth = factorial/armsigm1_truth/arms_truth.txt, scorer fixtures/scripts/igm_score.py):
- choose took the stage-B states at 8 of 11 nodes (margin > 3 SE); r, r22, r112 stayed EM (identical).
- Dominant arm eta^2 sharper at r1 .418->.445, r2 (EM mixed a2+a4; choose a4_th .103->.271), r11 .216->.260, r12 .407->.420,
  r122 .293->.342; equal at r111 (a5_th .155->.195); WORSE at r121 (a3_ph .211->.173). Other-arm leakage down at 6 / 7.
- Junk split r21: choose's junk group 9,675 vs EM 7,571, both pure.
- Extremes (10 % angle tails, group >= 3000 at purity >= 0.5): 0/20 at every single split, both arms.

**Why:** first test of choose on a continuous-motion system (IgM has no discrete conformers: 201 combos x 5 arm angles + junk).

**How to apply:** score IgM by eta^2 / tail purity / junk, not per-conformer F1 (combo ~800 < minpop); a node-level replay is minutes.
See [[choose-routing-fixes-closed-dilution-2026-10-28]], [[fast-loop-tests-size-to-the-question-2026-10-28]].
- FINAL IgM full trees (200k armsigm1, maxdepth 4, 2026-10-29): recipe14 (route A + choose) 83 leaves, NMI 0.441, 31,369/40,000 junk in
  7 pure groups, 3/20 arm extremes, eta^2 0.29-0.57 on every arm angle; recipe12b 26 leaves, NMI 0.334, 26,053 junk in 1 group, 0/20,
  eta^2 0.11-0.46; old m69M 58 leaves NMI 0.343. Route A + choose clearly better on IgM.
