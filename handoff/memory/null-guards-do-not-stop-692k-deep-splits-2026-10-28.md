---
name: null-guards-do-not-stop-692k-deep-splits-2026-10-28
description: "build-7 692k trees ran NO split guard (old max-of-20); replaying 12 deep K*=2 nodes under maps+0.03 floor (+core qnull, +pairs, +node-local count) keeps 11-12 of 12; deep rho 0.07-0.15 > closed pair 0.047 -> null guards cannot shorten the tree"
metadata:
  type: project
---

2026-10-28, adaptive split/merge session (battery6, /scratch/cluster_scratch/atifao/qnull_port/battery6, qnull_table6.py).
- TRAP: site_states/build_sst7_it3crop_emb (route-A branch on 8c1b66f31) has no split-test port. Its trees split with contrast
  unset, threshold = max of 20 shuffles, no z test, no floor. Check a tree's split.log for 'z_crit' before reading its K*.
- Replay of base r3222/r3221/r3122/r2321/r222/r3121/r232/r111/r21 + route-A r322/r332/r522, seed 1:
  - a (maps + rho 0.03 @ 460): all keep K*=2, z 4.1-20, rho 0.036-0.20.
  - b (+ QNULL=core + GUARD=pairs) and c (b @ 2N/3000-1): refuse only r2321 (z_q 1.46/1.52).
  - m33 closed pair r1111221 kept in every arm (rho 0.047, z_q 2.18/1.90).
  - build-7 closed splits (factorial scorer: base r11 -> r112, route A r11 K*=4) kept in every arm at z 12.7-19,
    z_q 7.9-13. The deep nodes sit at z 4.1-10, z_q 1.5-4.7: a z_q ~5 bar separates them here but refuses the m33 pair, and
    the closed purification chain below r112 is untested.
- The deep splits beat the closed pair on rho (1.5-3x) and mostly on z_q, so a null guard that stops them drops closed first.
  Runtime needs a non-null stop rule (depth/minpop/budget).
Related: [[quality-matched-null-core-score-2026-10-28]], [[rnasep-closed-state-must-be-retained]],
[[tree-runtime-goal-5-7h-at-700k]], [[contrast-switch-kills-real-rnasep-splits-2026-10-27]].
- TOP-DOWN z_q >= 5 walk (battery7, qnull_walk7.py), base tree at 8 h: 27 pass, 18 refused -> leaves 11-52k (median 23k,
  60 % of particles); 88 of 145 nodes saved. Both closed lineages pass (r11 9.69, r112 7.88, r1121 8.36; r12 7.71, r121 7.66,
  r1211 11.70, r12112 6.23; lowest closed 6.23 vs highest refused 4.94) but no finished node is at the closed corner yet.
  m33 pair parent r111122 fails z at 460 under maps (3.45) -> the within-closed conf19/18 pair is lost by the operating
  point itself, not by the bar.
- z_q >= 5 top-down FAILS on known-answer trees (battery8, qnull_walk8.py): qtoy2_bp stops r12 (2.08) above the closed-arm split
  r1211 (12.82) -> c18 purity 0.93 -> 0.27; OLE recipe4 17 equiv units -> 8 (small states + all 4 distinct sibling pairs
  folded under r21 2.27 / r122 / r1122). z_q is NOT monotone down a lineage: never use a top-down quality-margin stop.
