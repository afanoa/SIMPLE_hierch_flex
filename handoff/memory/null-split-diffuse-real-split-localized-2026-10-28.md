---
name: null-split-diffuse-real-split-localized-2026-10-28
description: "KEY: the qtoy2/OLE null splits differ between their stage-B maps DIFFUSELY (top-5 % share L 0.13-0.19, L/L_ref 0.2-0.3); every real split is LOCALIZED (0.45-0.96, ratio 0.7-1.5) incl. the m33 closed pair at the null's rho; corrected correlation and dB don't separate"
metadata:
  type: project
---

2026-10-28, adaptive split/merge session, map_sep.py (scripts/hetero_tree_math_experiments/split_merge_replay), outputs
/mnt/beegfs/elmlund/afan/hetero_tree_20261028/route9/sep.
- L = share of the envelope-matched, noise-subtracted residual energy (two stage-B state maps from even/odd unfil halves,
  40-10 A, soft sphere, sigma-3 smoothing) in the top 5 % of voxels; L_ref = the same for the consensus density.
- Nulls: recipe14 qtoy2-matched root 0.19 (0.3), OLE-matched root 0.13 (0.2).
- Real: route-A r111 0.45 (0.7), m33 closed pair 0.71, r1211 closed/open child 0.85, qtoy2 root 0.73-0.84, qtoy2_bp
  closed arm 0.91 / junk 0.82, OLE r1 0.87, build-7 closed lineage 34 pairs 0.74-0.96. Robust to q 2-10 %, sigma 2-4.
- Not separating: attenuation-corrected correlation D, the B-factor difference.

**Why:** the null difference spreads over the whole molecule; real states concentrate it in a region.

**How to apply:** a guard on the stage-B half maps (refuse/merge at L/L_ref < 0.5) is the candidate fix for the failing
null gate. Confirm on more null seeds first (only 2 independent nulls; recipe13 = recipe14 stage B).
Related: [[split-test-accepts-quality-split-on-null-2026-10-27]], [[null-guards-do-not-stop-692k-deep-splits-2026-10-28]],
[[rnasep-closed-state-must-be-retained]].
- CONFIRMED + BUILT (same day): null seeds 2-4 (recipe14, production env, root only) all six roots L/L_ref 0.2. Guard in state_split
  (locality_guard_on_r14.diff, md5 1653fc1b; default on, SIMPLE_SPLIT_LOCALITY=off disables; reads ../stageB halves; hard 40-10 A
  band - the soft bp() band put route-A r111 at 0.54 instead of 0.66). Handed to the factorial session for recipe15.
- CONFIRMED 2026-10-28 ~20:30 on 6 more null seeds (qtoynull/olenull s2-s4, recipe14 root only): ratio 0.2 everywhere (seeds 2-4 did
  not even split; a qtoynull s3 replay split was merged by the guard at 0.203); seed-1 qtoynull split 0.274 -> merged. Fortran guard
  (locality_guard_on_r14.diff, state_mergeup + exec_state_split, default ON) keeps every positive: qtoy2 root 0.97-1.08, route A r111
  0.66 (weakest), r1211 child 1.24-1.25, qtoy2_bp closed arm 1.26 / junk 1.37, OLE r1 1.48. Bar 0.5. Built as recipe16 = recipe14 + guard.
- 2026-10-28 ~22:00 RECIPE16 (recipe14 + guard) gates_env: qtoynull PASS (split test K*=2, guard 0.205 -> merged, K* 1, one group);
  olenull PASS (K*=1); qtoy2 root 4 states kept (guard 0.965-1.082); r111 route A kept (0.660). FIRST BUILD WITH BOTH NULLS PASSING.
- 2026-10-29 COARSE RE-TEST LEAKS: a 40-20 A re-test inside the cross-half EM band (EM shells 2..kb, each half labelled from its OWN
  data against the other half's maps) is NOT held out: T ~0.92-1.00 for nulls and positives alike, ALL 8 null seeds K*=2 (shuffle
  null 0.02-0.24). Production (~15-10 A, outside the EM band) is clean. Proposed fix being built (adaptive, env switch): fit the
  K-group model on ONE half in 40-20 A, score the OTHER half's per-particle log-likelihood gain vs consensus, shuffle null within
  view x defocus strata; band may then coincide with the fit band. Bars: nulls K*=1 without guard, positives kept, more power on
  m33 closed pair (prod z 3.66 vs 3.71) and r2121212 (prod 0.98).
- RECIPE16 FINAL (2026-10-29): guard changes nothing on any real gate (qtoy2 15 leaves, r111 7, r1211 5+5 identical to recipe14;
  real ratios 0.72-1.37) and fixes qtoynull (one 60k leaf; 55 min / 19.9 core-h vs recipe14's false split 117 min / 35.8). qtoy2
  leaves NMI(structure) 0.685, NMI(quality) 0.001.
