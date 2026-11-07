---
name: pfcrt-real-poses-do-not-reproduce-2026-10-30
description: "KEY: real PfCRT poses do not reproduce between refinements of the same particles (pfcrt_core 2_refine3D vs 3_refine3D: angle diff quartiles 93/129/156 deg, 96 % > 30 deg; vs 13_solve3D medians 126-131; random ~126). Their 4.9-5.8 A FSCs come from lp-set merged-reference refinement, not gold-standard"
metadata:
  node_type: memory
  type: project
  originSessionId: 9643ce3d-fc7a-4dc1-860a-59db10867e7c
  modified: 2026-10-30T03:32:43.449Z
---

2026-10-30, micelle session, pose_compare.py (hetero_tree_20261028/fixtures/pfcrt_micelle/scripts; validated: fixture vs truth 1.41 deg,
self 0), same particles row for row, best global frame rotation removed.
- pfcrt_core 2_refine3D (neigh, masked PCG, 6.1/4.9 A) vs 3_refine3D (greedy full box, 6.6/5.8 A): frame offset 8.8 deg; per-particle
  angle difference quartiles 93/129/156 deg, 99.9 % > 5 deg, 96 % > 30 deg. Either vs dibyajyoti's 13_solve3D (8 A): medians 126-131 deg.
  Uniform random rotations give ~126 deg: the runs agree no better than chance.
- Both passed lp= with the deposited model as reference, so halves were matched against the merged map
  ([[refine3d-lp-set-scores-merged-map-trap]]): their FSC does not measure pose accuracy.
- The decisive test is pose agreement between independent restarts: Hans's 10 PfCRT posterior restarts (pose_marginalization.md,
  4.03 A FSC 0.143, 6.2-7.1 A FSC 0.5) live in /home/elmlundho on the BOX, unreadable by the session.

**Why:** the user (2026-10-30) asked whether PfCRT states can be separated given a high-resolution map and good poses.
**How to apply:** don't calibrate pose accuracy from an lp-set FSC; treat PfCRT tree results as uninterpretable until restart-to-restart
pose agreement is shown. Related: [[micelle-fixture-pose-scatter-verdict-2026-10-29]].
- TRAP (2026-10-30): dibyajyoti's PfCRT_smpd_1p3/picking_map.mrc (= pfcrt_core/ref/emd20806_1p3_box240.mrc) holds 1.3 A content
  but a 1.035 A header (cell 248.4 A for box 240): SIMPLE takes smpd from its parameters, but Chimera / header-reading tools mis-scale it
  by 1.256. Our 16.9k map vs EMD-20806 also differ by an x/z axis swap.
- 16_selection (44,527 = 36 of 50 classes from 15_solve2D on the 56,809) is the user's cleaner set; 8,291 picks match deposited
  EMPIAR-10330 particles (sel16/match.txt). Being aligned against the 16.9k even/odd half maps (runs A/B) for pose validation.
- EMPIAR-10330 VERDICT (2026-10-30): NO STATES. Deposited poses verified gold-standard (3.53/3.30 A in the 6UKJ region, no alignment).
  Six fixposes + SIMPLE_TREE_HALVES roots ({all 16,905 / agree 8,601} x {no mask, 3-class support, 6UKJ}, minpop 250/130): K*=1 or a
  50/50 split at or below the null floor (z 4.83). The one above it (all, 6UKJ, z 5.13) tracks no nuisance and is absent from the agree
  half (z -0.65): noise. Fixed-pose roots are deterministic, so SIMPLE_SEED replicates are uninformative. Matches the user's "rigid".
- 16_selection (2026-10-30): greedy+neigh alignment failed (11 % unmatched within 10 deg); refine3D_auto vs independent 16.9k halves passed
  at the bar (A vs B 40 % unmatched within 10 deg; matched vs deposited 39.7 %; ~2.7 deg per run). Split maps (4,699 each, agreeing,
  no alignment): non-deposited picks uniformly worse (FSC, consensus corr, dB +24 A^2), population-specific difference equal in TM and Fab
  (Fab can't move), belt difference also in the same-population control, rigid fit 0.0 deg -> worse particles, no other state.
  Maps: ~/Downloads/pfcrt_sel16_{deposited_matched,non_deposited}_4699ptcls_*.mrc. Roots on the 19,176 agreeing picks (fixposes + halves,
  minpop 288): K*=1 for no mask (z -0.35), 3-class support (z -0.05) and 6UKJ (z 2.63, z_crit 3.38). 16_selection VERDICT: NO STATES.
