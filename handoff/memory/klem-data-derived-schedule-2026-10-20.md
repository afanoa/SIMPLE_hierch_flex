---
name: klem-data-derived-schedule-2026-10-20
description: "user 2026-10-20: derive the frequency march from the gold-standard data, not 20/14/12 by depth; node FSC gets COARSER with depth (fewer per state); refine3D_states re-classifies states (labels discarded by dtree_v13); filt_mode=uniform = FSC-driven matching lp; dtree_v14 sixc8d/sixc8s"
metadata:
  type: project
---

User: "can't we derive a non-arbitrary schedule now that we have gold standard trees" and: on real RNase P, refinement
separates states in the manual workflow (test it and see).
- refine3D_states pose_policy=local RE-CLASSIFIES states (prob_neigh over states): state overlap 0.62 at iteration 1
  -> 0.97 at 6 (not converged); dtree_v13 used only its poses. Its log has per-state FSC 0.5/0.143 per iteration.
- Measured node resolution (final iteration, mean over states): real gs2 r.1 (16k/state) 10.5-10.8 / 8.1-8.5 A;
  synthetic depth 1 (4k/state) 12-14 / 9.7-10; depth 2 (1.7k/state) 14.5 / 10.2 -> data say COARSER deeper (fewer per
  state), the opposite of the hand march 20 -> 14 -> 12. Stage-2 matching lp 8.08 A was finer than the FSC 0.143 (9.7).
- SIMPLE: filt_mode=uniform with lpstart/lpstop = matching lp estimated each iteration from the even/odd maps
  (estimate_lplim3D, simple_matcher_refvol_utils.f90) -> lpstart/lpstop become bounds.
- dtree_v14.sbatch: node_refine2.sbatch (filt_mode=uniform, bounds LPSHI 25 .. LPSTOP 6; writes refined/res.txt
  "r05 r143 lp_last" and refined/states_refined.txt), dtree_board6.sbatch (BANDRES=r05:r143 of the same half's parent
  refinement -> lps r05 x3, sqrt x3, r143 x3; else the depth default), CANDSRC=refined: split-test candidates = the
  refinement's re-classified states (>= MINGRP). Runs 02:00: sixc8d (54012860, board candidates), sixc8s (54012861,
  refined-state candidates), root reused from sixc8 (so depth-1 boards keep the depth default band).
- v14 runs (sixc8d/sixc8s) INVALID: refine3D_states rejects filt_mode=uniform (fsc|nonuniform_lpset|none) and the B
  project was corrupted -> every refinement fell back. v15 (dtree_v15.sbatch, dtree_board7 private project copies,
  node_refine3: default NU filter, lpstart = parent's FSC 0.5, lpstop bound 6): works; NU matching handoff picks 6.0 or
  8.1 A; depth-1 FSC 0.5/0.143 12.2-17.3 / 8.7-10.9 A, depth 2 13.8-16.1 / 9.6-10.8. Depth-2 boards ran 12.2 -> 10.3
  -> 8.7 A (data) instead of 12 -> 8. RESULT sixc8d2: 5 states, NMI 0.31/0.28, closed 0.51/0.40, arm 0.56, poses 3.84
  = hand schedule within run-to-run noise (+-0.03): the schedule is not the lever here (depth-2 splits mostly K* 1).
- sixc8s2 (v15 + CANDSRC=refined): leaves BETTER (NMI 0.27/0.24 vs 0.23/0.21 with board candidates), final 5 states
  0.32/0.29, closed 0.53/0.46, arm 0.61 = same band as every variant (final pass evens out); depth-2 refined states
  mostly < MINGRP 1500 -> no candidates (K* 1); one depth-2 refinement (3.7k, K 4) failed.
Related: [[klem-dtree-per-node-refinement-2026-10-20]], [[f692e-validation-cross-half-2026-10-19]].
