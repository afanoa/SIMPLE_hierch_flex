---
name: blind-greedy-realign-fails-rnasep-2026-11-03
description: "TEST 2026-11-03: one refine=greedy pass (and greedy + neigh x5) does NOT recover poses of held-out RNase P particles against a fixed consensus (median 120 deg, 4 % within 5 deg); scores and half-map FSC cannot detect it; adding 240k greedy-aligned discards drops the mean 4.9 -> 7.8 A"
metadata:
  node_type: memory
  type: project
  originSessionId: a0ebc4a3-03d0-403c-a0f5-633af1f0fc21
  modified: 2026-11-03T17:43:18.595Z
---

User asked (2026-11-03): align thrown-away / extra-dataset particles to a nice volume. Is one greedy pass enough, and should the tree use the consensus volume or a new reconstruction? Dir: batch2 `/mnt/beegfs/elmlund/afan/realign_test` (realign.sbatch, posetool.py, mapcmp.py, rowcopy.py, ctrl_summary.txt, disc_summary.txt). Build: math_reworks/build_lf. Keys: lp=8 trs=10 volrec=no pgrp=c1 mskdiam=232.

- **Control**: 20k of the 692k kept particles, poses zeroed, reference = map of the other 672k:
  - greedy x1: median 120 deg off, 3.9 % within 5 deg;
  - greedy + neigh x5: median 110 deg, 13.9 % within 5 deg;
  - neigh from the true poses: median 4.7 deg, so the reference and frame are fine. The global search is what fails, as on 10180 and PfCRT ([[refine3d-shc-virgin-field-trap]]).
- Wrong poses score HIGHER (corr 0.378 vs 0.372). The 20k map from greedy poses reaches FSC 0.143 at 8.08 A, which is exactly lp. Maps from the true poses reach 13 A. That is reference bias, so FSC and score cannot validate realignment.
- **Discarded 240k**, greedy x1 vs the 692k consensus, then all 932k reconstructed: FSC 0.143 4.88 -> 7.81 A, corr 0.91 to the 692k map.
- TRAP: discarded rows carry ptcl3D col 48 = 0 (a weight). Copying only poses and state left them out of the reconstruction. Copy whole records (rowcopy.py).
- The tree takes no vol1. Its root rebuilds the mean (S_MEAN reconstruct3D) from the input poses.
- Untested: alignment inside a refinement that rebuilds the reference, with the 20k blind control kept as the yardstick.

**Joint refinement (refine3D_auto, job 54921399, `realign_test/joint`).** 100k particles: the 20k control starting from their greedy x1 poses, plus 80k kept particles with consensus poses.
- Snapshot at iteration 12 of 50: the 20k are a median 115 deg off, 8.6 % within 5 deg (3.9 % at the start). The main loop is a local prob_neigh search, so it cannot undo errors that large.
- The registration pass (greedy over all particles) moved 42 % of directions by more than 4 deg. Resolution went 4.6 -> 6.9 A, then partly recovered to 5.7-5.9 A.
- The 80k drift a median 4.8 deg; 17 % move more than 30 deg. The neigh-from-true control also had 16.5 % beyond 30 deg, so about a sixth of the consensus poses are unstable under any re-search.
- Iterations sample only ~3.6k particles each (about 5.5 min per iteration).
