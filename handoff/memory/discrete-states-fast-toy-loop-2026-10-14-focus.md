---
name: discrete-states-focused-route-2026-10-14
description: "Fully automatic discrete-state route validated at 10 deg (13 % minority, SNR .03): focused flex (pcg_mskfile = moving-region mask) -> gmm placement (TV 0.097, 0.951) -> masked-deviation E-step, fixed placement fractions, x0.5 scale -> 0.963 flat over 5 its; whole-map flex fails <= 10 deg. RNase P focused runs per region launched"
metadata:
  type: project
---

2026-10-14 ~00:45 (reseed26, ~/SIMPLE_master1013 uncommitted). Recipe that works where the whole-map route fails:
1. flex_pca ... pcg_mskfile=<soft mask of the moving region> lp=10 neigs=8 npreimages=2 state_placement=gmm
   (focused basis: latent d' 4 deg 0.20->0.64, 6 deg 0.22->1.26; 10 deg placement view-balanced TV 0.097, acc 0.951).
2. fractions = the placement's mixing proportions (SIMPLE_SCORE_PIS), not EM.
3. E-step loop: reconstruct3D rec_backend=pcg state_refs=deviation m_estimator=flex with SIMPLE_DEV_MSK=<same mask>
   + score_states (true cross-half, reseed25 fix) with SIMPLE_SCORE_TEMP=0.5. 10 deg: 0.9635 0.9636 0.9630 0.9628
   0.9627; maps: state 2 r 0.997 to the truth-label map, difference r 0.84 at 0.75 amplitude (confined to the mask).
Without the mask the 10 deg loop drifts (0.983 -> 0.962 fixed fractions; 0.985 -> 0.939 EM); at 6 deg masked decays
slowly (0.911 -> 0.880), at 4 deg the minority dies. Scripts toyexp/focusloop.sbatch, estep_loop7.sbatch (SEEDDIR).
RNase P: masks rnasep_regions/region{1,2,3}_focus.mrc (G10cov_B3 cores > 0.5, dilated 7.3 A, cosine edge); runs
rnasep_coupled/focus/region{1,2,3} (every third refined particle, ~230k; focus_rnasep.sbatch; jobs 53853314-6).
Related: [[discrete-states-fast-toy-loop-2026-10-13]], [[multibody-flex-means-region-pca]].
