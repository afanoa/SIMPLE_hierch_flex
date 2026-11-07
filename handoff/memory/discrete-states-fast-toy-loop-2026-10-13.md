---
name: discrete-states-fast-toy-loop-2026-10-13
description: "Fast (1-2 min) known-truth loop on toySE latents+images: shared-cov GMM on the particle flex latent >> flex placement and >> any averages-only method; GMM seed + soft cross-half image E-step = 0.97 (SNR .03) / 0.88 (SNR .01); averages-only seeds lose the minority; K selection open"
metadata:
  node_type: memory
  type: project
  originSessionId: 34b66a02-03d2-42df-91f8-915afbf80120
  modified: 2026-10-14T01:50:31.639Z
---

User 2026-10-13 ~21:40: "you need to work on toys" and "1.5 h is too long for a test development loop; we need small
tests". Loop = Python on existing toySE outputs (~/toyse_snr03, ~/toyse_snr01: 20k, 4 states 13/14/15/17 at
10/40/30/20 %, refined poses, flex particle embedding run/flex/flex_pca_embedding.bin, pose-mean 2D sub-classes
view_tangent_labels_posemean_refined.txt, truth toy/gt_particles.txt). Scripts in the session scratchpad toys/:
admixKd.py (shared-cov GMM, subspace, admixture), avgs_only.py, kselect.py, viewsplit.py, gmm_seed.py, avg_seed.py,
ptcl_posterior_seed.py (repo ptcl_posterior.py + --seed-npz; POST_OUT env for the npz).

| method (accuracy, SNR .03 / .01) | |
|---|---|
| flex's delivered states (XD + GMM AUTO fold) | 0.71 / 0.41 |
| k-means / archetypes on sub-class means (averages only) | 0.49-0.40 / 0.46-0.43 |
| shared-covariance GMM, 12 starts, raw particle latent (8 comps) | 0.93 / 0.77 |
| + soft cross-half image E-steps (consensus + Wiener-shrunk deviation refs, fixed fractions, lp 15) | 0.97 / 0.88 (state 13 recall .94) |
| averages-only seed + same E-steps (SNR .01) | 0.67, state 13 never recovered |

- Sub-classes are mixtures (median dominant share 0.68 / 0.60); the E-step cannot create a state its seed lacks.
- Admixture (sub-class-weighted priors) adds little; between-sub-class subspace is NOT the state subspace.
- K: BIC over-splits (prefers K=7); parent co-occurrence fails. SOLVED on toySE by VIEW BALANCE after the E-steps:
  true states have TV(view distribution vs global) 0.02-0.08, fragments 0.24-0.39; rule = GMM at KMAX=6 -> 4 soft
  E-steps -> merge the state with TV > 0.15 into the partner whose union is most balanced -> repeat: ends at K=4 at
  both SNRs (0.971 at .03; 0.833 at .01). Script states_pipeline.py (2-3 min per run).
- Synthetic d'=0.8 latent: 8-D EM fails, 1-D fit along the averages' direction recovers separation/fraction
  (averages give direction, particles give scale); admixture MLE is misspecified at weak d' (theta_g overfit).

**How to apply:** develop placement/K/E-step in this loop first; Fortran targets = shared-cov GMM placement in flex
and a SOFT refine=states with consensus+deviation refs; validate on a ~10-min Fortran toy, not the 1.5 h chain
(toyexp.sbatch chains oc20_s03/oc30_s03 run in the background as end-to-end checks). Related:
[[flex-on-expanded-cavgs-plan-2026-10-13]], [[flex-placement-not-latent-loses-small-states]],
[[coupled-2d3d-realistic-regime-result-2026-10-10]], [[softem-toyoc-minority-recovery-2026-10-08]].

**Open/closed toy (toyOC, 13 % minority, 20 deg swing, refined poses):** SNR 0.1 (20k) and SNR 0.03 (80k,
toyexp/oc20_s03): shared-cov GMM 0.99-1.00, +E-steps 0.999, view-balance merge K=4 -> 2 at both lp 8/20. Flex's own
placement 0.53 (XD macro-clusters override npreimages: 37k/32k majority split, minority 28 % pure). BUT this toy is easy
per particle (minority d' 4.3 in the latent vs RNase P S1 ~0.8): harder toys oc6_s03 / oc4_s03 (toyexp2.sbatch, AMP env)
pending. On oc20_s03 the 2D expansion DID make pure minority sub-classes (89 % of minority particles in >= 80 % pure
sub-classes) but they were tiny (median 19 members) and the pop >= 20 cut before solve3D_cavgs dropped them.

**Fortran (snap reseed21, uncommitted ~/SIMPLE_master1013):** flex `state_placement=gmm` = tied-cov GMM from the placed
targets + 12 deterministic k-means++ restarts on <= 100k subsample, best loglik, final fit on all, no AUTO fold, no
respawn; npreimages honoured exactly (>= 2) for gmm. oc20_s03 placement-only (infile= embedding cache, seconds):
default 0.532, gmm K=2 0.992. Also `report_state_view_balance` after placement (log only): K=2 TV 0.006/0.038; K=3
fragment 0.648; K=4 fragments 0.525/0.647. Next Fortran: port score_states (flex-softem simple_commanders_prob.f90,
~270 lines; deps exist in master1013 except shsrch new_legacy) writing the state_weight_set, cross-half refs.

**Fortran soft-EM drift (2026-10-13 late, reseed22-24, toyexp oc4_s03/oc6_s03 = 4/6 deg swing, 13 %, SNR .03, 80k):**
truth-seeded loops (reconstruct3D m_estimator=flex -> score_states cross-half, lp 10) start at the image-evidence
ceiling (oc4 0.91, d' 2.1 in R2-R1; oc6 0.95) then collapse to 0.69-0.77 by it 3-4 in EVERY arm: gridding, PCG own
maps, PCG state_refs=deviation, EM fractions, fractions fixed (SIMPLE_SCORE_PIS). Decomposition of R2-R1 vs truth:
the state separation stays (844 -> 802) but the per-view midpoint SD grows 58 -> 158 -> 408 and within-view SD
390 -> 682 = Einstein-from-noise reinforcement of view-specific parts of A2-A1 by near-hard assignment (median top
posterior 0.99+). Chi-square scale c = meff/median(best R) is 1.4-1.6x overconfident already at it 1 (truth-fitted
logistic slope x0.64-0.73; per-particle c_i does not help: R_best spread = chi2). Deviation refs do NOT help (the
reinforced artefact is reproducible across halves, so FSC weighting keeps it). Flex particle latent at lp 20 / lp 10
does NOT carry the 4 deg state (held-out LDA d' 0.15 / 0.20 vs 2.1 in image space): the seed must come from elsewhere.
Next: SIMPLE_SCORE_TEMP (factor on c) and SIMPLE_SCORE_CALIB=lda (tied-cov mixture on c(R_s-R_1), fixed fractions).

**SELF-TERM BUG (found 2026-10-13 ~23:40):** with lp set, read_mask_filter_refvols loads the MERGED map into both
even and odd matching slots (refine3D convention under lp; without lp it blends eo below 20 A), so score_states'
"cross-half" (SIMPLE_SCORE_XHALF) never was: every E-step scored particles against maps holding their own image
(~5 nats pull to the current state for a 10k-particle minority at box 96 / lp 10). Fixed (reseed25) by the runtime flag
params%l_refs_eo_independent (set in score_states; skips the lp merge and the low-res eo blend). refine=states and any
lp-set refine3D multi-state classification have the same self-term. With true cross-half: oc20 (20 deg) gmm seed stays
0.9993 for 5 its; oc4 truth-ref ceiling is only 0.89 (baseline 0.87; the 0.91-0.92 before was self-term), oc6 0.93;
every oc4/oc6 loop still collapses by it 3-5 (best x0.5 + fixed fractions 0.892 0.891 0.824 0.705 0.635); deviation
refs shrink the difference to zero (one state by it 4). lda score-space fit underestimates the scale (0.26) -> collapse.

**Amplitude ladder (2026-10-14, reseed26, true cross-half, toyexp oc{4,6,10,13,20}_s03):** whole-map particle latent
(lp 20, 8 eigs) LDA d' 0.15 / 0.22 / 1.19 / 2.65 / 4.3. gmm placement finds the state at 13 deg (TV 0.076, 0.942) and
20 deg, view fragments at <= 10 deg. Seeded E-step: 13 deg 0.942 -> 0.990 -> 0.993 ... 0.992 at it 5 (maps: state 2
r 0.999 to the truth-label map, difference map r 0.974, amplitude 0.89); 10 deg truth-seeded EM drifts 0.985 -> 0.939
by it 5 (fraction -> 0.16), gmm-seeded 0.774 -> 0.661. Focused flex (pcg_mskfile = moving-region mask) raises latent
d' 0.20 -> 0.64 (4 deg), 0.22 -> 1.26 (6 deg) but no clustering resolves a 13 % mode at d' ~1.3 (view-centring helps
little). RNase P particle latent d' ~0.8 sits between the 6 and 10 deg toys. RNase P averages latents (c16vc/c100vc
arms, gmm placement K 2-4 from the caches, gmmplace/): only big-sub-class S3/S5 outliers (6 / 34 averages, S3 x3.9)
and view groups (TV 0.55-0.72); no S1 / S2 / S4 group. The user's own rnasep_estep / alternation jobs (coupled build)
use a Python E-step, not score_states.
