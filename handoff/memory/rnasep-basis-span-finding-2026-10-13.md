---
name: rnasep-basis-span-finding-2026-10-13
description: "RNase P steps 1-2 (2026-10-13): cls_expansion sub-class arm dipoles vs the flex_pca 3D span. Envelope support hides half; sphere basis predicts the smooth part; the real, split-half-reproducible residual (0.30 vs null 0.04) is the 47-17 A fine structure of the same displacements, not a global 3D mode: the 3D basis is resolution-limited (~30-40 A effective), 2D averages resolve to ~17 A. Next: finer-band basis, then M-step on sub-class averages"
metadata:
  node_type: memory
  type: project
  originSessionId: 5b04614c-23e0-4650-8600-03b464a7f576
  modified: 2026-10-13T23:37:44.436Z
---

Plan `doc/implementation_notes/coupled_2d3d_rnasep_incremental_plan_2026_10_13.md`, results
`coupled_2d3d_rnasep_step1_basis_span_2026_10_13.md`, manifest `coupled_2d3d_rnasep_manifest_2026_10_13.md`,
scripts `scripts/cls_expansion_proto/subclass_pool/{basis_span,basis_span_all,resid_rec,clshalf_compare}.py`,
batch2 `rnasep_coupled/{basis_span_all,clshalf}` (jobs 53829943, 53833307).

- Diagnostic: per sub-class direct |CTF|-weighted deviation (basis band 117-16.7 A) regressed on the projected
  consensus (view), the projected eigenvolumes (span) and the members' stored z (latent); span null = randomly
  rotated basis (the only valid one inside a view class). G10 (automask support): span = null, latent ~6 %,
  dipoles straddle the projected support. Sphere 692k basis (`full692k/1_flex_pca`, 16 comp): latent 17-19 %,
  corr 0.41-0.53; all-class pose-mean expansion: E_latent 0.34, eo_resid 0.38.
- The residual is real (independent split-half expansions matched 0.30 vs 0.04 null), not alignment, quality or
  background (all ~null), and NOT a global low-rank 3D mode (resid_rec K=1-6, in-sample even->odd transfer 0.00-0.04
  vs 0.44 for a synthetic mode). Figure `rnasep_coupled_figs/clshalf_pairs.png`: the 3D prediction is the smooth
  version of the 2D dipoles; the residual is their fine stripes. Per shell: span leaves 0.3 at 117-59 A, 0.6 at
  47-23 A, 0.9 at 21-17 A; eigenvolume power falls 50x from 59 to 20 A.
- Reading: the 2D side's gift is RESOLUTION of the motion, not a new direction; the tied-mixture placement cannot
  use it. Order now: finer-band sphere basis (lp 10-12, neigs 16-24) -> M-step on sub-class averages (needs poses
  to 2-3 deg at 20 A) -> guided split -> tied-mixture A/B last.
- Local model built (night): sub-class tracking across the sphere (`subclass_track.py`: Hungarian matchings between
  near-view classes, in-plane angle = minus the e3 difference, labels propagated from a seed): cycle consistency
  0.72-0.80 (chance 0.07-0.1), seed-robust, split-half 0.32 vs 0.21, reaches 48/87 classes. The global latent does
  NOT index the fine structure even within 12 deg (local decoder 0.05 of the residual). Group maps (one binary
  selection + reconstruct3D nstates=1 per group; import_particles with stktab resets states, oritab route dead):
  tracked groups carry ~2x the reproducible between-group map variance of latent k-means groups at 20-8 A; all maps
  view-limited (cFAR 0.001) until the tracking reaches the whole sphere.
- Formulation (doc/implementation_notes/coupled_2d3d_math_formulation_2026_10_13.md): the two latents couple in the
  COARSE subspace (shared/private factor model per view); the fine structure couples only through the state volumes:
  block-permuted regression (unlabeled sensing, r-local permutations) + permutation synchronisation + heterogeneous
  MRA; solved by alternation (pose-refine each sub-class average against ITS state map, deviation-score, Hungarian
  per class with the latent prior, reconstruct from all particles with composed poses). Overnight driver job
  53843745 (`rnasep_coupled/alternation`, 3 iterations + no-refine control + import gate); whole-sphere tracked
  maps = V^0 (`track_rec2`: tracked 2x the latent arm's reproducible between-state variance at 117-30 A).
  Latent per particle cannot carry states (LDA 0.203 vs 0.186 chance). import_particles+stktab resets states;
  poses are kept (gated).
- Fast tests (night): synthetic truth with real geometry/noise: assignment + pose refinement are sound (labels
  0.70 -> 1.00 in one iteration, poses 4.3 -> 1.1 deg with TRUE maps); the Fourier-grid V-step from averages is
  the weak link (0.7 vs truth even with 87 views) -> V-step on the cluster from all particles. Proxies that
  evaluate with the same labels on both halves are circular; the CROSS-HALF metric (assign from even, reconstruct
  from odd, and vice versa; `assign_cross_metric.py`) is the honest fast metric: best score = band 30-12 A + beta
  0.1 (0.64/0.63/0.51 vs flat 100-12 beta 0: 0.60/0.55/0.44). Cluster it1: refinement helps (main > norefine),
  flat-score reassignment loses the fine bands; v2 driver (score 30-12, gain 0.01) launched.
- After midnight: SOFT assignment (Sinkhorn responsibilities, tau 0.05) is the largest single gain: cross-half
  between-state variance 0.52/0.55/0.45 (hard) -> 0.86/0.85/0.69 (soft); realised on the cluster by splitting each
  sub-class's particles among states in proportion (driver v5). Negatives: plain spectral permutation
  synchronisation < greedy tracking; fused GW on the latents matches neighbouring classes at chance (2D latent
  geometry) / coarse only (3D latent): fine correspondence needs the images. Per-particle E-step vs current maps
  weak (0.18 vs 0.11 chance). Drivers v1-v5 + estep + lp10 ran overnight 2026-10-14 (`rnasep_coupled/alternation*`).
- Held-out-class prediction of 3D modes is weak even with all 112 views (synthetic 0.05-0.07): judge by the
  in-sample even->odd transfer.

**Traps found:** stacks are phase-flipped (`ctf=flip`): weight by |CTF|; cls_expansion_cavgs.mrc renders are
written origin-at-corner (fftshift; block k = class-map ids k*100+1..); stored z and the pc files differ by a
rotation of the near-degenerate subspace (fit a linear map on ~1500 particles); this tree's `pool_place.py` reader
was 8 bytes off for cache version >= 4 (fixed); projection convention at box 80: conv 0, shift sign -1, half-pixel
scale; cls_expansion labels of a state-selected 932k project index the 932k rows; batch2 python 3.12 has numpy 2.1
(3.6 default has 1.14); 10 % of pose-mean sub-classes have members spread > 60 deg (filter by median angle).
Related: [[rnasep-2d-subclasses-pooled-vs-particle-3d-latent-2026-10-13]], [[rnasep-s4-four-arm-support-test-2026-10-08]],
[[cls-expansion-vs-flex-pca-latent-diagnosis-2026-10-09]].

**lp 10 sphere basis (2026-10-14, job 53834775, 24 comps, 4 h basis + accumulation):** pooled-latent residual over
deviation per shell 0.67 / 0.39 / 0.39 / 0.71 / 0.96 at 117-59 / 47-34 / 29-23 / 21-17 / 16-12 A (16 A basis: 0.69 /
0.50 / 0.75 / 0.99 / 1.00); the remaining residual reproduces at 0.53-0.70. The finer band extends the 3D model of the
arm to ~20 A; nothing below 16 A. Absolute-deviation span test (`span_test.py`, rotated-basis null): 16 A sphere basis
explains 0.55 of the 117-17 A deviation energy (null 0.30), per shell to 26 A, zero beyond 21 A.
