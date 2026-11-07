---
name: not-averages-flex-heterogeneity-evidence-2026-10-14
description: NOT expanded averages DO carry heterogeneity flex resolves -- shown only with local view centring (kNN 25) + sub-class-split halves vs a NOT-matched homogeneous null (r 0.87/0.79 vs 0.11/0.02 at 30-15/15-10 A); even/odd gates and shuffle controls are route artefacts (latent = f(view))
metadata:
  node_type: memory
  type: project
  originSessionId: 34b66a02-03d2-42df-91f8-915afbf80120
  modified: 2026-10-15T04:32:25.827Z
---

2026-10-14/15 (autonomous). NOT route:
- kernel100_expansion (126 parents x 100 sub-classes = 25,200 even+odd averages)
- particle-mode solve3D on the averages from the human 80S, 7.8 A (not/s3d420_init80S)
- flex on the averages: flex_cavgs=model neigs 10 lp 12 nkern 2 equal_occ K5 view-centred (not/flexavg_80S)

**Artefact mechanism (decisive).** The latent of averages is a function of viewing direction, even with true poses.
- States along the latent are view groups.
- Their reconstructions differ in coverage and anisotropy, which reproduces between ANY halves sampling the same views.
- So these all fire on homogeneous data:
  - even/odd latent agreement (toyh1 0.72-0.87);
  - the end-state difference;
  - parent / block shuffle margins;
  - view-balanced end states;
  - sub-class-split halves.
- Homogeneous true-pose sim, uniform views, sub-class split: r 0.70/0.69/0.66.

**Valid test.**
1. Labels from LOCALLY VIEW-CENTRED latents: each row minus the mean latent of its 25 nearest rows in projection
   direction, excluding its own sub-class; equal-occupancy bins along PC1 (knnstates.py).
2. Halves split BY SUB-CLASS (eo = sub-class parity; subhalf.sbatch).
3. A homogeneous null MATCHED to the dataset: averages simulated from its consensus at its own poses, counts and SNR
   (SNR from the sub-class half repro: 0.61 -> 1.6), then the same pipeline (notnull.sbatch).

| (r at >30 / 30-15 / 15-10 A) | flex labels | kNN-centred |
|---|---|---|
| NOT-matched null | 0.86 / 0.27 / 0.18 | 0.63 / 0.11 / 0.02 |
| NOT | 0.96 / 0.88 / 0.80 | 0.93 / 0.87 / 0.79 (5.2 / 3.7 / 2.9x noise) |

- NOT carries real heterogeneity at 30-10 A.
- Caveats: the null has exact poses and Gaussian noise. Real pose errors explained at most 0.26 of the latent in toyh1.
- Maps: not/flexavg_80S/subhalf/knn25/knn; local ~/Downloads/NOT_states_localviewcentred_box210.
- kNN centring is weak where views are sparse: with uniform views, 25 neighbours span about 9 deg and the fine view
  dependence survives.
- Spliceosome, re-tested with the sub-class split (2026-10-15, avgmode_eval.sbatch) looked positive -- but see the
  strictest test below: NOT detectable once halves share no particles.
  - Real states r 0.42/0.38/0.37 (34-class full box) and 0.58/0.39/0.34 (200-class box 160).
  - Matched nulls r ~0.0 in every band (0.16 at >30 A for 200-class).
  - Biological check, frame-free: within parents, a sub-class's averages-mode path position tracks the B-run
    particle-flex state make-up of its members (34-class: B4 r -0.26, z -8.9; B2 +0.17, z +5.0; B1 +0.14, z +4.5;
    200-class weaker, B3/B4, z -4.7/+4.3).
  - The earlier "spliceosome no" was the flawed even/odd metric.

**Strictest test (2026-10-15): fully independent halves (crosstest.sbatch) vs the matched null under the same split.**
- Half A = even-particle averages of odd sub-classes; half B = odd-particle averages of even sub-classes.
- Needed because kernel siblings share particles: the sub-class split alone leaks shared noise.
- Results, r at >30 / 30-15 / 15-10 A, real vs null:
  - NOT: 0.79 / 0.75 / 0.67 vs 0.2 / 0.04 / 0.0 -> strong.
  - RNase P from the S2 starting map: 0.72 / 0.58 / 0.35 vs 0.5 / 0.38 / 0.1 -> passes, but ab initio poses
    (solve3D_cavgs) give 0.42 / 0.37 / 0.14 vs null 0.43 / 0.39 / 0.09 -> NOT detectable. The S2-seeded pass is
    possibly reference imprinting.
  - RNase P biological check (known S1/S2/S4 selections of the 932k project, row-matched):
    - the 10-component averages latent predicts CV R2 0.01-0.03 of the sub-class S shares;
    - 30 components (lp 12): R2 0.13 / 0.08 / 0.07; lp 8 is worse;
    - 30 components placed on nkern 5: the path tracks S1 vs S2/S4 (S1 r -0.26, z -28).
    - So the state signal lives beyond 10 components. Ab initio 30-component latent: CV R2 0.26 / 0.15 / 0.13,
      predictive direction even/odd 0.8.
    - Strict test of 30 components placed on 5 (rnp30.sbatch): ab initio real at noise vs null 0.5-0.8 / 0.39 / 0.1;
      from S2 real 0.60 / 0.37 / 0.23 vs null 0.48 / 0.34 / 0.14 -> NOT convincing.
    - The latent carries the states; unsupervised placement does not turn them into reproducible maps.
- Overnight verdict: only NOT passes the strict test.
- **Knob-free placement (2026-10-15 05:20):** choose placement components whose even/odd agreement EXCEEDS the
  matched null's (excess > 0.15); the components that agree as much in the null are view structure.
  - RNase P ab initio 30 comps: comps 1-13 view, 14-30 excess -> path tracks the known states S1 r -0.28 (z -29) --
    best -- but maps still at the null (blend limit).
  - NOT: comps 4, 5, 6, 8, 9, 10 -> strict test real 0.89 / 0.79 / 0.66 vs null 0.2 / 0.05 / 0.01 (stronger than
    nkern 2).
  - Maps ~/Downloads/NOT_states_auto_components_box210. Not in Fortran yet.
  - Spliceosome (34 and 200 classes): noise -> NOT detectable. The sub-class-split "yes" was kernel-sibling noise.
  - toyc1: correct yes. toyh1: correct no; its null is not at noise (SNR 5.3 view artefact), so always compare
    against the dataset's own null.

**In flex (reseed33, local branch ~/SIMPLE_master1013, NOT pushed): SIMPLE_FLEX_VIEW_KNN=k.**
- simple_flex_pca_application centre_latent_by_view_knn subtracts the mean latent of the k nearest rows in projection
  direction (antipodally folded, self excluded).
- It is applied AFTER the latent deconvolution, so resumes are centred too.
- TRAP: a resume (infile=) adopts the cached deconvolved block, so centring raw z before the deconvolution does nothing
  on a resume (reseed32 was invalid).

NOT vs the matched null, nkern 2, equal_occ K5, sub-class-split halves (r at >30 / 30-15 / 15-10 A):
- NOT: 0.87 / 0.85 / 0.79 (3.8 / 3.5 / 2.9x)
- null: 0.38 / 0.10 / 0.00 (1.5 / 1.1 / 1.0x)
- With all components the null keeps r 0.53 at >30 A, so use nkern 2.
- Kernel-weighted maps: not/flexavg_80S/averages_mode_maps; local ~/Downloads/NOT_states_averages_mode_box210.

Related: [[averages-basic-sim-placement-failure-2026-10-14]], [[view-imbalance-inflates-between-state-metric]],
[[averages-route-real-data-2026-10-14]].
