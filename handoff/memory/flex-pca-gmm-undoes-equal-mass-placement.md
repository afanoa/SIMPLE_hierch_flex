---
name: flex-pca-gmm-undoes-equal-mass-placement
description: The default tied-covariance GMM refit silently destroys equal-occupancy state placement on continuum data; skipping it turned an 88%-in-one-state RNA run into clean sextiles
metadata: 
  node_type: memory
  type: project
  originSessionId: 126275fa-8460-444b-8ba7-1e0140968435
  modified: 2026-08-13T17:31:59.625Z
---

Found 2026-08-13 while diagnosing the afm RNA continuum failure. `state_axis>0` in
`build_covariance_state_weights` (simple_flex_pca_model.f90) DOES place equal-occupancy slice
targets — quantile placement already existed. But two later steps undid it:

1. The kernel was scored with the full-`nk` posterior quadratic form, so off-axis (noise)
   directions set a `chi2(nk)` bandwidth instead of the target spacing.
2. `l_gmm` (default ON) then ran the tied-covariance mixture, which **re-fits the means**. On a
   continuum with one dense mode every component slides into that mode, erasing the sextiles.

Measured on 932k RNA particles, `npreimages=6 state_axis=1 lp=12`:

| | run 30 (before) | run 31 (after) |
|---|---|---|
| occupancy | one state 88%, rest 0.9-25% | max 16.96% (ideal 16.67%), all six 151-156k |
| frac_in_support | 0.014-0.962 | 0.333 all six |
| bandwidth | 5.32 from chi2(ncomp) | 13.8-36.1 from path-spacing, auto-wide at sparse ends |
| unassigned | — | 1.29% |

**Fix (implemented):** the `axis>0` branch now sets `ppath=z(:,axis)`, `tpath=tcen(axis,:)` and
`l_relpath=.true.` so distance and bandwidth live on the axis the slices were cut on; and
`l_relpath` skips the GMM refit unless `SIMPLE_COV_GMM=1` forces it back for A/B.

**Why it matters beyond this dataset:** the GMM was added to fix Epanechnikov coverage gaps, but
overlapping path-spacing kernels cover 98.7% on their own — so on equal-mass placement the GMM is
pure downside. Also retires `min_neff` there (path spacing `2N/nstates` dominates it) — one knob
gone toward [[flex-pca-must-be-knob-free]]. NOT yet checked against discrete synthetic data
(Ribosembly/IgG use `axis=0` k-means, so they take the untouched branch — but verify before
trusting). Run 31 = 1401 s end to end. Volumes not yet scored; occupancy being right is necessary,
not sufficient. See [[afm-rna-flex-pca-needs-continuum-readout]].
