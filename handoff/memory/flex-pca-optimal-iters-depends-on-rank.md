---
name: flex-pca-optimal-iters-depends-on-rank
description: TRPM4 neigs=20 collapsed 10->2 reproducible components over 4 EM iters while neigs=10 kept all 10; the optimal probe-iteration count falls as rank rises
metadata: 
  node_type: memory
  type: project
  originSessionId: 7e382bce-970c-4cdb-bc32-4a1a1c1801a9
  modified: 2026-09-16T15:25:59.470Z
---

On TRPM4 (162,326 ptcls, box_crop=64, lp=14, PCG + solvent mask, 2026-09-16)
`neigs=10` delivered ALL 10 components >=0.85 (median 0.918) at `n_probe_iters=4`,
but `neigs=20` at the same 4 iterations delivered only **2 of 20**:

| iter | S(t) | >=0.85 of 20 | engine's own n>=0.9 |
|------|--------|--------------|---------------------|
| 1    | 0.8166 | 10           | -                   |
| 2    | 0.7209 | 3            | 13/12               |
| 3    | 0.6050 | 3            | 7/6                 |
| 4    | 0.5455 | **2**        | 4/4                 |

Iteration 1 was the BEST basis at rank 20 (10 of 20 >=0.85, 13 of 20 >=0.84).
So rank 20 is not simply "too high for the data" — the representation carries
>10 axes initially and the EM iterations destroy them, faster at rank 20 than
at rank 10 because there are more free parameters per particle to overfit.

**Why:** this refines [[flex-pca-em-iterations-overfit-cross-half]], which recorded
"4 iters optimal" — that calibration is rank-dependent, not universal. As rank
rises the optimal iteration count FALLS.

**How to apply:** do not reuse `n_probe_iters=4` when raising `neigs`. The
untested configuration that would actually settle TRPM4's true dimensionality is
`neigs=20 n_probe_iters=1` or `2`. Never read flat eigenvalues alone as evidence
of a high-dimensional continuum: TRPM4's rank-20 spectrum was flat (204.97 ->
114.09, only 1.80x) AND non-monotone (6>5, 10>9, 13>12, 15>14, 18>17) — permuting
order is the noise-floor signature. The deconvolution agreed independently:
signal fraction (deconvolved/observed) was 1.035-1.081 for z8-z20 vs 1.344 for z1,
i.e. no measurement-error inflation left to remove on the tail axes.
See also [[flex-pca-basis-rank-over-specified]].

State placement survives the bad basis and so is NOT a check on it: all 5 states
came out hard one-hot (`neff == population`, sum = 162,326 = N) with bands
7.96-8.86 A, but mass spread degraded to 3.55x (51,894 vs 14,598) against rank
10's 19.2%, and the 9% minority state's res(0.5) was 15.60 A vs 8.9-10.5 A.
