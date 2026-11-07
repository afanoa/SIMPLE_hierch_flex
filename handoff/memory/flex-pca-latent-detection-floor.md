---
name: flex-pca-latent-detection-floor
description: "Latent-space structure tests on IgG are blind by 35-106x: they need 10-30% of observed variance, conformational signal is 0.28% — so clustering/multimodality results on these embeddings are uninformative either way"
metadata: 
  node_type: memory
  type: project
  originSessionId: 75235947-f21e-4c4b-ab6c-e947c2f45046
  modified: 2026-08-05T16:44:56.386Z
---

Measured 2026-08-05 with `scripts/flex_latent_structure_test.py` (GT-frame ANOVA
for the noise scale, planted controls on an amplitude LADDER, statistics run in
both the eigenvalue subspace and the between-frame signal subspace). Results in
`reports/flex_structure/`.

**The detection floor**, as a fraction of observed latent variance, against a
covariance-matched null at 3 sd:

| statistic | threshold | vs measured signal |
|---|---|---|
| arch, quadratic-in-PC1 | f >= 0.1 | blind by 35x (probe), 8x (c64_p0) |
| arch, principal-curve dVE | f >= 0.1, often NONE | unusable, see below |
| mixture, deconvolved BIC | f >= 0.3, or never | blind by 106x |

Measured conformational signal: **0.282 %** (probe, reproduces
[[flex-pca-latent-is-997-percent-noise]] to 3 s.f. via an independent code path),
1.32 % (c64_p0).

**Consequence: multimodality on these latents cannot be answered, and neither can
absence of curvature.** Every frame/run combination returns CANNOT RULE OUT. This
also independently explains why flex_pca's hard state labels carry zero GT
enrichment — state assignment is running ~100x below its own detection floor, so
it cannot be doing anything but clustering noise.

**Why binning cannot rescue it:** bin-averaging 1000 particles cuts noise 1000x,
which naively should expose 0.28 % signal. It does not, because the bins are
ordered by OBSERVED PC1, which is 99.7 % noise — so bin order is nearly random
with respect to the true conformational coordinate. Any statistic that sorts on a
noisy latent axis inherits this. Sort on something independently reliable instead.

**Retire the principal-curve dVE statistic.** It has no power here and is
sometimes INVERTED: on c64_p0's signal frame, planting an arch at f=0.3 gives
dVE 0.5756 against a null of 0.5992 — more planted structure, lower score. Any
past conclusion resting on dVE (including the original "no arch" call in
[[flex-pca-latent-has-no-arch]]) is weaker than it looked.

**RETRACTED: the "~70x non-conformational structure" claim.** It came from
reading the real quadratic statistic against that run's null, and THAT NULL WAS
BROKEN — see below. Corrected picture:

- **Per-particle metadata explains essentially nothing.** Regressing the top PCs
  on the ptcl3D fields (`scripts/flex_pc_metadata_regression.py`, results in
  `reports/flex_pc_metadata/`) gives R^2 for PC1 of 0.0001 defocus, 0.0000
  orientation, 0.0001 shift, 0.0000 halfset, 0.0000 quality — **0.0003 for all
  blocks jointly**. Top univariate |r| is mean defocus at 0.020. Residualizing
  all metadata removes only 16.5 % of the quadratic statistic. So the latent is
  NOT driven by defocus, pose, shift, halfset, or any recorded quality field.
- **A genuine excess does survive a correct null, and it is NOT conformational.**
  Against a Gaussian copula null (exact real marginals, Gaussian dependence) the
  quadratic statistic is z = +4.6 in the eigenvalue subspace and z = +23 in the
  between-frame subspace; the corrected script's own matched null independently
  gives +3.0 and +40 on the same data, so two null constructions agree.
  Reproduces out-of-sample (select subspace on half A, measure on held-out B:
  z = +25.2).

  **But a permuted-label control matches it exactly: z = +26.0.** Choosing the
  "signal" subspace from SHUFFLED frame labels — zero conformational content —
  produces the same effect as the genuine one. So the nonlinearity is a generic
  property of LOW-VARIANCE latent directions, not of conformation, and the
  eigen-vs-signal frame difference is about variance rank, not signal. Marginals
  are near-Gaussian (excess kurtosis 0.09-0.46, |skew| < 0.07), so heavy tails
  are not it either. Origin is the estimator's output geometry; nothing found in
  this whole arc supports a nonlinear CONFORMATIONAL manifold.

  Method note: always run the permuted-label version of any GT-selected subspace.
  The split-half control alone did NOT catch this — held-out reproduction looked
  like validation, and only label permutation exposed it.

**THE NULL BUG (important, cost a wrong conclusion).** `arm_null` built synthetic
arms as `fill(per-component leftover variance) + per-particle noise`. The noise
term's mean covariance is NOT diagonal in the PC frame, but the real data's PC
covariance is diagonal BY CONSTRUCTION. So every synthetic arm carried
component-to-component correlation up to |r| = 0.21 that the real data cannot
have, plus 37 % excess variance on PC2. A correlated null makes the bin-mean
profile nearly deterministic, collapsing the quadratic-increment baseline to
0.0001 — about 100x below the ~d/n ~ 0.01 expected from fitting 3 parameters to
100 bins — which manufactured z = +826.

Rule this generalises to: **when the real data is expressed in its own PC frame,
any synthetic comparison arm must match the FULL covariance matrix, not
per-component variances.** Fixed by `match_noise_mean` (symmetric whitening so
mean_i Psi_i equals the ANOVA within-frame covariance exactly) plus a bisection
for the largest noise scale keeping `target - struct - beta*noise` PSD. Arms now
match to <1 % variance and <=0.01 off-diagonal.

Planted-vs-null THRESHOLDS are less affected, since both arms shared the flawed
construction — but real-vs-null z-scores from the original run are void.

**Method trap worth keeping:** a threshold read as "smallest f that beats the
null" latched onto an isolated low-f fluctuation and produced a spurious
STRUCTURE DETECTED verdict on c64_p0. Require MONOTONICITY — smallest f such that
it and every larger f clear the null. Fixed in the script.
