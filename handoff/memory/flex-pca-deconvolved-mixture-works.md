---
name: flex-pca-deconvolved-mixture-works
description: "A measurement-error mixture (extreme deconvolution) on the neigs=5 IgG embedding is the first method to beat chance — but only with the per-particle noise inflated ~13x, and likelihood cannot pick that factor"
metadata: 
  node_type: memory
  type: project
  originSessionId: 98926916-3af3-4bb7-9587-824545bc5ef2
  modified: 2026-08-05T18:55:55.391Z
---

2026-08-05. `scripts/flex_pca_mixture_deconv.py`. Extreme deconvolution
(Bovy/Hogg/Roweis) on the neigs=5 IgG embedding:
`z_i | k ~ N(mu_k, Sigma_k + c*Lambda_i^-1)`, Lambda_i^-1 = stored per-particle
posterior covariance (NO sig2_eff factor — see
[[flex-pca-latent-is-997-percent-noise]]).

Scored on the only honest GT variable, frames 0-19 vs 20-99 (prevalence 0.200):

| K | model | enrichment | z |
|---|---|---|---|
| 2 | plain GMM | 1.01x | +1.45 |
| 2 | deconv c=1 | 1.02x | +1.60 |
| 2 | **deconv c=13.03** | **1.13x** | **+9.81** |
| 5 | plain GMM | 1.09x | +4.47 |
| 5 | **deconv c=13.03** | **1.18x** | **+11.38** |

**First method in the whole investigation to beat chance from the embedding side**
— every flex_pca kernel state under every ordering/truncation sat at 0.89-1.33x
with no significant z. Calibrated runs also give a significantly DEPLETED
complementary component (z = -5.4), as a real two-population split should.

**Two things the controls prove.** (1) The NOISE MODEL is the active ingredient,
not the mixture: plain GMM on identical coordinates gives z=+1.45, adding the
per-particle noise term gives +9.81. (2) The CALIBRATION is essential: c=1 (the
model's own uncalibrated Lambda_i^-1) gives z=+1.60, indistinguishable from plain.
Only inflating ~13x unlocks it.

**SOLVED (partly) 2026-08-05: the split-data replicate calibrates c without GT.**
flex_pca ALREADY solves each latent twice from disjoint checkerboard halves of the
particle's own Fourier samples (`cov_herm_inner(..., half)`), uses the pair only
for one correlation, and discards it. Added `write_zhalf_replicates` in
simple_flex_pca_columns.f90, gated on **SIMPLE_COV_ZHALF** — purely additive, no
numerics change when unset. Estimator: c = Var(z1-z2) / model prediction of the
same quantity (= diag(2*A_half^-1), A_half=(A_full+P)/2), which cancels the prior
shrinkage a raw Var(diff)/4 would leave.

`c_split = 6.54` GT-free vs `c_gt = 13.03` oracle. Mixture at K=2: plain GMM
z=+1.40, c=1 z=+1.37, **c_split z=+5.91**, oracle z=+9.81. So a GT-free
calibration clearly beats none but reaches ~half the oracle z and a third of its
|corr|. Per-component c_q = [13.6, 11.95, 1.05, 1.38, 1.29] is NOT uniformly
better — wins at K=5 (+7.38 vs +6.05), loses at K=2 (+3.11 vs +5.91).

**FULLY SOLVED 2026-08-05 by adding the two half-set BASIS replicate.**
`heldout=yes` fits two bases on disjoint halfsets; `probe_external_basis`
(SIMPLE_COV_PROBEBASIS/PROBEN=0/PROBEEIG, with infile= skipping the refit) embeds
all 100k particles through each; align frames from the eigenVOLUMES, never from
the latents being differenced. z_A and z_B share the image so estimation noise
cancels and Var(zA-zB) ISOLATES the basis term. No Fortran change beyond the
SIMPLE_COV_ZHALF export.

  c_split (estimation noise) = 6.54  (62 %)
  c_basis (basis error)      = 4.04  (38 %)
  c_full                     = 10.58     vs oracle c_gt = 13.03

Mixture: K=2 c_full z=**+9.74** vs oracle +9.81 (99 %); K=5 c_full z=**+11.94**
vs oracle **+11.40 — the GT-free value BEATS the oracle**. That is the signature
of c_gt being an OVER-estimate: it is within-frame variance, which counts genuine
conformational spread inside an MD frame as error. True optimum ~10.6; the
residual 19 % is ANOVA over-counting, not estimator deficiency. **No ground truth
is needed to calibrate the mixture.** Script:
`scripts/flex_pca_basis_replicate.py`.

**Biggest finding underneath: THE BASIS IS UNSTABLE.** flex_pca's own diagnostic
gives mean principal-angle cosine **0.383** between two 5-D bases fitted on
disjoint halves of the same 100k particles. That is why basis error is 38 % of
per-particle error, and it makes basis stability — not the E-step — the next
target.

**TRAP: `heldout=yes` yields a cache that CANNOT be deconvolved.** It leaves
`latent_second` ZERO for halfset-A particles, and calls
embed_latents_with_contrast twice so flex_pca_zhalf.bin is overwritten with one
halfset. Take precision/split-term/mixture from a NORMAL-path run; take only the
basis term from the held-out run. Inverting those zero blocks corrupts everything
silently.

**Why the split-data term alone undershoots:** both halves go through the
SAME basis volumes, so the replicate is blind to basis error. c_gt/c_split = 1.99
says **roughly half of flex_pca's per-particle latent error is BASIS error, not
estimation noise** — pointing at the basis, not the E-step. Note c_q tracks
split-half rho exactly: the model's error estimate is right for well-measured
components (rho 0.83-0.86 -> c_q 1.05-1.38) and wrong by ~13x for z1/z2
(rho 0.42-0.49), i.e. most confident where least entitled, and z1 drives state
placement. Scripts: `scripts/flex_pca_calibrate_c.py`,
`scripts/flex_pca_mixture_with_c.py`.

**Likelihood still CANNOT pick c** —
in-sample ll is -25.6337 at c=1, -25.6595 at c=13.03, -25.6662 plain, so
likelihood PREFERS the useless uncalibrated model. Viable route: re-estimate each
particle's latent from two disjoint halves of its own Fourier samples; their
difference measures per-particle error directly. c is a single SCALAR, a far
easier problem than the per-component ranking that
[[flex-pca-signal-ranking-is-blocked]] shows is unavailable.

**Honest scale:** 1.18x and |corr(resp,GT)|=0.08 is overwhelming at n=1e5 but
practically modest — it does not yet separate the populations usefully. What it
establishes is that the information IS recoverable by POOLING (the only route left
per the per-particle-SNR argument) and that recovery is gated on noise calibration,
not on the embedding or the basis.

Implementation notes: the noise is genuinely heteroscedastic (per-particle
posterior variance CV 0.23-0.50 per component, Frobenius spread 38 % of mean), so a
constant-noise shortcut is invalid. Use BLAS matmul, not einsum — an einsum version
did not finish in 600 s where the matmul version takes ~3 min for all six fits.
