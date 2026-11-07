---
name: flex-embedding-precision-is-map-sampling-precision
description: "TRAP: flex_pca_embedding.bin's precision block is the MAP SAMPLING precision Q = A G^+ A (map_sampling_precision), so Cov(z_MAP|x) = Q^-1; flex's own noise calibration/deconvolution read it as A (D = Q - P) and mest_charts' z_ML un-shrinking reads it as G"
metadata:
  node_type: memory
  type: project
  originSessionId: c6b87431-7325-4b60-ad85-a72f7e8559d4
  modified: 2026-11-07T17:30:31.950Z
---

Found 2026-11-07 (spectrum-stop session) by reading merge_master2's source. In fit/simple_flex_pca_embed.f90 every latent solve
stores `map_sampling_precision(Gtil, prior)`, which is Q = (Gtil + P) Gtil^+ (Gtil + P), never A = Gtil + P. So:
- the noise of z_MAP is Q^-1 (scaled by flex's half-split a or a_comp);
- Var(z_MAP) = Var(R x) + E[Q^-1] balances, and generalized eigenvalues of these two give per-direction SNR without needing G;
- G cannot be recovered from Q and P alone: Q/P = mu + 2 + 1/mu has two roots mu and 1/mu, and typical nodes sit at mu of about 0.3-3.

Misreadings:
- run/simple_flex_pca_state_service.f90 passes it to calibrate_noise_scale and deconvolve_latent, whose noise_and_projection takes
  it as A (D = precision - prior). The scalar a absorbs part of the mismatch (a = 0.4-0.55 everywhere).
- mest_charts.py builds z_ML = inv(prec)(prec + Pi) z_MAP, i.e. it treats prec as G. Its "observed / deconvolved / measurement" budget
  therefore fails at low SNR (olenull 5/12/3).

**How to apply:** use Q^-1 times a_comp as the measurement covariance of z_MAP, and do not un-shrink with it. If flex's deconvolution
or noise scale ever matters for a decision, flag the A-vs-Q mismatch to the user as a possible flex fix (not done). Related:
[[spectrum-stop-rule-verdict-2026-11-07]], [[tree-charts-coordination-gap-verdict-2026-11-06]].
