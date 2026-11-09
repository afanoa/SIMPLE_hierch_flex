---
name: sqrt-weight-on-planes-tempers-flex-estep-2026-11-09
description: TRAP: per-particle weights applied as sqrt(w) on the flex plane cache temper the E-step (|z| scales with w at flex's SNR), collapse low-weight particles at the latent origin and break the k-center state placement; Bishop-Tipping weights belong in the M-step sums only
metadata:
  type: project
---

Branch nodelevers 7d51f644d (SIMPLE_FLEX_WEIGHTS) puts sqrt(w) on planes / transfer at `plane_store_whiten_bands`; the cache
serves both E- and M-step. Measured 2026-11-09 (node_levers/handoff SOFT_LOOP.md): median |z| .8 at w < .2 vs 24 at w >= .9 in
the r13 -> r131 fit (corr(w,|z|) .69; qtoyflex .73); the ~30k soft extras form a spike at the origin, the GMM merges 5 of 8
k-center targets into one 65 % state, 4 states fall under the 2,000-n_eff floor, 948 hard members end unlabelled (state 0,
excluded from stA0 only). The soft arm's r131 split (x1.37-1.42 vs hard x1.00) follows this four-state seeding, not a better
latent (AUROC .716 vs .713).

Control (members at 1, extras at their weights): all 8 states kept, still x1.27-1.34, extras still collapsed (|z| .8 vs 19.8):
the state count is not the cause; the tempered extras in the fit are.

**Why:** at 15 A / box 64 a single particle's latent information A is far below the unit prior, so the posterior mean
(wA + I)^-1 w b scales with w; Bishop-Tipping (eqs. 27, 31-36) weights only the M-step log-likelihood terms.

**How to apply:** for a soft (weighted) flex, put w_i inside the M-step accumulations (basis / mean normal equations, noise
moments, state n_eff sums) and leave the plane cache unweighted; never read a weighted fit's state partition as evidence of soft
membership. Related: [[flex-pca-must-be-knob-free]], [[soft-carry-no-gain-on-fixtures-2026-11-09]].
