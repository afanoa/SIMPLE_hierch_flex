---
name: flex-pca-states-stratify-by-neff
description: "Under the default GMM refit, delivered flex_pca states differ mainly in effective sample size, not conformation — noise level tracks neff at r=-0.99, and SIMPLE_COV_GMM=0 fixes it"
metadata:
  type: project
---

Measured on EMPIAR-10049 (RAG SEC, 108,544 ptcls, 2026-08-17). With the default tied-covariance
GMM refit, the 16 delivered states had `neff` spanning 20,698-89,731 and solvent noise spanning
4.8x, with **corr(log solvent_sd, log neff) = -0.987 (16 states) / -0.998 (6 states)** and
`sd*sqrt(neff)` flat to 1.64x. The states were stratified by effective sample size, NOT structure.

The inversion is the giveaway: state 11 had **416 particles but neff=88,096** (`frac_in_support`
0.984 - a kernel broad enough to span everything) while state 15 had **15,804 particles but
neff=17,513**. A 38x smaller state carried 5x the effective count, so its map is smooth and looks
like the consensus.

`SIMPLE_COV_GMM=0` (equal-mass placement, kernel instead of the GMM refit) fixes it: neff
4,418-5,270 (1.2x), solvent noise flat, corr weakened to -0.473. Same lever as
[[flex-pca-gmm-undoes-equal-mass-placement]].

**Why:** any per-state metric — deviation energy, spread, pairwise correlation — is confounded by
this. A "deviation energy" ranking computed under the default GMM is just ranking 1/sqrt(neff);
the three highest-deviation states were the three LOWEST-neff ones.

**How to apply:** print `neff` and solvent-shell sd per state and check
`corr(log sd, log neff)` BEFORE interpreting any state comparison. If it is near -1, the states are
noise-stratified and the comparison is meaningless — rerun with `SIMPLE_COV_GMM=0`.
See [[flex-pca-state-reproducibility-gap-test]].
