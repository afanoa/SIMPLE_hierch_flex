---
name: flex-pca-signal-ranking-is-blocked
description: "Ranking flex_pca components by signal fraction is worth ~2x but is BLOCKED — both GT-free estimators in the code fail, and split-half rho (r=-0.77 vs ICC) is currently steering state placement"
metadata: 
  node_type: memory
  type: project
  originSessionId: 75235947-f21e-4c4b-ab6c-e947c2f45046
  modified: 2026-08-05T17:17:36.444Z
---

Measured 2026-08-05 on `~/igg_flexpca_probe` (100k particles, 48 components),
per-component signal fraction from GT-frame ANOVA as truth.

**The gain is real.** In the whitened metric — which is what matters, since
[[flex-pca-probe-iters-and-whitening]] says whiten before any distance and
`wcomp` already standardizes to 1/var — mean ICC over the selected subspace:

| select top-k by | k=5 | k=10 | k=20 |
|---|---|---|---|
| eigenvalue (current) | 0.542 % | 0.348 % | 0.256 % |
| split-half rho | 0.146 % | 0.380 % | 0.268 % |
| TRUE ICC (ceiling) | 1.046 % | 0.753 % | 0.532 % |

So ~2x at k=5, confirming the earlier estimate — **but only with an oracle**.

**Both GT-free estimators in the code FAIL.** Spearman vs true ICC over 48
components: eigenvalue +0.056, `component_reliability_proxy`
(spread/posterior-variance, used on the resume path) **-0.009**, and — over the
10 components whose values are logged — split-half `rho` **-0.770**.

**Why split-half rho fails is structural, not a bug.** `cov_herm_inner(half=)`
splits each particle's Fourier samples on a checkerboard parity
(`modulo(h/pf + k/pf, 2)`), which is a sound independent-noise split. But
split-half reliability answers "is this component reproducible from THIS
particle's own data", and any per-particle NUISANCE (residual pose error,
contrast/scale, CTF misfit) is perfectly reproducible and perfectly
non-conformational. Reliability != signal fraction. Conformational signal needs
variation BETWEEN conformations relative to spread WITHIN one, which a
same-particle split cannot see.

**This corrects fix (2) in [[flex-pca-latent-is-997-percent-noise]]**: "split each
particle's Fourier samples in two and embed both, the difference being a direct
noise replicate that gives a GT-free ICC" — that is already implemented
(`zhalf`), and it does NOT give a GT-free ICC, for the reason above.

**Live defect:** `comp_rho` is not diagnostic-only. It is passed to
`diffusion_kcenter_targets` and `reliability_path_targets`
(`simple_flex_pca_model.f90:768,777,787`), so state-target placement is currently
steered by a quantity anti-correlated with conformational signal. Placement is
otherwise standardized (`wcomp` = 1/var, `SIMPLE_COV_STDZ`), so this is the main
route by which the bad ranking reaches results.

**Caveat on the -0.770:** `rho` is logged only for `min(ncomp,10)` components
(`simple_flex_pca_columns.f90:2988`), and those 10 are the leading ones by stored
order, so it is a biased subset and n=10. Raising that bound — or better, adding
`rho` to the embedding cache — is a one-line change and is the cheapest way to
confirm on all 48.

**Conclusion: do NOT reorder by either available proxy** — it is worse than
eigenvalue order. The reordering is trivial once a valid estimator exists; the
estimator is the blocker, and no same-particle statistic can be one.
