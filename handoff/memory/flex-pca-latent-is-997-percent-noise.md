---
name: flex-pca-latent-is-997-percent-noise
description: "Measured on 2026-08-05 — 99.7 % of flex_pca's IgG latent variance is per-particle estimation noise, the stored precision understates it by ~1e6, and eigenvalue order is uncorrelated with signal"
metadata: 
  node_type: memory
  type: project
  originSessionId: 98926916-3af3-4bb7-9587-824545bc5ef2
  modified: 2026-08-05T16:52:25.837Z
---

ANOVA over the 100 IgG GT MD frames (1000 particles each share a conformation, so
their spread IS estimation noise), run `~/igg_flexpca_probe`, 100k particles,
48 components:

- total latent variance 2.243e4; **within-frame 99.72 %**; between-frame
  (all conformational signal, bias-corrected) **0.28 %**.
- per-particle conformational SNR = **0.0028**; ~**353 particles** must be averaged
  for unit SNR. That part is physics at SNR 0.01, not a defect.

**There is NO sig2_eff scale bug — an earlier version of this memory claimed one
and was wrong.** The stored matrix is a posterior PRECISION `A`, and the main path
(`embed_latents_with_contrast`, simple_flex_pca_columns.f90:2759) builds
`A = (a^2/sig2)*G + diag(1/Gamma_q)` — already the posterior precision, so
`Cov[z|y] = A^-1` with **no** sig2_eff multiplier. Applying one (taken from a
comment on the differently-scaled `prepare_projected_latent_estep_part` path)
manufactures a spurious factor of 5e5. Correct numbers: model posterior variance
`sum diag(A^-1)` = 2901 = **12.93 %** of observed latent variance, vs measured
within-frame noise 99.72 % — an understatement of **7.71x**, not 1e6x. And 7.7x is
not obviously a defect: A^-1 is the error CONDITIONAL on basis and pose being
right, while within-frame spread also absorbs basis error, pose error and CTF
misfit. **Do not "fix" this scale — editing it would introduce a bug.**

Consequence: the measurement-error mixture model
`z_i|k ~ N(mu_k, Sigma_k + Lambda_i^-1)` is **available today, not blocked** —
though Lambda_i^-1 should be inflated by a calibration factor (~8x here).

**Consequence: eigenvalue order ranks noise.** Spearman(eigenvalue, ICC) = **+0.056**
over 48 components. z1 has the LARGEST eigenvalue (31 % of the mass), an ICC of
0.0030, ranks **17th of 48** by signal, and has AUC 0.489 (below chance) for the
real contrast — and it is the axis the trajectory states are placed along. Top-5 by
eigenvalue captures 21.3 % of signal; top-5 by ICC captures 41.2 %.

**Downstream:** every state's GT enrichment is 0.96-1.01x against prevalence 0.200,
i.e. the hard labels carry ZERO conformational information; kernel neff is
22.6k-38.1k (23-38 % of the dataset) when ~353 particles suffice, so ~80x
over-averaged; 75 % of particles are unassigned (state 0). LDA over all 48 whitened
components, fit-even/score-odd, gives AUC 0.5275.

This also explains why flex_pca's state maps differ REPRODUCIBLY (ceiling +0.948)
while being conformationally identical — that difference is entirely
non-conformational, which is why a conformation-specific ENM basis scores most
negative on it. See [[flex-deformation-model-falsified-on-igg]].

**Ranking components by signal fraction was TESTED on 2026-08-05 and does NOT fix
the state stage.** Permuting/truncating the cache and replaying with `infile=`
(which re-runs only the state stages) across four configurations gave GT enrichment
0.89x-1.33x, none statistically significant. Ranking raises the fraction of
available signal captured (41.2 % vs 21.3 % for top-5) but the available signal is
0.28 % of variance, so the best-5 subspace is still ~98 % noise per particle.

**The binding constraint is per-particle SNR, not component order.** Honest
per-particle discrimination is AUC 0.5275 (near chance) while GROUP means separate
strongly — the embedding has population-level information and ~zero per-particle
information. Kernel regression weights each particle by its own z, so it is
inherently per-particle and cannot use population information. That is an
architectural mismatch, not a tuning problem. Unit SNR needs ~353 particles, so the
right operating point pools GROUPS of ~200-400; current kernels pool neff
20k-98k (50-250x too many) and hard assignment pools 1.

Two mechanical gotchas found: `nkern` does NOT narrow the kernels (the log says
"bandwidth floor from chi2(ncomp)" — computed from ncomp, not nkern), and k-means
collapses on the featureless noise blob, putting 97-99 % of particles in one state
under every ordering.

**No GT-free ranking exists today.** `component_reliability_proxy`
(simple_flex_pca_model.f90:2268) computes the structurally right quantity
rho=ratio/(1+ratio) with ratio=spread/posterior_variance, but measured against true
ICC it scores Spearman **-0.013** — worse than plain eigenvalue order (+0.056) —
and its top-5 captures 8.6 % of signal vs 21.3 % for eigenvalue. It inverts
`precision` CORRECTLY (an earlier claim that it omits a sig2_eff factor was wrong).
It fails on DYNAMIC RANGE: its rho spans only 0.862-0.946 across all 48 components,
so there is nearly nothing to rank on, and it reports ~0.89 where true ICC is
~0.003 because A^-1 understates real per-particle error ~8x.

**REFIT AT neigs=5 (2026-08-05): per-particle SNR improves 6.27x, discrimination
does not.** Everything else matched (ncols=256, n_probe_iters=3, box_crop=64,
lp=15). SNR 0.00284 -> 0.01783; particles for unit SNR **352 -> 56**; best ICC
0.0172 -> 0.0297; mean ICC 0.0026 -> 0.0117; noise share 99.72 % -> 98.27 %;
runtime 2106 -> 1248 s. BUT honest even/odd LDA AUC 0.5275 -> **0.5240**
(unchanged) and the state stage still fails — all GT enrichments ~1.0x, state 5
still 86,631 particles. Reason: concentrating signal into 5 clean coordinates also
discards the 43 directions holding ~59 % of it, so total discriminative
information is flat. **The win is that the POOLING requirement fell 6.3x — run any
mixture/pooling estimator on the neigs=5 embedding, not the 48-component one.**
Also: the reduced solve DOES scale with neigs (279.7 s -> 8.8 s); an earlier note
saying d_tilde=ncols makes it neigs-independent was wrong.

**Highest-leverage item is NOT a scale fix.** It is: calibrate A^-1 against
observed within-group spread and use it in a POOLING estimator — the
measurement-error mixture model is the concrete form, needs no code fix to become
available, and addresses the real constraint (population-level information exists,
per-particle information does not). Second option: refit with far fewer components
in the E-STEP (neigs=5); note d_tilde in the reduced solve is set by `ncols`, not
`neigs`, so that stage does not speed up.

Instruments: `scripts/embed_diagnosis.py`, `scripts/embed_noise_correct.py`.
Score IgG on the BINARY frames 0-19 vs 20-99 contrast, never on frame index as a
continuous coordinate — see [[flex-covariance-tau-metric-is-uninformative]], whose
sub-null ceiling this very likely explains.
