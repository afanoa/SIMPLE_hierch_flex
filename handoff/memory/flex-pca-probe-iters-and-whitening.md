---
name: flex-pca-probe-iters-and-whitening
description: n_probe_iters nearly doubles flex_pca basis quality and is now ON by default at 5 (commander, not the struct); and the latent must be whitened by sqrt(eigvals) before any distance is computed
metadata: 
  node_type: memory
  type: project
  originSessionId: 9671e865-e13f-4298-838f-d873182e5df9
  modified: 2026-08-12T18:36:04.367Z
---

Two settings dominate every other flex_pca tuning knob. Measured on Tomotwin-100
(2026-08-04, 100 species x 1000 particles, SNR 0.01) against a ground-truth ceiling.

**1. `n_probe_iters=3` — worth more than everything else combined.**

**It is now ON by default at 5** (verified 2026-08-12). The default lives in the *commander*
(`simple_commanders_flex_pca.f90`, a `cline%set` before `params%new`), NOT in the `parameters`
struct, which still declares `0`. Reading the struct alone tells you "off" and is wrong. It is
also an upper bound, not a count: the probe stops early at principal-angle cosine >= 0.97.
The measurements below were taken when it was still off by default.

| box48 lp26 ncols256 | probe off | `n_probe_iters=3` |
|---|---|---|
| basis capture vs achievable ceiling | 42.2% | **80.5%** |
| species separation | 76.1% | **94.6%** |
| species below 50% recall | 2/100 | **0/100** |

A full ladder of `d_tilde` (64/96/128/177/220), `ncols` (128/256), `box_crop` (32/48/64) and
solver ablations topped out at **49.2%** efficiency. One default-off flag reaches 80.5%.
Cost: `n_probe_iters` extra full-data passes (~400 s each at box_crop=64), so ~1.5x runtime.

**2. Whiten by `sqrt(eigvals)` before computing any distance.** Latent eigenvalues span
~61x, so unweighted Euclidean drowns the signal in the noisiest components: 61.8% -> **86.0%**
separation, species below 50% recall **43 -> 0**. Free — eigenvalues are in the cache.
Necessary only when probe is OFF: with probe on, Euclidean (95.0%) slightly *beats* whitened
(94.3%), because the 61x spread was itself a symptom of poor conditioning.

**Other parameter behaviour:** `d_tilde` has an interior optimum near 177 (96 -> 42.2%,
128 -> 45.5%, 177 -> 49.2%, 220 -> 34.8%) — both directions hurt. Lower resolution helps
(box48 > box64 by ~5 points). `ncols` 128 -> 256 is worth ~4 points. `SIMPLE_COV_DIRECTSVD=1`
is catastrophic (29.2% efficiency, 22.8% accuracy) — the reduced solve's span-and-re-estimate
step is essential.

**How this was measured:** `scripts/flex_pca_basis_capture_gt.py` — captured between-species
variance as a ratio against the best achievable subspace of the same dimension, with a random
subspace as floor. `corr(efficiency, accuracy) = 0.923` over 7 configurations, so it is a
validated predictor and not just descriptive. Note this metric was RETIRED in the IgG work;
it works here because it is a ratio against a *measured ceiling* with a demonstrated link to
a downstream outcome.

See [[tomotwin100-project-setup]], [[flex-pca-embedding-cache-format]].
Full record: `doc/for_developers/ideas/tomotwin100_investigation_log.md`.
