---
name: flex-covariance-basis-defect
description: "flex_pca basis capture — the 0.001-vs-0.449 defect was REVERSED on 2026-08-02 (0.719 vs 0.686); capture is retired as a scoring metric and the real gap is per-particle conditioning"
metadata:
  node_type: memory
  type: project
  originSessionId: 0f67a542-3f13-4d0e-b048-ed1b3d3ee8c8
  modified: 2026-08-03T15:29:24.934Z
---

**SUPERSEDED HEADLINE.** This note previously recorded (2026-08-01) that the covariance basis
captured **0.0010** of the true conformational direction against RECOVAR's **0.4489**, i.e. that
the basis contained no signal. That is no longer the measurement.
`doc/for_developers/ideas/flex_covariance_state_and_directions.md` (2026-08-02, table at the top)
re-scores the same run 31 files and reports the reverse, in SIMPLE's favour:

| | SIMPLE (box64, d̃=128) | RECOVAR box128 |
|---|---|---|
| GT-subspace basis capture (lp=15, null 0.0009) | **0.719** | 0.686 |
| best single eigenvolume | **0.865** | 0.844 |
| frame-averaged latent capture | 0.968 | 0.969 |

That doc also explicitly **retires the metric** ("stop scoring basis work with capture"), so
`scripts/flex_basis_capture.py` is no longer the instrument to rank configurations by.

**What the real gap is now:** not the subspace, but the *per-particle* estimator — SIMPLE's
coordinates carry ~3x less usable information (20% vs 60% of latent variance conformational;
1 vs 8 components resolved at SNR>1). It is conditioning and ordering of the basis, not
inference and not the covariance columns. The classic knobs are documented as exhausted:
column regularization, right-kernel width, `lp`, box 128, sharpening, `ncols`/`neigs`.

**Why this note still matters:** the durable lesson is the failure MODE, not the number. Four
"decisive" diagnostics were artifacts and three memories were written from them. Half-set
difference reproducibility is scale-invariant (0.87 coexists with state correlation 0.989); the
cached `precision` block is the MAP *sampling* precision `Pi = A Gtil^+ A`, not the posterior, so
law-of-total-variance and obs/post numbers from it are predicted exactly by `r(2-r)` and
`1/(1-r)`; `resid_energy` and `resid_mean_energy` are in different units.

**How to apply:** be suspicious of any flex_pca "improvement" measured only by state
correlation, difference reproducibility, or swing coverage — all three move the way scatter
moves them. Confirmed 2026-08-03: `NORMW` and `EXCLUSIVE` both looked like wins on those metrics
while every point of gain tracked a proportional LOSS of flexible-mask half-map FSC-AUC
(0.8068 -> 0.7803). Demand a noise-robust metric that separates conformational spread from
scatter. See [[igg-test-project-setup]] and [[flex-covariance-tau-metric-is-uninformative]].
