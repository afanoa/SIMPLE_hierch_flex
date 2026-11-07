---
name: mac-lsqr-vstep-proxy-invents-soft-gain
description: the Mac Fourier-grid lsqr V-step (states from a few sub-class averages) rewards soft assignment by regularisation; real particle reconstructions on batch2 showed soft (tau 0.05) as the WORST variant on every band
metadata:
  type: feedback
---

The Mac proxy for the alternation (states solved by least squares on the Fourier grid from a handful of sub-class
averages, `assign_cross_metric.py`) gave soft Sinkhorn assignment the largest cross-half gain of the night
(0.52/0.55/0.45 -> 0.86/0.85/0.69). On batch2, with real reconstructions from 20-80k particles per state, the same
soft assignment (v5, tau 0.05, mean max responsibility 0.51) was the worst of all variants on every band
(0.31/0.29/0.15/0.09 vs 0.43-0.48 coarse for hard variants, tracked start 0.36/0.58/0.38/0.23).

**Why:** the lsqr V-step is under-determined (few averages per state), and soft weights regularise it; a real
reconstruction has nothing to regularise, and mixing just halves every state-minus-mean deviation while the half-map
noise stays, so the even/odd SNR of the deviations collapses. The metric (reproducible between-state variance) is the
deviation SNR and penalises softness by construction.

**How to apply:** never rank assignment variants on the Mac lsqr proxy; use it only for ordering hard variants, and
confirm any claimed gain with the cluster reconstructions before reporting it. See [[rnasep-alternation-verdict-2026-10-14]].
