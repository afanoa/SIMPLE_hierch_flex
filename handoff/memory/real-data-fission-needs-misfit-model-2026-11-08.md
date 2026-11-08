---
name: real-data-fission-needs-misfit-model-2026-11-08
description: "TRAP (RNase P 2026-11-08): sat5's fission test accepts per-particle misfit splits (log-likelihood / amplitude) even with nuisances profiled and groups stratified; needs a per-particle band-wise sigma2 in the whitening and the fission draw"
metadata:
  node_type: memory
  type: project
  originSessionId: 1218bb2b-a21a-4d25-9e48-c40ddd89ceff
  modified: 2026-11-08T05:19:37.004Z
---

The tree-latent decisive rerun on RNase P, 2026-11-08 (batch2 `handoff/tree-latent/STATUS.md` addendum 3) changed two things:
- the amplitude and the 12 nuisances were profiled out of the E-step;
- sat5's four groups were formed inside 3 x 3 strata of the y1 amplitude and the y1 log-likelihood.

The 692k root's split was still DISCRETE (z_mid 43 -> 20).
- It tracks log-likelihood (|AUROC - .5| .42) and amplitude (.30).
- Its 76k child is closed-depleted (.007 against .013).
- Random-split controls stay at z_mid 0-1.
- Profiling cost the coordinates only .005 of closed AUROC.

**Why:** fission (y1 = y + Z, y2 = y - Z) cannot remove a structured per-particle residual, because it is part of y and sits in both
copies. A per-particle noise level that differs from the model's sigma^2 also correlates the copies, by e / (2 + e). The rnasep-misfit
session found such a residual: band-wise residual power SD about .29, reliability .96, clustered by micrograph.

**How to apply:**
- On real data, give any fission-tested split a per-particle band-wise sigma^2, in the E-step whitening and in the fission draw.
- Report each accepted split's association with log-likelihood and amplitude before reading it as a state.
- Toy fixtures do not show this trap.

Related: [[tree-latent-verdict-2026-11-07]], [[data-fission-merge-test-verdict-2026-11-07]].
