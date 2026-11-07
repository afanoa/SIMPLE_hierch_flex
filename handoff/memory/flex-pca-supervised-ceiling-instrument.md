---
name: flex-pca-supervised-ceiling-instrument
description: "The supervised ceiling (2-fold cross-fitted logistic regression on published labels over raw latents) tells you whether to fix the REPRESENTATION or the STATE STAGE; on 10076 the gap closed to 0.005"
metadata:
  type: reference
---

`ceiling.py` (standing script, vendored in doc/implementation_notes/flex_em_scripts/):
2-fold cross-fitted multinomial logistic regression trained ON published labels over the raw
latents, scored held out. No unsupervised state stage can beat a classifier handed the answer
key, so the delivered-vs-ceiling gap LOCALIZES the loss:

- delivered ~= ceiling => clustering is solved; only representation changes can help.
- delivered << ceiling => the state stage is leaving science on the table.

This inverted the Ribosembly-era "state stage is the bottleneck" read for 10076, then closed it:
plain GMM 0.1802/0.3357 (gap 0.080) -> MCFA mixture prior -> direct K=8 delivery **0.2620/0.3738
against ceiling 0.2671, gap 0.005**. Positive-contrast project best: 0.2681/0.3905.

**Why:** without it you cannot tell a bad basis from a bad partitioner, and the round wasted
effort in both directions before it existed. It is also what makes "derived-neigs is a bound"
legible — see [[flex-pca-derived-lp-is-two-branch]].

**How to apply:** requires published labels, so it works on 10076 / Ribosembly / CryoBench and
NOT on unlabeled real data (10049). Two traps measured: (1) trust AMI over ARI when the arms
disagree on K-granularity, per the scorer's own doc; (2) rerun drift is dataset-specific —
10076 drifts ~0.01 but **Ribosembly's rerun band is ~0.85-0.97**, so single-pair deltas there
prove nothing and any Ribosembly A/B needs n>=3 per arm or paired seeds. Extends
[[flex-pca-downstream-scoring-protocol]].
