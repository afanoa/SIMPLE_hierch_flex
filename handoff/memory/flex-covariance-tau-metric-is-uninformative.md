---
name: flex-covariance-tau-metric-is-uninformative
description: the IgG tau (frame-index) metric in flex_covariance_validate.py cannot decide anything -- its GT ceiling sits below its own null
metadata: 
  node_type: memory
  type: project
  originSessionId: 205ad3f4-4a2f-4250-8466-860b63608957
  modified: 2026-07-31T23:54:37.919Z
---

Measured 2026-07-31. `flex_covariance_validate.py --tau-mode igg` scores against
tau = floor((pind-1)/1000)/99. This metric is incapable of ranking two runs:

- its "frame-averaged best |corr|" is a MAX over ncomp components with only 100 frames.
  Permutation null (shuffle tau across frames, 2000 draws, ncomp=20): **mean 0.21,
  p95 0.29**. The script reports this statistic with NO null, and its printed advice
  ("frame-avg high ... the direction is usable") is therefore misleading.
- the achievable ceiling is corr(GT-volPC1, frame_index) = 0.11 -- BELOW the null.

So even a perfect method scores ~0.11 against a 0.21 noise floor. Measured: run81 0.1227
(p=0.97) and run86 0.1950 (p=0.58), both not significant.

Corollary for the per-particle number: tau is CONSTANT within a frame, so any genuinely
frame-dependent signal must survive frame-averaging (it cannot shrink). run81's per-particle
CV corr 0.0916 -- nominally 6.5x its overfit floor -- does NOT cohere at frame level, so it
is not conformational signal. It is not raw index leakage either (corr(z11, pind) = -0.036).

**Why:** two full 25-min runs were scored on this before the null was checked, and an
apparent "59% improvement" turned out to be noise-vs-noise.

**How to apply:** use `scripts/flex_gt_volpc1.py --gt-dir <100 CryoBench GT volumes>` as the
primary metric. Those GT volumes are NOT on this machine as of 2026-07-31 -- obtain them
before running any further flex_covariance A/B. See [[flex-covariance-recovar]].
