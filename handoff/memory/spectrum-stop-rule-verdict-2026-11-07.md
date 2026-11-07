---
name: spectrum-stop-rule-verdict-2026-11-07
description: "VERDICT 2026-11-07: counting flex-latent directions above a BBP noise edge is no stop rule (nulls continue on nuisance directions, N sets the 692k stop, RNase P closed spike <= 0.13); r+1 is a lower bound, not an upper one; the n/sqrt(N) rare-state law IS the BBP law (p_eff 0.16-0.41 of the masked Shannon count)"
metadata:
  node_type: memory
  type: project
  originSessionId: c6b87431-7325-4b60-ad85-a72f7e8559d4
  modified: 2026-11-07T17:30:24.075Z
---

Task from the coordinator (local_4b539ac3): can a node latent's eigen-spectrum decide how much to split, with r + 1 bounding the states?
The verdict and tables are on batch2 in hetero_tree_20261028/handoff/spectrum-stop/ (VERDICT.md, report.txt, nodes.tsv). The work dir
is math_reworks/spectrum (all_v2.sbatch reruns everything in about 6 min).

- **Statistic.** rho_j = generalized eigenvalues of (Cov z_MAP, diag(sqrt a_comp) E[Q^-1] diag(sqrt a_comp)) minus 1. Edge (a) is
  2 sqrt(g) + g, with g = (4pi/3) kstop^3 * support / N.
- **Nulls.** qtoynull 1.33 and the detlim null 1.03 sit above the 0.77 edge. Their direction is cross-component and Gaussian; it is
  not pose, shift, contrast, B, SNR, junk or view. micnull is 0.82 against 0.68. IgM fix and null twins have the same spectrum
  (shift-error tangents plus contrast).
- **692k.** rho_1 stays at about 0.8 from depth 3 to 7 while the edge grows as 1/sqrt(N), so N sets the stop. The closed spike has
  median 0.009 and max 0.13.
- **Rare state.** detlim conf18 is latent direction 2, below the nuisance direction. The BBP fit gives p_eff 1,000-2,600, so
  c = sqrt(p_eff)/s^2 is about 14-20, against the measured 16-17.
- **Site focus.** IgM T1 true-arm masks raise the arm spikes 15-140x (to 0.97-1.73), but the continuous null arm also spikes
  (0.54-0.63).

**Why:** the coordinator wanted the latent to play a principled role in splitting decisions.
**How to apply:** do not propose spectrum or rank counts as a stop or merge rule, or K <= r + 1 as a state bound. A spectral criterion
would need pose removal by per-particle tangent profiling (basis deflation fails: [[pose-tangent-deflation-verdict-2026-11-07]]), a full
half-split noise covariance and a focused site latent, and still cannot tell discrete states from motion. Related: [[rare-state-detection-limit-2026-10-29]], [[igm-focused-flex-mask-is-bottleneck-2026-11-07]],
[[tree-charts-coordination-gap-verdict-2026-11-06]], [[global-chart-coordinate-verdict-2026-11-06]], [[flex-embedding-precision-is-map-sampling-precision]].
