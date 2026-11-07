---
name: flex-pca-state-weights-are-gmm-not-kernel
description: "flex_pca state weights are tied-covariance GMM responsibilities (opt-out), NOT the Epanechnikov kernel whose log lines also appear; overlap ~2 of 3 states per particle after a 1e-3 responsibility floor"
metadata: 
  node_type: memory
  type: project
  originSessionId: 991559ff-d34f-451c-ad5a-ac1ba95b88c7
  modified: 2026-08-10T15:50:27.012Z
---

The weights that reach the flex_pca state reconstructor are **tied-covariance Gaussian mixture
responsibilities**, not Epanechnikov kernel weights. `l_gmm = .not. cov_env_flag_off('SIMPLE_COV_GMM')`
— on by DEFAULT, opt out with `SIMPLE_COV_GMM=0`.

**The trap:** the Epanechnikov loop still runs and still logs
(`>>> FLEX_PCA Epanechnikov kernel regression, posterior-precision distance; chi2 median=...`),
so a run log shows BOTH mechanisms. The kernel output is only a diagnostic plus `dist_out`, which
feeds `cv_select_bandwidths` and the `SIMPLE_COV_EXCLUSIVE` path; `gmm_state_weights` then
overwrites `weights` outright. Reading the Epanechnikov line as the assignment rule is wrong.

**Thresholding is real and load-bearing:** `RESP_FLOOR = 1.d-3` zeroes small responsibilities and
renormalises each particle to unit total mass. Not cosmetic — `insert_plane_oversamp_multi_scaled`
compacts live states per particle and assumes most are exact zeros. Dense responsibilities produced
NO state-reconstruction output after 35 minutes where the compact path takes ~1.

**Overlap:** the GMM logs `mean active states per particle` directly — measured 1.78 and 2.03 of 3
states on IgG 30k. This equals the sum of per-state `frac_in_support` in the occupancy table, so the
two reports are consistent views of the same thing. Hard `labels` are a separate argmax over these
weights and imply a disjointness the maps do not have (state 3: 436 hard-labelled particles but
neff=3157).

`SIMPLE_COV_NORMW` and `SIMPLE_COV_EXCLUSIVE` are both gated `.and. .not. l_gmm` — under the GMM
they are redundant, since it already renormalises to a partition of unity.

**Why:** two assignment mechanisms live in one routine and only one is active; the inactive one is
the louder in the log.
**How to apply:** to tell which ran, grep for `GMM tied-covariance responsibilities` /
`GMM mean active states per particle`, not for the Epanechnikov line. See
[[flex-pca-deconvolved-mixture-works]] for the deconvolution scale (`SIMPLE_COV_DECONV`, off by
default) that sits on top of this mixture.
