---
name: rnasep-lp25-flex-basis-negative-2026-10-14
description: RNase P lp=25 flex_pca (24 comps, sphere support) is NOT a better coarse-motion basis; its span margin over a rotated null is 0.01-0.08 at 117-59 A (lp10 0.2-0.4), its pooled PC1 is 81 % class identity (view), balanced windows at/below the random floor; main run stalled after the cache (cancelled)
metadata:
  type: project
---

Test D of the 2026-10-14 battery (`rnasep_coupled/tests/lp25`, notes section 19.4): flex_pca lp=25 neigs=24 on the
692k refined RNase P particles (working box 64). Span test of the absolute sub-class deviations: lp25 explains 0.89 /
0.87 / 0.69 of the 117 / 78 / 59 A shells but a randomly rotated copy explains 0.88 / 0.84 / 0.61 (lp10: 0.58 / 0.71
/ 0.68 vs null 0.38 / 0.49 / 0.40), dead beyond 40 A. Pooled per sub-class, its PC1 has class-identity share 0.81
(lp10 0.27); view-balanced windows along it reproduce at 0.30 / 0.14 (117-30 / 30-20 A) vs lp10 0.43 / 0.25 and the
random floor 0.23 / 0.16.

**Why:** a lower band buys smoothness, not structure: the leading low-resolution modes are generic (any rotation of
them fits as well) and pose nuisance dominates the latent.

**How to apply:** do not lower the flex band to chase large motions; the lp10 basis already holds the coarse part.
The main lp25 run stalled for an hour after writing its embedding cache (between the deconvolution step and the
deconvolved UMAP); the cache (version 5, 24 columns) can be pooled with pool_cache.py without the coordinates file.
Related: [[rnasep-2d-fine-3d-coarse-complementary-2026-10-14]], [[rnasep-per-window-poses-gold-standard-2026-10-14]].
