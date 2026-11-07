---
name: mask-core-null-calibration-no-go-2026-11-01
description: audit 6.1 / 6.2 prototype verdicts - a noise-calibrated mask core covers the whole molecule (between-state variance is significant everywhere); CV-SVD factors over-count without a share threshold
metadata:
  type: project
---
2026-11-01 (math reworks, mask_disc_proto.py on the base_qtoy2 / base_qtoynull / base_r131 root flex half maps):
- v = sum_j w_j (E_j - E)(O_j - O) recomputed exactly (corr 1.000 with mask_variance.mrc on qtoy2); its null sd from half-difference
  surrogates gives peak z 170-3500 and z > 3 over 75-100 % of the sphere on every root, null included. The flex states differ
  reproducibly by small amounts everywhere (scale / B / pose residuals), so a significance core does not localise: 6.1 is a NO-GO;
  the 0.2 x peak rule (localisation) stays.
- Weighted SVD of a component's voxel-by-state deviations: 2-4 factors reproduce across halves (r 0.9-0.999, n_eff 7-37) and 2-5
  singular values beat a noise floor, including 1-4 % factors; the 2-means disagrees in both directions (splits a 0.93/0.04 component,
  calls 0.77/0.20 one factor). Partial go: SVD should replace the 2-means but one share threshold on the singular values stays.
**How to apply:** do not propose noise-calibrated cores or reproducibility-only factor counts for the tree again; any SVD port keeps a
share threshold. Related: [[hetero-tree-math-reworks-setup-2026-11-01]], [[tree-design-feature-freeze-2026-10-25]].
