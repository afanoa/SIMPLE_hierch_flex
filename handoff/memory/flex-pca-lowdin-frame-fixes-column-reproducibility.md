---
name: flex-pca-lowdin-frame-fixes-column-reproducibility
description: "2026-09-20 toy: the Gram-eigenframe orthonormalisation re-rotates near-degenerate columns each EM iteration, differently per half, so matched-column cross-half cosines measure the FRAME; SIMPLE_COV_ORTHO=lowdin (X G^-1/2, columns keep identity) took rank 12 from 4/12 to 12/12 pairs >= 0.85 at unchanged prediction -- the 'more rank, fewer reproducible components' pattern is at least partly this artefact"
metadata:
  type: project
---

Measured on the region toy (snr 0.05, union window, 4 EM iterations), branch flex-region-pca 9027314f3:
rank 4 eigenframe 2 of 4 pairs >= 0.85 (0.78/0.81 for the two weak columns) -> Lowdin 4 of 4 (0.92-0.99);
rank 12 eigenframe 4 of 12 (0.58-0.97) -> Lowdin 12 of 12 (0.95-0.98). CV R^2 of the true angles and the
localisation index unchanged (a rotation). Run-to-run noise on this toy is 1-4 % in eigenvalues, far below.

**Why:** `orthonormalize_representatives` builds u_q = sum_i X_i V_iq / sqrt(lambda_q) from the Gram
eigenvectors (variance-sorted). Columns of comparable norm are a degenerate eigenproblem, so each half's
frame turns arbitrarily each iteration and the merge's matched cosines drop although the SUBSPACE
reproduces. The Lowdin frame U = X G^-1/2 is the closest orthonormal set to the columns as given.

**How to apply:** the TRPM4 neigs=20 collapse ([[flex-pca-optimal-iters-depends-on-rank]]) and the cnga1
/ 7.2 rank-collapse readings need re-measuring with `SIMPLE_COV_ORTHO=lowdin` before any rank rule is
built on them. Lowdin is automatic when SIMPLE_COV_L1 > 0; it falls back to the eigenframe whenever the
set is truncated (init: representatives > rank). Related: [[flex-pca-varimax-rotation-recovers-region-axes]],
[[flex-pca-sparse-l1-mstep-verdict]].
