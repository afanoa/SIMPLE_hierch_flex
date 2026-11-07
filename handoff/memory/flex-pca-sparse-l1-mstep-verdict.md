---
name: flex-pca-sparse-l1-mstep-verdict
description: "2026-09-20: the sparse (L1/ADMM) M-step (SIMPLE_COV_L1=t, commit 9027314f3) works mechanically but localises columns SLOWLY (rank 12: 0.52 -> 0.71 at 4 iters -> 0.78 at 8) because the M-step only shrinks the current frame and the frame turns through the E-step; t=2 thresholds the weakest motion out of the basis; everything it gains at t<=1 the varimax rotation gives exactly and for free -- keep it opt-in, do not tune it"
metadata:
  type: project
---

Design: ADMM lasso continuing the plain PCG solve (x-update = same CG with rho on the Tikhonov term,
z = soft threshold at t x column rms on its window, debias refit on supp(z)); Lowdin frame forced on.
Toy (snr 0.05): support kept 20/9/4 % of the window at t = 0.5/1/2; primal residual 0.6-2.6 -> 0.1 in
4 rounds; rank 4 t=1 localisation 0.63 -> 0.68 (single-column |corr| 0.64/0.77/0.91 -> 0.71/0.91/0.91);
rank 12 t=1 0.52 -> 0.71 (3 of 12 columns >= 0.8), 8 iterations 0.78 (6 of 12) with reproducibility
dropping (10 of 12). t=2: theta1 R^2 0.85 -> 0.39 (the weak motion is thresholded away).

**Why it stalls:** for fixed latents the penalised regression has no rotation freedom; the sparse
frame is reached only through E-step/M-step alternation (SPCA's Procrustes step is our
orthonormalisation), a few degrees per iteration, and iterations beyond 4 overfit cross-half.

**How to apply:** use [[flex-pca-varimax-rotation-recovers-region-axes]] for localisation; the L1 stays
as an opt-in experiment (knobs SIMPLE_COV_L1, _RHO, _ROUNDS, _DEBIAS). Its real by-product is
[[flex-pca-lowdin-frame-fixes-column-reproducibility]].
