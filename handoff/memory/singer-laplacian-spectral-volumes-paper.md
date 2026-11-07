---
name: singer-laplacian-spectral-volumes-paper
description: "Reference (user pointed to it 2026-10-15): Moscovich, Halevi, Anden, Singer, 'Cryo-EM reconstruction of continuous heterogeneity by Laplacian spectral volumes', arXiv 1907.01898 / Inverse Problems 2020: particle graph from covariance-PCA coordinates (16-px images) -> Laplacian eigenvectors -> volumes expanded in that basis (spectral volumes) by a generalized tomographic least squares; poses/CTF known, simulations only; = the fixed-graph cousin of klem"
metadata:
  type: reference
---

arXiv https://arxiv.org/abs/1907.01898 (Inverse Problems 2020). Method: (1) images downsampled to ~16x16, 3D mean +
covariance estimated, PCA eigenvolumes, per-particle PCA coordinates beta_s; (2) Gaussian-kernel or kNN graph on the
beta_s (no view-specific embedding needed because beta_s is already a 3D coordinate), graph Laplacian eigenvectors
phi^(l) = coordinates on the conformation manifold; (3) volume of particle s = sqrt(n) sum_l alpha^(l) phi_s^(l); the
"spectral volumes" alpha^(l) solve sum_s |y_s - sum_l (P_s alpha^(l)) phi_s^(l)|^2 by CG on the normal equations with
NUFFT projections, no explicit regularisation; (4) poses and CTF assumed known, white noise SNR ~0.033, two simulated
ion-channel sets (1-D spin, 2-D stretch, 10k images, 108^3): FSC improves with the number of spectral volumes, the
1-D case reconstructs better than the 2-D one; (5) stated limits: covariance step impractical beyond ~25 px,
downsampling limits the manifold resolution, perfect poses assumed, continuous heterogeneity only.

Relation to klem (2026-10-15): same decoder idea (volumes = smooth functions on a latent manifold) but the graph is
FIXED from the linear covariance latent (no E-step feedback); on toy10 the per-particle linear latent is too noisy
(15-NN 0.17-0.25) for a fixed graph, which is why klem re-estimates the latent from the images each iteration. Its
graph construction is used in klem as the --init spectral (kNN graph on the contrast-free lp 12 latent, eigenvectors
2-3 of the normalised Laplacian). Related: [[klem-kernel-latent-em-2026-10-15]].
