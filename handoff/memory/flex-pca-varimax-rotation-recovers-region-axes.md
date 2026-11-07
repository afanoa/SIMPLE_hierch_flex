---
name: flex-pca-varimax-rotation-recovers-region-axes
description: "2026-09-20: a post-hoc real-space VARIMAX rotation of the global flex_pca basis (rotation within the span, no labels, no refit, seconds) recovers one clean axis per region on the toy as well as the joint block fit (rank 4: single-column |corr| 0.91/0.94/0.97, own-group R^2 0.84/0.87/0.95, leaks <= 0.05); on RNase P the RANK decides: G10 cache -> region-leaning axes (deliveries 0.48-0.54 in-region), G30 cache -> region-CONFINED axes (columns 0.78-0.85, deliveries 0.61-0.73 in-region, concentration 0.78-0.83, e/o 0.98, 4.0-4.2 A) = the labelled region fits (0.68-0.81) without labels. SIMPLE_COV_VARIMAX=1 validated on fresh runs and cache resumes"
metadata:
  type: project
---

The global fit's SPAN already contains the region-localised columns; the Gram eigenframe (and any PCA
frame) hides them as mixtures. `scripts/flex_pca_verify/region_toy/varimax_score.py` (prototype) and
`src/main/flex/simple_flex_pca_varimax.f90` (pairwise Kaiser sweeps over the window voxels on the
unit-norm columns; Z' = Z D R, precisions R^T D^-1 A D^-1 R, Gamma' = diag(R^T D Gamma D R)).
Numbers: note section 7.12 of `doc/implementation_notes/flex_rna_symmetry_multibody_investigation_2026_09_19.md`.
At snr 0.01 the weakest motion is absent from a rank-4 span but present in rank 12, where the rotation
finds it (R^2 0.50 vs 0.55 for the union of labelled region fits): the earlier "rank budget" advantage of
region fits was RANK, not region support. The block fit keeps leak-free attribution at high rank
(surplus global columns are mixtures the rotation spreads) and the per-region delivery machinery.

**How to apply:** `SIMPLE_COV_VARIMAX=1` on a resume: `flex_pca infile=<run>/flex_pca_embedding.bin
state_axis=q npreimages=3 SIMPLE_COV_DECONV=0 ...` slices the ROTATED column q (the state stage never
reads the basis; the driver rotates z/precisions/eigvals from the cache + the eigenvolumes in the cache
dir, writes flex_pca_varimax_pc*.mrc and flex_pca_varimax_rotation.txt; a cache dir with that file is
already rotated). Next: toy resume of T_G4new2 must reproduce the Python numbers; then the RNase P G10
cache, 10 on-axis deliveries, scored label-free by rna_axis_score.py against B_regions (0.68/0.80/0.81
in-region, e/o 0.98). Related: [[flex-pca-lowdin-frame-fixes-column-reproducibility]],
[[multibody-flex-means-region-pca]], [[flex-region-pca-branch-state]].

**RNase P result (V10_ax1..10, 2026-09-20):** rotation of the ten G10 columns doubles the criterion (28 sweeps,
|diag R| 0.26-0.72); column energy per G10cov region: leaders col 10 -> r1 0.70, col 2 -> r2 0.70, col 7 -> r3 0.61
(index 0.41 -> 0.59, none >= 0.8). On-axis deliveries: two axes per region at 0.48-0.54 of the end-to-end difference
energy in-region (ax10/9 -> r1, ax2/4 -> r2, ax7/8 -> r3), concentration 0.69-0.73 (global axis 0.66), rel diff
0.09-0.125, e/o 0.98. About two thirds of the labelled region deliveries' localisation (C_regions 0.78-0.80,
B_regions 0.68-0.81); the RNase P continuum has no block-confined frame (region finder ratio 1.3-1.5 said the same).
A rotated resume never adopts the cache's deconvolved block (90c90a30b); use SIMPLE_COV_DECONV=0 for states-only reruns.

**Rank-30 result (V30_ax*, 2026-09-20 afternoon):** G30 cache (30 columns, run 53267818 completed) rotated in-job: criterion
6.9e-3 -> 3.5e-2 in 90 sweeps, |diag R| 0.02-0.57 (column numbers do NOT carry over from the Python prototype -- always rank
the job's own flex_pca_varimax_pc*.mrc); index 0.41 -> 0.675, 22/14/4 of 30 columns >= 0.6/0.7/0.8, 10/9/11 per region.
Leaders r1: 9, 8 (0.85, 0.84); r2: 11, 15 (0.78, 0.78); r3: 4, 3 (0.85, 0.81). Deliveries: r1 0.71/0.73, r2 0.62/0.61,
r3 0.71/0.72 in-region (control axis 1: 0.17/0.30/0.43); labelled fits 0.68-0.81. Recipe for a new dataset: global rank 30
+ SIMPLE_COV_VARIMAX=1, rank the rotated columns by region energy, deliver state_axis=q for the leaders; the block fit with
those regions as labels only when a region needs more rank or a cleaner window. Note section 7.14.
