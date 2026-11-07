---
name: multibody-flex-means-region-pca
description: "The user's definition of 'multi-body flex' (2026-09-19): split the molecule into regions, run the covariance PCA in each region, combine the bases -- NOT rigid-body generators; the region-PCA is the nonparametric form, rigid generators its 6-column special case, and a 2-region version + the combine already exist (SIMPLE_COV_COMPLEMENT block, SIMPLE_COV_COMPOSE union embedding)"
metadata:
  type: feedback
---

Stated 2026-09-19 after I had written up multi-body as rigid tangent generators (E1): "my interpretation of multi body
was to split the molecule into different regions and do pca in each of those regions and combine the basis somehow".

**Why:** the user wants to reuse flex_pca as is (learned basis per region, any motion shape incl. occupancy), with the
region split as the only new structure; rigid generators are a parametric special case, not the plan.

**How to apply:** frame multi-body work as region-wise PCA + union basis + joint embedding. Existing hooks:
`SIMPLE_COV_COMPLEMENT=k_out` = a two-block region PCA in the M-step (block 1 envelope window, block 2 sphere minus
envelope, deflation per block on its own window; `simple_flex_pca_pcg.f90:66-70, 2092-2129`); B blocks = per-column
window assignment, per-voxel coupled solve block-diagonal, E-step joint (handles leakage). `SIMPLE_COV_COMPOSE=<dirs>`
= the combine: loads each run's pc MRCs + probe meta, Gram-Schmidt (same-box copied verbatim, listing order), prior
variances rescaled, ONE embedding, no EM, then deconv/placement/states on the union latent
(`simple_flex_pca_em_compose.f90`). Zero-code RNA experiment: B runs with `pcg_mskfile=<region>` +
`SIMPLE_COV_COMPLEMENT=k`, then one COMPOSE run. Combine at the LATENT level (per-region on-axis states; second-level
PCA of the whitened union latent for coupled motions, as relion_flex_analyse does over body poses); never re-PCA the
union volumes. Rank per region and region count by the cross-half gate. Regions: dynamic-domain clustering of the
global covariance (patch co-variation), rigid-consistency of the finals, or seg8. Doc section 4.4b of
`doc/implementation_notes/flex_rna_symmetry_multibody_investigation_2026_09_19.md`.
See [[multibody-flex-design-2026-09-19]], [[flex-pca-basis-composition]], [[flex-pca-complement-block-and-step0]].

**Measured 2026-09-20 (toy harness, 3 regions x {6,9,12} deg, 40.5k ptcls):** region-wise fits beat the global at
matched rank on reproducibility (12/12 vs 4-5/12 pairs >= 0.85), readability (one clean axis per motion,
0.89/0.92/0.93 vs 0.55/0.74/0.92), per-region on-axis delivery (0.97-0.98 vs GT marginal, 92-93% of the difference in
the right region) and the low-SNR rank budget (theta1 at snr 0.01: union 0.56, global rank 4 0.14); lone region
fits LEAK the other regions' motions into their lower columns (region-1 block predicts theta2/theta3 at R^2 0.61/0.69
even in the joint embedding) -> the joint block fit on branch flex-region-pca. Region discovery (regions.py,
co-variation clustering) recovers the toy regions at 91% voxel agreement from the plain global basis; on the RNA the
covariance is smooth (ratio 1.3-1.5) so the split is a partition of a continuum. RNA region runs: 30/30 columns
>= 0.85 (means 0.95/0.99/0.99) vs G10 10/10 at 0.96.

**RNA readout 2026-09-20 (label-free):** per-region on-axis deliveries (3 states per region, halves) put 78-81% of the
end-to-end difference inside the region (global leading axis: 26/30/34%), top-1% concentration 0.81-0.84 (0.66),
even/odd 0.98, states at 4.0-4.1 A; union (`C_regions`) and joint block fit (`B_regions`) agree at 0.86/0.96/0.92. This
is the multi-body readout on the RNA. Recipe: label volume (`regions.py` -> `<tag>_labels_B3.mrc`), block fit
`flex_pca pcg_mskfile=<labels> neigs=10 rec_states=no`, then `flex_pca infile=<cache> rec_states=regions
pcg_mskfile=<labels>` with `SIMPLE_COV_DECONV=0` (else hours in the K ladder); score with `rna_axis_score.py`.
