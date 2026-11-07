---
name: rnasep-2d-fine-3d-coarse-complementary-2026-10-14
description: VIEW-BALANCED verdict on RNase P (2026-10-14 02:35): windows over sub-classes ordered within each class along the synchronised 2D-latent coordinate carry the FINE motion (20-12 A 2.1x, 12-8 A 2.8x the within-class random floor, nothing coarse), the pooled flex_pca lp 10 latent carries the COARSE motion (117-30 A 1.8x, 30-20 A 1.6x, nothing fine); the two latents are complementary by band; a map-learnt basis does not span the deviations; next lever = per-window poses along the fine coordinate
metadata:
  type: project
---

Design (`run_window_rec_bal.sbatch`, keep files from the balanced writer in the trajectory scripts): 8 disjoint
windows per family, each with the same share of every class's members, sub-classes ordered along the coordinate
inside each class; within-class random ordering = floor (0.23/0.16/0.07/0.05 per band); tilt-only partition =
pure view effect (0.64/0.68/0.43/0.25, see [[view-imbalance-inflates-between-state-metric]]). Results: synchronised
2D axis 1 0.23/0.26/0.14/0.14; pooled lp 10 latent PC1 0.43/0.25/0.08/0.06; unregularised span-coefficient PC1
0.33/0.19/0.08/0.05; state discriminant 0.32/0.19/0.07/0.04. The pooled flex latent beats my own regression
coefficients (regularisation matters). Notes section 16.3 of
`doc/implementation_notes/coupled_2d3d_math_formulation_2026_10_13.md`; figure `viz_window_families_balanced.png`.

Also closed tonight: a 3D basis learnt from 45 window/state maps (plain or cross-half PCA) explains 0.28 of the odd-
half sub-class deviations vs 0.53 for the 16 A sphere basis (worse in every shell); the lp 10 sphere basis spans the
sub-class deviations to ~20 A (0.60 at 29 A, 0.43 at 23 A, 0.28 at 20 A) and the deviations reproduce only to ~17 A.

**Why:** the 3D basis is fitted at consensus poses on particles and is resolution-limited to ~20 A; the 2D sub-class
coordinate is built from view-local averages and keeps the fine structure; pooling the 3D latent over 2D sub-classes
denoises it without adding resolution.

**How to apply:** couple at the LATENT level per sub-class: [pooled flex latent (coarse), synchronised 2D axes
(fine)]; build states/trajectories view-balanced with the random-within-class floor; to push the fine part into a 3D
basis, refine poses per window along the 2D coordinate first, then fit the basis to pose-refined windows. Never read
an unbalanced between-state number again. Related: [[rnasep-subclass-window-trajectory-2026-10-14]],
[[rnasep-sync-od-coordinate-2026-10-14]], [[rnasep-basis-span-finding-2026-10-13]].
