---
name: rnasep-2d-subclasses-pooled-vs-particle-3d-latent-2026-10-13
description: "RNase P: the 'rich 2D vs blob 3D' gap is averages vs particles. Fortran cls_expansion sub-classes (clsexp_circle_prior) are real (direct averages 11x null in band, 1.30x z 15 at 15-8 A); the G10 3D latent separates them pooled (sub-class means 18x shuffled) but not per particle (LDA 0.39 vs 0.34)"
metadata:
  node_type: memory
  type: project
  originSessionId: 34b66a02-03d2-42df-91f8-915afbf80120
  modified: 2026-10-13T14:25:04.774Z
---

2026-10-13, classes 486/477/187/350/339 of the 932k set (all members, 13-16k each; data
`/mnt/beegfs/elmlund/afan/rnasep_coupled/latent2d3d/cls_full`, scripts `scratchpad latent2d3d/fortran_vs_3d.py`,
figures `doc/implementation_notes/rnasep_coupled_figs/{fortran_vs_3d,pooled_3d}_cls486_full.png`). 2D = the Fortran
sub-class labels of `flex3sets/rnasep/clsexp_circle_prior` (100 per class); 3D = G10 flex_pca coordinates (10-D,
envelope pcg, consensus poses) of the same particles.

- Direct CTF-Wiener averages of the top-6 sub-classes from the images reproduce the Fortran conditional averages
  (arm displacements): class 486 11x shuffled in 200-10 A (biased band), 1.30x (z 15.5) in 15-8 A (beyond the fit
  band, unbiased). rand284k-only members (4-5k/class) gave 2.5-4x and 1.07-1.14x on all five classes.
- 3D latent on the same sub-classes: pooled between/within trace 0.25 vs 0.006 shuffled (z 204); sub-class means on
  the first two discriminant axes spread 0.291 vs 0.016; per-particle 5-fold LDA 0.39 vs 0.34 (majority 0.34). Same
  pattern in all five classes (z 30-77 on 4-5k).
- **Reading:** the 3D latent carries the motion the 2D sub-classes show, but per particle it is noise-dominated; 2D
  looks rich because it is shown as sub-class AVERAGES of hundreds of particles. Coupling = pool the 3D latent over
  the 2D sub-classes (tied placement, posterior-mean embedding), i.e. the coupled flex_pca `state_placement=subclass`.
- Open: shared-noise share of the pooled separation (needs a split-half or parity control); the per-particle 2D latent
  itself needs SIMPLE_FLEXCLS_DUMP from the coupled-2d3d tree (only on the RAG machine).
- A numpy per-class PPCA on these images put its leading direction on the background level (|CTF|~0.1 at the lowest
  shells inflates it); even with an 80-20 A band it reached 2x null vs Fortran's 11x: do not stand in for the Fortran
  model ([[use-fortran-outputs-not-numpy-reimplementations]]).

Related: [[rnasep-stacks-are-class-prealigned]], [[cls-expansion-vs-flex-pca-latent-diagnosis-2026-10-09]],
[[coupled-2d3d-realistic-regime-result-2026-10-10]].
