---
name: flex2d-must-be-alignment-free
description: user direction 2026-09-30 - flex2D heterogeneity must work WITHOUT alignments, orientations or 3D information
metadata:
  type: feedback
---
2026-09-30, after gold-standard alignment still gave "weird" classes: "we already have a 3D implementation that works
fine. I want this to work without alignments. Perhaps this model needs to be more sophisticated" and then "I meant we
shouldn't use orientations or 3D information". So no reference alignment, no 2D poses (not even abinitio2D's), no
abinitio3D/class-orientation lifting. Heterogeneity must come from rotation/translation-INVARIANT representations of
the particle images (polar/angular-Fourier invariants, steerable PCA, bispectrum-type terms) and a richer latent model
on them (mixtures, manifold structure).

**Why:** every failure mode of the current flex2D (reference bias/noise alignment, smeared classes, population drain,
view/state mixing in the latent) comes from its alignment loop; 3D heterogeneity is already covered by flex_pca.

**How to apply:** do not propose fixed-pose, class-orientation or 3D-lifting variants; judge alignment-free models on
the toys against truth (state separability of the invariant embedding, unsupervised purity/recall). CTF correction per
particle is allowed (not orientation information). Related: [[flex2d-findings-2026-09-21]].
