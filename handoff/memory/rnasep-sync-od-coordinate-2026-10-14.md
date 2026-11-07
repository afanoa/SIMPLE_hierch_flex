---
name: rnasep-sync-od-coordinate-2026-10-14
description: one global continuous coordinate from the per-class 2D latents by orthogonal-group synchronisation (image correspondences, Procrustes, spectral O(d) sync, 3D-latent gauge); validated on a synthetic continuum (rank 0.90); on RNase P it is view-free, separates tracked states 6x null, couples to the 3D latent within class, and its quantile bins reconstruct as an ordered ladder 1.5-2.7x a random-frame control but 2-5x below the tracked sub-class states (per-particle noise limit)
metadata:
  type: project
---

`scripts/cls_expansion_proto/subclass_pool/subclass_sync_od.py` (+ `sync_coupling_test.py`, `synth_sync_test.py`,
`run_bin_rec.sbatch`), results in `doc/implementation_notes/coupled_2d3d_math_formulation_2026_10_13.md` section 14.
Per class: sub-class mean latents, weighted PCA to d dims, unit RMS; correspondences between near-view classes from the
sub-class deviation images (the tracking's Hungarian; wide and antipodal tiers optional); orthogonal Procrustes per
pair with a permuted-correspondence null; spectral O(d) synchronisation on the largest connected component;
translations by least squares; gauge from the pooled 3D latent. Per-particle global coordinate for 329k (near-only,
44 classes) or 495k (tiered, 80 classes) particles.

Facts: the view graph with near pairs only (<= 12 deg) has 11 components on this data, so always check components;
the tiers connect everything but align worse (cycle consistency 0.61 vs 0.35). Synthetic known truth: sub-class
|Spearman| 0.90, particle |r| 0.33 of a 0.45 ceiling. Real: tracked-state F 1.59-1.77 vs 0.28 permuted, within-class
coupling to the 16-dim pooled 3D latent R2 0.12 vs 0.04 permuted (axis 1: 0.19 vs 0.05), latent k-means groups F 0.27
vs 0.13. Bins (8 quantiles): near44 0.18/0.16/0.08/0.10 vs random frames 0.10/0.10/0.07/0.06 per band
(117-30/30-20/20-12/12-8 A); extreme bins carry opposite reproducible deviations, middle bins none (a continuum);
tracked sub-class states 0.36/0.58/0.38/0.23. The 448k-particle density along the coordinate is one Gaussian blob.

**Why:** the per-particle 2D latent has SNR far below 1, so a per-particle embedding of this data cannot show modes
and its bins are diluted; the state structure lives at the sub-class (pooled) level and in maps.

**How to apply:** use the synchronised coordinate for pooling (bins, trajectories, as a prior for sub-class
assignment), not as a per-particle classifier; keep tracked sub-class states as the fine-band partition; a
cross-half test of the coordinate needs half expansions that dump particle latents (the clshalf runs did not).
Related: [[rnasep-alternation-verdict-2026-10-14]], [[rnasep-basis-span-finding-2026-10-13]],
[[rnasep-2d-subclasses-pooled-vs-particle-3d-latent-2026-10-13]].
