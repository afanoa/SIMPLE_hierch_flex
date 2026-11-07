---
name: flex-latent-density-basins-no-power-2026-10-22
description: "Basins (ToMATo, kNN density) of the flex_pca_coordinates.txt latent have NO power: swing30 4 kT double well, 5 % separate state and simsites discrete states all read 1 basin; TRAP: those z are XD posterior means, raw MAP z + stored precision sit in flex_pca_embedding.bin (layout inside)"
metadata:
  node_type: memory
  type: project
  originSessionId: f0f8499a-d0a6-450b-9f93-483f98f27006
  modified: 2026-10-22T05:34:38.995Z
---

2026-10-22 overnight, for the 2D-to-3D coordinating session (data-derived width, no path rule): `routes/toy55/land_latbasins.py`
(kNN log-density, k 100/300, D 10 or 3, ToMATo persistence on the symmetric kNN graph, half check).
- **Truth:** swing30 s30a node_r22 latents (subsampled by true t): sim / flat / broad / dw2 / dw4 / shoulder ALL give 1 basin; top
  prominences never exceed ~0.4 kT (noise level), so dw4 (3.86 kT) and the 5 % state are invisible. simsites focus latent (8 discrete
  site combinations, arm closed 8 %): 1 basin, ARI 0 with truth.
- **f692st r1, three runs (script, port1, port2):** 1 basin each. Agreement is trivial; it says nothing about states.
- **TRAP:** `flex_pca_coordinates.txt` z are POSTERIOR MEANS under flex's extreme-deconvolution mixture prior (they equal the trailing
  block of `flex_pca_embedding.bin`), not the raw MAP latents. Point-estimate densities are noise-blurred and shrunk; never count states
  from them.
- **flex_pca_embedding.bin layout (version 5):** magic 8 B, version i4, nptcls i4, ncomp i4, box i4, smpd r4, pinds i4[n],
  z f8[n,ncomp] (Fortran order), eigvals f8[ncomp], contrast / resid / resid_mean f8[n] each, precision f8[ncomp,ncomp,n], sig2_eff f8;
  then 'SIMPLFXD', noise_scale f8, flag i4, n i4, ncomp i4, deconvolved z, precisions, labels i4[n]. The size check is exact.
- **The stored precision is the estimator's sampling precision Q = A Gt^+ A** (`map_sampling_precision`, `simple_flex_pca_embed.f90`
  :247/:377), but `simple_flex_pca_deconv.f90` (`noise_and_projection`, `calibrate_noise_scale`) treats it as the posterior precision A
  (D = precision - P). For particle-components dominated by the prior, this reads a shrunk value as a precise one. Flag it to the user
  before relying on flex's deconvolved prior.
**How to apply:** for a data-derived state count use a deconvolved density (per-particle likelihoods: [[path-model-pose-marginalised-reweighting-2026-10-22]],
or XD through the per-particle noise), with truth calibration, never the latent point cloud.
