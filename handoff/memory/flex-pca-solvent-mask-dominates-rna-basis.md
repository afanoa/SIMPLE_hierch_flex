---
name: flex-pca-solvent-mask-dominates-rna-basis
description: "On the afm RNA data pc1/pc2 were solvent-shell noise because mskdiam=210 left a 70-105 A solvent annulus 3x the molecule's volume; the band+concentration eigenvolume diagnostic caught it and mskdiam is the lever"
metadata: 
  node_type: memory
  type: project
  originSessionId: 126275fa-8460-444b-8ba7-1e0140968435
  modified: 2026-08-13T17:42:04.962Z
---

Found 2026-08-13 on run `31_flex_pca` (932k RNA particles, clean equal-mass sextile readout, yet
frames showed only core "breathing"). Radial analysis of the eigenvolumes told the story:

- consensus mass: 98.7% inside 70 A radius (molecule ~140 A diameter)
- pc001: 51% of |pc|^2 energy at 55-105 A, running flat to the 105 A mask edge and stopping there;
  top-1% concentration of that peripheral energy = **0.02** (uniform = 0.01) -> structureless haze
- core modes pc004/pc008: concentration 0.24-0.28 -> proper lobes (the control)
- verdict: pc1/pc2 = solvent-annulus noise variance. mskdiam=210 (radius 105) on a ~70 A-radius
  molecule leaves ~3x the molecular volume of unmasked solvent, which is the largest-variance
  subspace and takes the top of the spectrum. Whitening cannot remove it (spatial, not spectral);
  the image mask (COV_MASK_IMAGES via mskdiam -> cov_image_mask_radius) is the only lever --
  the RECOVAR mask_images_in_H_B equivalent. The run-31 trajectory was therefore placed along a
  solvent coordinate -> frames = consensus + jitter ("breathing"), periphery silent. Retroactively
  explains run 29's latent geometry and the SNR-selection OUTER-half warning too.

**Diagnostic worth keeping (and automating):** radial band fractions + top-1% concentration of
each eigenvolume's energy, computed on the box-64 pc*.mrc in milliseconds -- solvent modes are
peripheral+diffuse (conc ~uniform), motion modes are peripheral+concentrated, core modes are
interior. Fits the [[flex-pca-must-be-knob-free]] constraint (cheap, latent/volume-tiny, could
auto-rank components or auto-flag a bad mask at runtime).

Run `32_flex_pca` = run 31 with mskdiam=170 (85 A radius, headroom for domain excursions past the
consensus envelope). Check after: pc1 peripheral concentration should rise toward lobed values or
the spectrum re-ranks; if flexible domains extend past ~85 A the user's 2D classes know the true
extent -- adjust mskdiam accordingly, runs cost ~23 min. See
[[flex-pca-gmm-undoes-equal-mass-placement]], [[afm-rna-flex-pca-needs-continuum-readout]].
