---
name: flex-pca-deconv-precision-mismatch-flag-2026-10-22
description: "FLAG (unverified, from the motion session 2026-10-22): flex_pca's deconvolution treats latent%precision as the posterior precision A, but embed stores the MAP sampling precision Q = A Gt^+ A; they differ for poorly measured components"
metadata:
  node_type: memory
  type: project
  originSessionId: 25d47be4-7d1d-48d8-abdf-97b1c3c946ab
  modified: 2026-10-22T07:38:48.801Z
---

Reported by the motion/occupancy session overnight 2026-10-22. Not changed; it awaits the user, who wrote flex_pca.

**The claim:**
- Embed stores the MAP sampling precision Q = A Gt^+ A (`simple_flex_pca_embed.f90:247/377`).
- The deconvolution (`simple_flex_pca_deconv.f90` noise_and_projection / calibrate_noise_scale) treats it as the posterior precision A.
- For well-measured components they agree.
- For components a particle barely sees, the prior-shrunk value is read as a precise measurement.
- Fix: store A or Gt next to Q.

**Ruled out (2026-10-22):** unmodelled outside-mask variance is not behind the 0.41-0.47 noise scale: under pcg_falloff it moved
only to 0.44, and it cannot act through sig2_eff, which is the images' calibrated high-frequency noise floor, not a model residual
[[flex-pcg-falloff-prior-pulls-solvent-shell-2026-10-22]].

**Related trap:** `flex_pca_coordinates.txt` z are posterior means under flex's XD prior, not raw MAP latents, so a density over them is
shrunk.

**How to apply:** verify against the code before any deconvolution-based occupancy or basin result is trusted, and before changing the
code ask the user, since it is their method. Related: [[flex-pca-authorship-and-headline-result]].
