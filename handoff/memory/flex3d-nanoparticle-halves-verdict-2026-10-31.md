---
name: flex3d-nanoparticle-halves-verdict-2026-10-31
description: "flex_pca on Pt fission halves (NP87, 71-17 p1/p3) adds nothing beyond the occupancy model; traps: 3-shell default band, deconvolved coordinates, pose-tied half agreement"
metadata:
  node_type: memory
  type: project
  originSessionId: 0254fd76-6e74-4609-9ade-f6bfe0782cb1
  modified: 2026-11-01T03:30:47.088Z
---

2026-10-31 (doc/implementation_notes/grfp_tseries_flex3d_2026_10_31.md, cluster pilot/flex3d, scripts fx3d_*.py):
flex_pca (master 6ddc4a0d2) fitted per raw fission half at lp 1.5 / 5 / 16 A finds nothing the ramp occupancy
model misses.
- At 1.5 A the axis gate drops 55/60 axes, and cross-half cos is <= 0.47.
- Site changes follow the fades weakly and never beat the ramp model against the published segment maps.
- The only off-site half-consistent change (NP87 at 16 A, outer shell) looks like shared background.

**Why:** settles the "does flex see something the site model can't" question for the GRFP pilot; do not rerun it
without independent poses.

**How to apply / traps:**
- Default lp=16 on a 57 A box = a 3-shell band (kstop=3), with per-frame Gram cond ~10^3.
- flex_pca_coordinates.txt holds the DECONVOLVED latents; they go rank one when one axis carries the signal. The MAP
  latents are only in the raw block of flex_pca_embedding.bin (v5 stream layout).
- This build's embedding holds contrast at 1 (contrast file all 1.000); the leak measured R² <= 0.01 here.
- Fission halves share the poses, so latent half agreement can be view-tied misfit (p3: 0.79, view-predicted
  0.4-0.5). Always test against a same-view predictor. Related: [[nanox-pt-trajectories-and-grfp-exp]],
  [[flex-pca-deconv-precision-mismatch-flag-2026-10-22]].
