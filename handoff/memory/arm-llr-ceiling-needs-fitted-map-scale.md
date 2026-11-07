---
name: arm-llr-ceiling-needs-fitted-map-scale
description: "TRAP: simulate_particles images are ~300x weaker than the map projections (LS scale 0.003); a likelihood ceiling without a fitted per-image map scale is pure map-power bias (one map takes 86 % of images, pooled 1000 still at chance). With the scale: 55 RNase P conformers d' 1.0 per image, IgM arms d' 1.5 at true poses"
metadata:
  type: project
---

arm_llr.py (toy55, 2026-10-19): scoring every conformer map against each image at its true pose with the map
amplitude as is gave chance accuracy and NEGATIVE pair d' (10th pct -2.6); llr_bias.py showed one map took 86 % of
the argmins. Cause: simulate_particles output is ~0.003x the projection amplitude, so |P|^2 dominates the SSE.
**Why:** the earlier "IgM arms are below the per-image SNR floor, only pooled evidence can see them" conclusion
came from this bug and is wrong. **How to apply:** always fit the per-image LS map scale (--lsscale 1, default
now) and print the most-chosen-map fraction; a ceiling where one map wins most images is a bias, not a result.
With the fit, true poses: 55 conformers d' 1.0 (acc 0.2/55, 10 pooled images 84 %), two_50k pair 0.85 per image,
IgM arms3 100 conformers acc 0.45, d' 1.5. Refined poses (angerr ~12 deg mean): 55 -> d' 0.46 (1000 pooled 54 %),
pair 0.70 per image. Related: [[klem-synthetic-arms-first-result-2026-10-17]], [[klem-null-leaves-sort-by-alignment-quality-2026-10-19]].
