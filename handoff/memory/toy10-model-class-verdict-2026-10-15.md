---
name: toy10-model-class-verdict-2026-10-15
description: "Truth-map model-class test on the RNase P toy10 conformers (no particles): a smooth deformation field per conformer explains 63/62/51/40 % per band vs 35/36/36/33 % for a 5-component linear volume basis; deformation space is ~2-D (volume space 3.5-4.7-D); regularisation stiffness decides; composition (domain deleted) needs occupancy, added density needs an additive component; deformation and composition axes separate"
metadata:
  type: project
---

`model_class.py` on the 10 docked Degenhardt conformers vs the image-weighted consensus (box 160, 1.464 A, 200 A
mask), unexplained fraction of V_i - V_c per band 117-30 / 30-20 / 20-12 / 12-8 A (batch2 routes/toy55/model_class_s3,
_s6, _c3; 2026-10-15):
- linear leave-one-out span of the other conformers' differences: K=1 0.70/0.75/0.76/0.80, K=2 0.65/0.68/0.70/0.75,
  K=5 0.65/0.64/0.64/0.67 (K=all is exact only because the consensus is the weighted mean). Volume-difference PCA
  effective dimension 3.5 / 3.9 / 4.2 / 4.7 per band (PC1 share 0.47 -> 0.39): the landscape is NOT low-rank in
  volume space, which is why flex_pca's few linear components only "breathe".
- multi-scale symmetric demons deformation of the consensus (global smooth field, fluid sigma 5 px, elastic sigma 3 px):
  0.37/0.38/0.49/0.60 mean; best conformers 0.19-0.29 at low res (37, 27, 48, 49), worst = the two CLOSED conformers
  18 and 51 (0.60-0.68) = the largest rearrangement. Displacements p90 3-7 A, max 5-12 A. Flow-field PCA across the
  10 conformers: effective dimension 2.6, mode 1 = 57 %, modes 1+2 = 77 %: the motion is ~2-D in deformation space.
- stiffness decides: elastic sigma 6 px leaves 0.65/0.61/0.62/0.65 (displacements 1.4-3.2 A), sigma 8-14 px 0.7-0.9.
  The true displacement field is sharp at domain scale; a global smooth field must be weakly regularised (or
  piecewise), otherwise it collapses to breathing like the linear basis.
- composition (model_class_c3, variants from make_variants.py): conformer 49 minus the mobile domain: deform alone
  0.64/0.64/0.57/0.61, deform + multiplicative occupancy 0.08/0.22/0.45/0.58 (occupancy p5 0.06 finds the hole);
  conformer 46 plus a surface blob: 0.43 -> 0.41 (occupancy in [0,1] cannot ADD density: an added component needs
  its own additive/occupancy term). The deformation mode-1 score of each variant equals its parent's (+0.62 vs +0.65,
  -0.37 vs -0.37): deformation (conformational) and occupancy (compositional) axes separate cleanly.
- residual after deform+occ (0.34-0.59) = non-smooth/local-resolution differences between deposited maps.

**How to apply:** a latent model for RNase P must be nonlinear in volume space: deformation-type latent (1-2 modes)
plus per-component occupancies; a linear volume basis needs >= 4-5 components per band and still misses 2/3.
Verify any new latent against this table before particle runs. Related: [[toy10-latent-night-2026-10-15]],
[[rnasep-toy55-design-2026-10-15]], [[no-region-basis-for-rnasep]].
