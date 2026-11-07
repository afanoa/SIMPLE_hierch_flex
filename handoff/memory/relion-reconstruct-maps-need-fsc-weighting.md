---
name: relion-reconstruct-maps-need-fsc-weighting
description: "user 2026-11-05 'these relion states just look like noise': raw relion_reconstruct maps are unweighted to Nyquist; deliver relion_postprocess default output (FSC-weighted, B=0, unmasked postprocess.mrc) instead"
metadata:
  node_type: memory
  type: feedback
  originSessionId: a25e6f6d-c94c-4ba6-92f7-c1fb7553af83
  modified: 2026-11-05T06:25:16.637Z
---

2026-11-05, OLE per-state maps. I delivered Cpose/<L>/full.mrc (relion_reconstruct, all particles) as the "raw" map. These carried
96-99 % of their power beyond 6 A and had skewness ~0. The user: "these relioon states just look like noise". The RNase P maps the user
liked were run_class001.mrc from relion_refine, which RELION already weights by the FSC.

**Why:** relion_reconstruct has no FSC weighting or regularisation, so with 1.5-9k particles high-frequency noise swamps the molecule.
**How to apply:** for maps reconstructed outside relion_refine, deliver `relion_postprocess --i h1 --i2 h2 --mask m` with no
--auto_bfac/--adhoc_bfac/--low_pass: postprocess.mrc is then FSC-weighted (sqrt(2FSC/(1+FSC)), zero past FSC = 0), unsharpened and
unmasked (skew +1.7 to +4.6). That counts as "raw"; never a global-B map ([[deliver-maps-unsharpened-for-judging]]).
Check a map's skewness or high-frequency power before handing it over.
