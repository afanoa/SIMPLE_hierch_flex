---
name: site-dictionary-global-basis-verdict-2026-11-07
description: VERDICT NO (2026-11-07): localised (varimax/ICA) modes of the seam-free tree-wide basis cannot replace site discovery/ranking; single-mode closed AUROC qtoy2 .69-.74 vs span .92, r131 .59 vs site factors .65-.69; the closed-arm site's mode is the arm swing
metadata:
  type: project
---

Tested 2026-11-07 (f2r_sitedict.py, L/sitedict*/, handoff/f2r-global/VERDICT.md "Site dictionary"): the q = 30 seam-free basis
(every charted node's flex modes, input poses) rotated for spatial localisation inside the sphere (varimax; orthogonal FastICA),
z_rot = z R exact. Modes become compact (half the energy in 0.2-0.4 % of the voxels), site-resident (0.6-0.8 of a mode inside one
tree site mask) and reproducible across trees (qtoy2 .95-.99; r131 closed mode .64-.82 varimax). On qtoy2 the closed-arm mask's
mode is the closed-best mode, but one localised mode reads closed at AUROC .69-.74 against the 30-mode LDA .92, and a one-bit-per-
site dictionary label reaches NMI .27-.34 (z_MAP) against the tree's leaves .75-.78. On r131 no localised mode passes .60 (span
.63-.65 z_ML / .70 z_MAP; per-node site factors .65-.69; root chart .62) and the mode living in the closed-arm site does not
discriminate closed particles (.50-.53): the site's dominant mode is the arm swing, the open/closed occupancy is spread over
several modes. A z_MAP single axis is not a local readout (the fitted prior pools the span into every axis).

**Why:** a site's state needs several directions inside the site (what site factors and focused flex fit per node); a global
coordinate's localised axes each hold a slice of it, so site discovery and ranking cannot become a lookup.
**How to apply:** do not retry site dictionaries / localised global modes as a site layer substitute; keep per-node site discovery
+ site factors / focused fits. Related: [[global-chart-coordinate-verdict-2026-11-06]], [[flex-to-refinement-coupling-verdicts-2026-11-06]].
