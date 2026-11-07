---
name: rnasep-24-refine3d-auto-quarter-unrefined
description: "RNase P 24_refine3D_auto project (rnasep_states/afm_cryoem.simple = rnasep_regions copy, G10's input): 240,104 of 932,496 particles (25.7 %) never refined (e3 = x = y = 0, e1/e2 copied in batches); they made fake islands in 2D sub-classes and the G10 3D latent; use flags_refined.txt (692,392)"
metadata:
  node_type: memory
  type: project
  originSessionId: 34b66a02-03d2-42df-91f8-915afbf80120
  modified: 2026-10-13T15:14:45.851Z
---

Found 2026-10-13. In `rnasep_states/afm_cryoem.simple` (byte-identical to `rnasep_regions/afm_cryoem.simple`; the
user's "24_refine3D_auto, best poses"; map `24_refine3D_auto/rec_final_state01.mrc`, FSC 0.143 ~3.9 A) a quarter of
the particles carry an unrefined pose: in-plane angle and shifts exactly 0, projection direction shared by groups
(7,912 distinct directions). Most likely refine3D_auto's fixed per-iteration sample never reached them. They are spread
evenly over the refine3D_states states (25-28 % each), so they are not the states 3/5 junk sinks, and refine3D_states
(25_*) and S1/S2/S4 (rnasep_hires, refined from 25 finals) did refine them.

Effect when the 24 poses are used directly: in class 350 the 3,379 unrefined members formed a separate island in BOTH
the per-particle 2D latent (pose-mean residual is wrong at a wrong pose) and the G10 3D latent; across five classes
the sub-class/view NMI fell from 0.08-0.21 to 0.008-0.026 and the 3D latent's view encoding from 0.36-0.94 to
0.08-0.13 once they were excluded.

**How to apply:** for any flex/coupled work on the 24 poses, deselect them first:
`/mnt/beegfs/elmlund/afan/rnasep_coupled/flags_refined.txt` (932,496 rows, 1 = refined, 692,392) with
`selection oritype=ptcl3D infile=...`; or refine them before using the full set. Check `e3 = x = y = 0` counts on any
project that went through refine3D_auto with fractional sampling. Related: [[rnasep-stacks-are-class-prealigned]],
[[rnasep-2d-subclasses-pooled-vs-particle-3d-latent-2026-10-13]].
