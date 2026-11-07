---
name: integrin-10345-not-a-tree-target-2026-11-03
description: "aVb8 integrin EMPIAR-10345 processed end to end; user verdict: flex/tree is not the right tool for it (continuous arm motion + compositional heterogeneity)"
metadata:
  node_type: memory
  type: project
  originSessionId: f001000f-0802-4c45-a7e2-0dffbd7e0651
  modified: 2026-11-03T05:10:29.124Z
---

2026-11-03: EMPIAR-10345 (84,266 ptcls, box 300, 1.345 A/px, deposited CTF) on batch beegfs EMPIAR_downloads/10345, project
/mnt/beegfs/elmlund/afan/integrin10345/int10345 (refine3D_auto c1 4.39 A; tree_mu_s1 7 leaves; 6_refine3D_states 8.0 A avg;
10_flex_pca arm-masked 4 states). User: "I dont think flex tree is the right tool for this." RECOVAR (Gilles & Singer, Fig. 4)
finds residual compositional heterogeneity first, then a continuous wide-range arm motion with no clusters.

**Why:** the motion is a continuum on thin legs; the tree's splits were small local sites and the states looked undramatic.
**How to apply:** do not pick 10345 again as a tree/flex showcase; use it only as a continuous-method (RECOVAR/cryoDRGN/3DFlex) reference.
TRAP seen here: flex_pca with lp beyond the working Nyquist (box_crop=64 on a 403 A box) writes past the PCG wrap table and segfaults;
size box_crop so lp sits inside Nyquist. Related: [[tree-continuum-fixture-and-verdict-2026-10-29]], [[tmem16a-empiar10124-target-2026-11-01]].
