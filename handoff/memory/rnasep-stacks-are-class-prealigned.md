---
name: rnasep-stacks-are-class-prealigned
description: "RNase P 932k batch2 projects (rnasep_states, rnasep_cons/all, rand284k) point at stacks already rotated/shifted into the 4_abinitio2D class frame; ptcl2D e3=x=y=0 is CORRECT; re-applying the original alignment double-rotates (discs)"
metadata:
  node_type: memory
  type: project
  originSessionId: 34b66a02-03d2-42df-91f8-915afbf80120
  modified: 2026-10-13T14:25:04.752Z
---

Measured 2026-10-13 (`/mnt/beegfs/elmlund/afan/rnasep_coupled/prep2D`): the batch2 RNase P projects read the 10
consolidated phase-flipped stacks `/mnt/beegfs/elmlund/rnasep/22_denoise_project/denoise_raw_particles_partNN.mrcs`
(ctf=flip, box 160, 1.464 A). Their ptcl2D keeps the 4_abinitio2D class labels (500 classes, 112 populated, run in
`afan/atifao/reextract/4_abinitio2D`, final `cavgs_iter044.mrc` box 160) but e3=x=y=0 for all 932,496 rows, and that
is right: make_cavgs on rand284k as-is matches the original averages at ncc 0.995 (40-10 A). Transferring e3/x/y from
`reextract/6_selection` (same particles, same order: pind/xpos/ypos/defocus/class identical; stkind/indstk differ
because of the restacking) gives rotational discs. The box-320 `cavgs_*_partN.mrc` files there are padded partial
sums, not a different box.

**How to apply:** use the 4_abinitio2D classes directly from these projects (no fresh solve2D, no alignment restore);
cls_expansion needs nothing else from 2D. `make_cavgs` in the flexref5 build needs `nparts` (it prints the error and
exits 0). Related: [[rnasep-2d-subclasses-pooled-vs-particle-3d-latent-2026-10-13]].
