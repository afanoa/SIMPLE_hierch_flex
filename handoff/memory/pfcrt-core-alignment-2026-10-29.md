---
name: pfcrt-core-alignment-2026-10-29
description: "PfCRT 13_solve3D particles (56,809) aligned to deposited EMD-20806 inside its automask: greedy global pass + neigh x10, PCG with pcg_mskfile; /mnt/beegfs/elmlund/afan/pfcrt_core, job 54626045"
metadata:
  type: project
---

2026-10-29, the user: align the selected PfCRT particles to the deposited model with a solvent mask, to align on the core.
- Dir `/mnt/beegfs/elmlund/afan/pfcrt_core` (README, PfCRT.simple copy of 13_solve3D, ref/, scripts/core.sbatch, 1_/2_refine3D).
- Reference `PfCRT_smpd_1p3/picking_map.mrc` = EMD-20806 at 1.3 A box 240; mask = default automask of it (lp 8, 1 dilation, 6 vox
  edge): 0.50 M A^3, extent 139 x 108 x 101 A (vs 0.12 M A^3 for PfCRT + Fab alone: generous or residual micelle).
- refine3D keys: pgrp=c1 mskdiam=200 lp=8 box_crop=128 rec_backend=pcg pcg_mskfile=ref/automask3D.mrc nparts=10 nthr=8;
  step 1 refine=greedy maxits=1 vol1=masked deposited map (global search sidesteps docking + hand); step 2 neigh maxits=10.
- RERUN job 54627628 with a 6UKJ MODEL mask (atoms + 3.9 A, 5-vox cosine edge, 0.19 M A^3; user: "looks good", wants micelle +
  solvent removed, not per-helix). Density masks of EMD-20806 (default automask, thres, 8/5 A at 2-6 sd) all keep the micelle ring:
  the deposited map has a dense micelle. 6UKJ coords -> our box: voxel = (xyz + 110*1.035)/1.3 - 87, array [x,y,z] order.
  SIMPLE automask thres= did not tighten it. TRAP: stale vol1_even/odd next to vol1 get reused -> clear them.
- Facts: pcg_mskfile only constrains the PCG reconstruction (refused only by refine3D_states_cavgs); refine3D_auto has no mask key.
Related: [[pfcrt-tree-run-2026-10-27]], [[tree-user-mask-opt-in-not-default]], [[refine3d-shc-virgin-field-trap]].
- 3_refine3D = FULL-box greedy pass (user: box_crop 128 was "ridiculous", it capped maps at 4.9 A Nyquist; full box 5.8 A, 7 min).
- Flex: 4_flex_pca (job 54634235, recipe17b build, tree root-flex keys + pcg_mskfile=6UKJ mask, 10x8).
- Tree: recipe17b_umask (= recipe17b + tree_mskfile.diff: opt-in key clips half_variance's sphere, so every focus/site/step mask, and
  sets the root flex pcg_mskfile; gate 15/15, key in simple_args) -> pfcrt_core/tree_umask_s1 (job 54634739; recipe17b gate keys,
  nparts 5, minpop 3000, nodepar 4, kway_maxdepth 2). Diff: hetero_tree_20261028/recipe17b_umask/tree_mskfile.diff.
- tree_umask_s1 leaf maps = micelle blobs (user: "doesn't look like pfcrt"; Downloads copies deleted). Cause: pcg_focus consensus
  background solved over the mskdiam SPHERE (micelle) + unmasked node means/plain leaf maps; masked leaf 01 vs greedy corr 0.58 (same hand).
- Fix recipe17b_umask2 (= umask + mean_leaf_mask.diff): new key pcg_bgmskfile (pcg strategy background support, default sphere);
  tree user_mask_keys(): reconstruct3D/refine3D without a mask get pcg_mskfile=tree_mskfile, pcg_focus ones get pcg_bgmskfile.
  Gate 15/15. Rerun pfcrt_core/tree_umask2_s1 (job 54656774); verify a leaf vs the greedy map before delivering.
