---
name: pfcrt-tree-run-2026-10-27
description: "PfCRT hetero_tree test (seed 1, m8c90t, OLE-scale keys) from dibyajyoti's 13_solve3D (56,809 ptcls) in /mnt/beegfs/elmlund/afan/pfcrt_tree, controller job 54398962"
metadata:
  node_type: memory
  type: project
  originSessionId: b2360bb0-fdc3-4e91-8fc3-27dd729625f7
  modified: 2026-10-27T19:58:22.371Z
---

2026-10-27, the user: "try the tree flex with this molecule" from dibyajyoti's selection + solve3D; picked **13_solve3D** (from
11_selection = the selected 10_solve2D classes; 1 state, mskdiam 200, 8.2 A) over 18_/20_/21_solve3D (16_selection lineage).
- Run dir `/mnt/beegfs/elmlund/afan/pfcrt_tree`: README.md, PfCRT.simple (copy + update_project slurm/csbdevel/32 GB/user_email=nil),
  scripts/tree.sbatch (resumable), HETERO_TREE_OUTPUT, 1_hetero_tree. Root 56,809 particles, 36 Bonferroni split tests, z_crit 3.01.
- Keys: mskdiam=200 pgrp=c1 minpop=3000 maxdepth=20 nparts=5 nthr=8 nodepar=3 site_layer=no oneflex=yes sta0_npeaks=2; junk_screen
  left at the code default (the "no" rule is RNase P only). Snapshot ~/SIMPLE_stree/snap/m8c90t.
- KILLED on the user's ask at 1 h 33 min (controller + 5 step jobs, by ID), during the site-7 test: model/masks not good enough yet.
- Root: sites 1-6 no split (site 1's 5,009-ptcl pole 1.27/0.14 SE). No step masks solvent or micelle: mean/flex/recB0 = 200 A
  sphere only; fflex + stA0 use node_r/mask.mrc (variance focus mask) = 2.92 M A^3 = 70 % of the sphere, only 13 % of it above
  1 sd of the mean density, and it misses ~20 % of the dense molecule. A protein-only mask would need EMD-20806 / a model map.
- dibyajyoti is actively running solve3D in that tree (21_, 22_ on 2026-10-27): never use a dir still running.
- TRAP: with a projfile, qsys_env takes the partition from the project's compenv, not SIMPLE_QSYS_PARTITION -> update_project the copy.
- TRAP: never `grep -r` across beegfs run dirs over ssh: it ran >2 min on the login node; killed it by PID.

Related: [[pfcrt-clsexp-run-2026-10-27]], [[tree-launch-fast-defaults-2026-10-25]], [[site-layer-off-for-test-runs-2026-10-27]].
