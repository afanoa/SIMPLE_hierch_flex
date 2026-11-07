---
name: reconstruct3d-negative-e1-degrades-map
description: plain reconstruct3D on the cryoSPARC-imported 10028 project (52% of e1 negative) gives 4.5/8.3 A; the SAME rotations re-gauged to e1,e3 in [0,360) give 3.3/4.0 A; always re-gauge or refine before a control
metadata:
  type: project
---

Found 2026-09-30 on `/mnt/beegfs/elmlund/afan/10028/10028_simple/RIBO` (105,247 ptcls, box 360 @ 1.34 A, poses from
the cryoSPARC import via W_fam/p.simple). `reconstruct3D pgrp=c1 mskdiam=320 sigma_est=global nparts=10` on the raw
project: e/o 4.5 A (0.143) / 8.3 A (0.5), mobile half 6.2/9.5. The identical project with e1/e3 rewritten through
euler2m -> m2euler (numpy, `~/toy3/patch_project.py`, e2/x/y bit-identical, mean |d e1| 188 deg = the gauge change
only): 3.3/4.0 A, 3.4/4.0. So a negative e1 (55,036 of 105,247 rows) is mishandled somewhere in the reconstruction path
(the toy projects, written by simulate_particles, never had negative angles; refine3D re-gauges, which is why the
W_prop eqrec maps were fine at 3.77 A).

**How to apply:** never use a plain reconstruct3D of an imported-pose project as a control; re-gauge the Eulers first
(patch_project with CAL=0 does exactly that) or reconstruct after one refine3D iteration. Worth a code look: where
does reconstruct3D / the orientation class assume e1 in [0,360)? Related: [[simple-star-pose-import-traps]],
[[toy3-body-frame-test-system]].
