---
name: not-sampling-fix-reruns-2026-10-05
description: not 2026-10-05 Hans-matching reruns: new solve2D (300 cls) + user's 151-class selection -> solve3D nsample=25000 balance=cavg on master 311b8817a; flex2D control on same build; dirs/jobs/selection facts
metadata:
  type: project
---
Root `/mnt/beegfs/elmlund/afan/not_20261004/`. Builds: `~/SIMPLE_m1004` (3a4d57819), `~/SIMPLE_m1005` (origin/master 311b8817a).
- `solve2D/` job 53672432: solve2D ncls=300 mskdiam=400 on all 804,598 ptcls (m1004), RC=0 14:26. Final = cavgs_iter031.mrc
  (box 420, smpd 1.3); last stage ends with a terminal greedy all-particle pass + full-size make_cavgs (~40 min after "iter 30").
- User selected on the iter-25 preview (artifact https://claude.ai/artifact/Jhw153zT1eYnnmJXuNuvRC, click-to-reject grid,
  scratchpad page.py/mont.py pure-python since the Mac python has no numpy). Selection carried by index: iter25 vs iter30
  self-corr min 0.953 median 0.999. Kept 151 classes / 508,321 ptcls (63%); rejected ranks 131-300 (median 222).
- `ab3d_solve2Dsel/` job 53682283: `selection oritype=cls2D deselfile=rejected_classes.txt` -> solve3D nsample=25000
  mskdiam=400 pgrp=c1 (balance=cavg default) on m1005, started 14:28. Single run (user: "I dont need two replicates").
- `ab3d_flex2D_cavg/` job 53675638: control, flex2D project (all 300 cls, NO selection) same solve3D on m1005. balance=cavg
  grouped 234/300 classes (78% of ptcls) into group 1 -> group 1 gets 5% of each draw; singleton outlier groups oversampled
  (one 1066-ptcl class sampled every iteration). Grouping the cavgs takes ~9 min before iteration 1.
Hans's reference run: 285 classes, 431,798 ptcls, 25,005/iter, 20 groups, group 1 42.8%, 3.1 A.
Related: [[not-sampling-fix-reruns-2026-10-04]], [[chain-steps-inside-one-job]], [[map-cavgs-selection-ranked-stack-trap]].

RESULTS (all RC=0, finished 2026-10-05 evening; FSC0.143 / FSC0.5 / cFAR):
- ab3d_solve2Dsel (user solve2D sel, 151 cls / 508k): 3.29 / 4.44 / 0.22  (largest cavg group 23%, visits 1.02-75.6)
- ab3d_flex2Dsel_v2 (user flex2D sel, 79 rejected -> 221 cls / 609k): 3.37 / 4.55 / 0.31 (first 42-class list was
  superseded by the user's corrected list; v1 dir ab3d_flex2Dsel stopped via STOP file before solve3D)
- ab3d_flex2D_cavg (flex2D, NO selection, 300 cls): 6.35 / 9.93 / 0.05 -> unselected input + balance=cavg = collapse.
Compare Oct-4 old-code run on old selection 3.71/5.15/0.07; Hans 3.1 A.
