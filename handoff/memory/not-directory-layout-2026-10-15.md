---
name: not-directory-layout-2026-10-15
description: "NOT /mnt/beegfs/elmlund/afan/not layout after the user's cleanup: steps sit side by side in not/ (1_abinitio2D, 3_selection, kernel100_expansion, 4_solve3D_cavgs), not nested in kernel100_expansion; everything else in _archive_2026-10-15"
metadata:
  type: project
---
2026-10-15. User: "1. abinitio2D 2. 3_selection 3. 100 cls expansion 4. solve3D_cavgs is all I want in the directory",
"it should not be done inside kernel100 expansion should be outside", "please run on multiple parts".
- not/1_abinitio2D, not/3_selection (126 classes), not/kernel100_expansion/1_cls_expansion (kernel, ncls=12600, master 5e1e6aeee)
- not/4_solve3D_cavgs: solve3D_cavgs on ALL 12,600 kernel averages (user: "dont want to do the selection"), run 36's
  line (b57538b7b, nstates=3, lpstart_ini3D=20 lpstop_ini3D=8, mskdiam=400) distributed nparts=10 nthr=8, job 53894463;
  run.sh + slurm log inside; mkdir=yes outdir=4_solve3D_cavgs
- 4_solve3D_cavgs result: 4434/211/7955 (state 2 emptied). User picked state 3 ->
  not/5_selection (selection oritype=cls3D infile=5_selection/state3.txt, 1 = cls3D state 3; selection state=3 would
  skip rows without updatecnt) -> not/6_refine3D_states_cavgs nstates=3 (defaults, nparts=10 nthr=8), job 53904793.
  User: no particle reconstruct3D; refine3D_states_cavgs seeds poses against the project consensus (still 4's state-1
  map) and rebuilds the flex mean from the selected averages itself.
- refine3D_states_cavgs needs a registered `vol` (state 1); solve3D_cavgs registers its maps only as `vol_cavg`
  (4_solve3D_cavgs: 3 vol_cavg, no vol) -> first try failed. User: "make it accept vol1". Patch (uncommitted, not
  pushed): refine3D_states_cavgs takes vol1= as the consensus (captured + deleted from cline before the state
  controller's vol1..N refusal; preflight box/smpd check) + UI entry + how2 line; lives in batch2 ~/SIMPLE_m1015
  (git diff there), built as snap/b57538b7b_vol1. 6_refine3D_states_cavgs rerun with vol1=4_solve3D_cavgs/
  rec_final_state03.mrc, nstates=3, nparts=10 nthr=16, job 53904819 -> OOM-killed after 18 min: its flex_pca
  initialisation runs in the MASTER, shared memory, at box_crop=420 (cavgs policy), >=160 GB on 15,910 rows. User: one
  fat node, 1 part, nthr=32, 500 GB -> job 53905208 on cn001 (2 TB node); refine3D_states_cavgs REQUIRES nparts on the
  cline (nparts=1), SIMPLE_QSYS=local.
- not/36_kernel_cavg3D_master kept as reference (its input stack is in the archive's 34 dir)
- not/_archive_2026-10-15: everything else (24, 34, 35, cls_expansion test dirs, root logs/scripts, cancelled chain
  attempts); nothing deleted, the user deletes it.
**How to apply:** new steps go in not/ with `mkdir=yes outdir=<N>_<prg>` (plain mkdir=yes would number past 36_ to 37);
next step runs on 4_solve3D_cavgs/proj.simple; the user decides the next step - do not chain extra steps unasked. See [[persist-project-lineage-no-new-projects]].
