---
name: not-sampling-fix-reruns-2026-10-04
description: not 2026-10-04: abinitio3D on the user's 2D selection with latest master (sampling fix) vs flex2D-seeded control; dirs, builds, jobs
metadata:
  type: project
---

User 2026-10-04: "the sampling issue seems to be fixed" -> rerun the user's abinitio3D (their class selection,
not/3_selection/not.simple; original cmd `abinitio3D nthr=16 nparts=10 mskdiam=400 pgrp=c1`) on the LATEST code, and as a
control run the newest flex2D seeded from abinitio2D, then abinitio3D from that.
Root `/mnt/beegfs/elmlund/afan/not_20261004/`:
- `ab3d_usersel/` (run.sh): latest master built in a new worktree `~/SIMPLE_m1004` (origin/master 3a4d57819, Release ARCHOPT,
  build_m1004.sh, build/build.done); submitted automatically after the build.
- `flex2D/` (run.sh, job not_flex2D 53649978): snapshot `~/SIMPLE_flex2D/snap/d1a20eff6` (simple_exec + bin/ workers),
  copy of not/1_abinitio2D/not.simple (300 cls, 804,598 ptcls), `flex2D cls_init=prev ncls=300 nstates=4 mskdiam=400`,
  env FLEX2D_SHORTLIST=8 FLEX2D_POLISH_FREEZE=1 (recipe from impsamp/is_flex.sh), nparts=10 nthr=8.
- `ab3d_flex2D/` (to do when flex2D finishes): abinitio3D on the flex2D project, same args as ab3d_usersel, latest master
  (user: MUST use the latest binary; if master cannot read the flex2D project, stop and report, no fallback). No manual selection: flex2D deselects junk itself.
Related: [[impsamp-toy-campaign]], [[not-ribosome-flexpcg-plan-2026-10-02]], [[judge-2d-partitions-by-heterogeneity]].

**Release-4 traps hit (2026-10-04):** (1) abinitio3D was RENAMED solve3D (also abinitio2D->solve2D, cluster2D->refine2D,
5b8ad9f06, Oct 3): `prg=abinitio3D` -> "Program not recognized" (ab3d_usersel/ is that dead attempt). (2) the sampling fix
(f9e1f0ada view-balanced sampling) is OPT-IN: solve3D sets partition=no unless given; user chose `partition=yes`
(nclust=20 default; log "VIEW PARTITION: clustering 126 class averages into 20 groups"). (3) pre-Release-4 projects die in
the sigma2 bootstrap with "nptcls_stk not present in stack field; map_ptcl_ind2stk_ind" -> run `prg=fix_projfile` (in place,
reads stack headers) on a COPY first (ab3d_usersel_r2/ is that dead attempt). Working run: `ab3d_usersel_r3/run.sh`
(cp -> fix_projfile -> solve3D partition=yes mskdiam=400 pgrp=c1 nparts=10 nthr=16), job 53650117, past iteration 1.
The flex2D control project (snapshot build, old format) will need fix_projfile before its solve3D too.
RESULTS: user-selection solve3D partition=yes (job 53650117, 6 h, done 2026-10-05 01:09): FSC0.143 3.71 A / FSC0.5 5.15 A,
cFAR 0.07 (very anisotropic; the Sep-30 abinitio3D on the same selection was 3.16 A, cFAR 0.21). flex2D control (53649978)
done 22:44 (3 h 40, RC=0). Control solve3D = `ab3d_flex2D/run.sh` (cp flex2D/not.simple -> fix_projfile -> solve3D
partition=yes, same binary), job 53672158, started 10:13 on 2026-10-05.

**Matching Hans's run (his email, 2026-10-05):** same algorithm as ours (cc + average linkage, 20 groups) but his inputs
differ: NEW 2D (solve2D) with 285 selected classes / 431,798 ptcls, largest group 80 classes 42.8 %, and 25,005 sampled per
iteration (5.8 %, ~1,323/group, likely nsample=25000) vs ours: old Sep-29 abinitio2D, user's 126 classes, default 10,000/iter
(500/group) -> 3.71 A cFAR 0.07 vs his 3.1 A. Step 1 running: `not_20261004/solve2D/run.sh` (cp 8_extract -> fix_projfile ->
solve2D ncls=300 mskdiam=400 nparts=10 nthr=8, latest master), job not_solve2D 53672432 (started 10:31). solve2D has NO auto
class selection: the USER selects (asked to be pinged when solve2D is done; session cron d7b10568 pings). Step 3 after the
selection: solve3D partition=yes nsample=25000 mskdiam=400 pgrp=c1 on the latest binary. The flex2D-control solve3D
(53672158, default 10k) kept running.
