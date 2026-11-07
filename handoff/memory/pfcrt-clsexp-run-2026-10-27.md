---
name: pfcrt-clsexp-run-2026-10-27
description: "PfCRT cls_expansion on dibyajyoti's 10_solve2D (100 cls, ~85k) runs in /mnt/beegfs/elmlund/afan/pfcrt_clsexp, ncls=5000, 40 cores, job 54393327, build ~/SIMPLE_proc aea5b1846"
metadata:
  node_type: memory
  type: project
  originSessionId: b2360bb0-fdc3-4e91-8fc3-27dd729625f7
  modified: 2026-10-27T18:06:53.261Z
---

2026-10-27, the user: cls_expansion on PfCRT, "maybe 50 classes per subclass", in my own directory, never touch dibyajyoti's tree, 40 cores.
- Source (read only): `/mnt/beegfs/elmlund/dibyajyoti/SIMPLE_Projects/PfCRT_smpd_1p3/10_solve2D/PfCRT.simple` (solve2D 100 classes,
  mskdiam 200, box 240 at 1.3 A; all 100 classes, not the 11_selection subset). Its stack/cavg paths are absolute, so a plain cp works.
- Run dir `/mnt/beegfs/elmlund/afan/pfcrt_clsexp`: README.md, PfCRT.simple (copy + update_project slurm/csbdevel/32 GB/user_email=nil),
  scripts/clsexp.sbatch, 1_cls_expansion. `ncls=5000 mskdiam=200 nparts=4 nthr=8` + 8-core master; job 54393327 (54393283 aborted at
  my own path check). DONE 13:50, 270 s, NORMAL STOP: 3374 subclasses (minpop 25 caps a class at pop/25), 20 groups, none kept whole.
- User judged step 1 "pretty bad"; step 2 = same + 2D solvent prior (~/SIMPLE_stage1 498ff04e7, subcls_avg=kernel) in
  2_cls_expansion, job 54394842, 281 s, 3374 subclasses; assignments/FRCs differ from step 1 (prior active). Verdict pending.
- TRAP: grepping `stk=` in print_project_field output also matches `nptcls_stk=`; anchor on whitespace.

Related: [[pfcrt-poses-are-the-blocker]], [[persist-project-lineage-no-new-projects]], [[cls-expansion-stage1-ready-2026-10-12]].
