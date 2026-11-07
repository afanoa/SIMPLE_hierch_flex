---
name: empiar13553-and-not-abinitio2d-runs
description: "2026-09-29 abinitio2D runs for EMPIAR-13553 (1.34M ptcls, single 550 GB stack) and the unpublished ~8M \"not\" set; working dirs, parameters, job ids"
metadata:
  node_type: memory
  type: project
  originSessionId: b8bc4431-3743-4b17-aafd-6d31e354ea19
  modified: 2026-09-29T15:04:36.567Z
---

Submitted 2026-09-29 on batch2 csbdevel; user said never remove anything from these dirs.

- **EMPIAR-13553**: `/mnt/beegfs/elmlund/afan/empiar13553/` (first job 53430891 FAILED; retry 53431026 on `1_import_particles/empiar13553.simple`). Source
  `/mnt/beegfs/elmlund/EMPIAR_downloads/13553/data/particle_stack_neg.mrc` + `deftab.txt` (dfx/dfy µm, angast deg),
  1,341,961 ptcls, box 320 @ 1.04, cs 0.01 kv 300 fraca 0.07, imported as-is (no extra inversion, per user).
  `import_particles stk=` references the stack by abspath, it does NOT copy it. abinitio2D ncls=300 mskdiam=0
  (tight boxes; mskdiam=0 falls back to the max mask) nparts=10 nthr=8, 64 GB/worker.
- **not (unpublished; the user said ~8M but the 8_extract project holds 804,598 ptcls, 19,602 stks, box 420 @ 1.3; one strongly preferred view)**: `/mnt/beegfs/elmlund/afan/not/` (job not_ab2D, 53430892), a COPY of
  `/mnt/beegfs/elmlund/processed_test_datasets/not/8_extract/not.simple` (reboulcf's; never write it).
  The user's exact line: ncls=300 nparts=8 nthr=18 mskdiam=400 (deliberately large) job_memory_per_task=64000.

Both run as hand-written master sbatch scripts (update_project compenv -> slurm, then abinitio2D in-job), which avoids
the script=yes header patches in [[simple-script-yes-slurm-submission]]. The user asked me to watch the account's core total
([[batch2-shared-account-cpu-cap]]): about 530 once both sets of workers are up.

**Trap (cost one failed job):** `new_project dir=X` followed by `import_particles projfile=X/p.simple` (default mkdir=yes)
writes the POPULATED project to `X/1_import_particles/p.simple`; the root p.simple stays an empty stub, and a downstream run
on the root dies with "no particles found" in exec_calc_pspec_distr. Point downstream at the numbered-dir project, or pass
mkdir=no to import_particles.
