---
name: flex-nparts-spawns-slurm-workers
description: "flex_pca with nparts>1 on a project whose compenv is slurm submits every part as its own sbatch job (simple_<projname>, cpus=nthr, 125 GB) -- 3 'jobs' at nparts=10 became 17 jobs / 544 cores on 2026-09-24, breaking the user's per-job and job-count limits; use nparts=1 when a core/job cap is given"
metadata:
  node_type: memory
  type: feedback
  originSessionId: 2590d10f-7aff-4aff-b2cc-d9d85f962ad3
  modified: 2026-09-24T16:35:26.243Z
---

2026-09-24, cnga1 region arms (`flex_region/v2_*`, region binary, nparts=10 nthr=32 inside a 32-cpu sbatch): each paired
probe/embedding round submitted 10 worker jobs named `simple_cnga1_test` (32 cpus, 125 GB each) through the project's
compenv; 120 worker jobs total, **peak 14 concurrent workers = 448 cores** on top of the 3 masters. The user had asked for
<= 80 cores per job and <= 3 jobs at a time. The `Submitted batch job` lines are in FLEX_PCA_OUTPUT (and cyl45_n10's
single-node 53-min / 36 GB run had looked in-process, which is why I missed it).

**Why:** the user polices batch2 usage (shared account CPU cap, see [[batch2-shared-account-cpu-cap]]); silent fan-out
breaks their explicit limits and crowds out other sessions.

**How to apply:** when the user sets a core or job-count cap, run flex_pca with `nparts=1` (in-process, nthr = the job's
cpus), or check `grep -c "Submitted batch job" FLEX_PCA_OUTPUT` a few minutes after launch; count worker jobs in
`squeue -u atifao` (name `simple_<projname>`) before claiming a job count. Related:
[[simple-script-yes-compenv-trap]], [[simple-run-parallelism-preference]].
