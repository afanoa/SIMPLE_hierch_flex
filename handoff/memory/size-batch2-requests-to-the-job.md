---
name: size-batch2-requests-to-the-job
description: "batch2 csbdevel: size --mem/job_memory_per_task/walltime to the actual job, never copy NOT's 128 GB / 1-day values; most free-memory nodes are often DRAINing, so big requests pend for hours"
metadata:
  node_type: memory
  type: feedback
  originSessionId: 230e38d3-a42d-4cd7-ab95-1ddcf5c084f4
  modified: 2026-10-06T18:26:18.403Z
---

2026-10-06, thyroglobulin bootstrap_rec3D (33,742 ptcls, box 448): I copied 128 GB/task and a 1-day
limit from the NOT scripts. The job sat PENDING (Priority, est. start +4.5 h) because the csbdevel nodes
with unallocated memory were all DRAIN; only ~60 GB slots were open. Then I cancelled and resubmitted
three times, and one resubmit died in 1 s because my `grep "NORMAL STOP" update_project.log` guard
was wrong (update_project prints no NORMAL STOP line). The user: "this is stupid".

**Why:** every resubmit cycle is user waiting time; the requests must be right on the first sbatch.

**How to apply:** before the first sbatch, check open slots (scontrol show node: CPUAlloc/AllocMem/State,
skip DRAIN) and size to the job: a box-448 distributed reconstruction needs ~10-20 GB per part,
so `job_memory_per_task=40000`, master `--mem=48G`, `-t 04:00:00` and `walltime=14400`, which can
backfill. Guard steps on `$?`, not on log strings. Related: [[keep-test-loops-minimal]],
[[batch2-shared-account-cpu-cap]], [[thyroglobulin-challenge-dataset]].
