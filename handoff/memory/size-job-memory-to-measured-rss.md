---
name: size-job-memory-to-measured-rss
description: "don't request 128 GB per task by habit; size job_memory_per_task / --mem from measured MaxRSS (sacct) with ~1.5-2x headroom; cls_expansion workers peak ~17-20 GB at 100 subclasses/class on NOT"
metadata:
  node_type: memory
  type: feedback
  originSessionId: b8bc4431-3743-4b17-aafd-6d31e354ea19
  modified: 2026-10-06T18:35:37.558Z
---

User 2026-10-06: "we should lower memory per task. I dont think we need 128gb memory for this at all."
Measured (sacct MaxRSS, NOT 452k ptcls, box 420): cls_expansion ncls=100 workers 17-20 GB; failed/empty workers 1.6 GB.

**Why:** 128 GB per task blocks scheduling on the shared csbdevel queue for no reason.

**How to apply:** before writing update_project job_memory_per_task or #SBATCH --mem, check sacct MaxRSS of a previous
run of the same program/size and request ~1.5-2x that (cls_expansion: 32 GB). If no measurement exists, say so and
pick a moderate value, then measure after the first run.
Related: [[batch2-shared-account-cpu-cap]].
