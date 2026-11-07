---
name: keep-cluster-busy-queue-next-job-first
description: "user 2026-10-22 (via coordinator): never let batch2 sit idle between steps; queue the next job BEFORE analysing the last one, and use the full cap (~400 csbdevel cores, high nodepar)"
metadata:
  node_type: memory
  type: feedback
  originSessionId: 58de2c8f-ad53-44fb-b8a2-05aadff16ec4
  modified: 2026-10-22T17:52:24.722Z
---

2026-10-22 the user noticed the cluster idle for minutes after two analysis jobs finished while I wrote up results, and
asked (relayed by the coordinating session) to launch the next run immediately and "from here on keep the next job queued
before analysing the last one, so the cluster doesn't sit idle between steps"; use csbdevel up to about 400 cores, nodepar
as high as is safe.

**Why:** wall-clock on the shared batch2 cap is the bottleneck of the overnight campaigns; analysis can overlap with runs.
**How to apply:** when a job finishes, first submit the next planned job (or the obvious follow-up), then read and report the
finished one. Size distributed trees for the cap (e.g. hetero_tree nparts 6 x nthr 8, nodepar 4 per tree). A single sbatch
cannot exceed one node's cores on csbdevel (-c 128 is rejected; 96 works). Related: [[use-the-full-cpu-cap-for-tests]],
[[keep-test-loops-minimal]], [[batch2-400-core-limit-and-worker-timeouts]].
