---
name: use-the-full-cpu-cap-for-tests
description: "test jobs on batch2 are distributed at 80-100 cores EACH (nparts=10-12 x nthr=8), several in parallel to fill the 642-core cap; never a shared-memory-only path for testing"
metadata:
  node_type: memory
  type: feedback
  originSessionId: 03e8b343-4213-4208-890c-56346e7de4b2
  modified: 2026-09-30T18:30:56.340Z
---

2026-09-30, testing cls_split pca_mode=flex. The user objected twice: first "why is testing so slow,
should be distributed with a lot of cores" (I had serial 16-thread arms), then, when I proposed
640-core single arms, "this is stupid... I would limit each job to 80-100 cores (of course
distributed). We cannot develop if you make these sorts of decisions."

**Why:** the cap is shared across the user's own jobs; one fat job blocks them, one thin serial job
wastes hours. A shared-memory-only code path (my pooled fit) is a defect for development because
it cannot be tested at speed.
**How to apply:** every test arm = distributed master with nparts=10-12, nthr=8 (80-100 cores);
launch several arms at once rather than one big one; check `squeue -u atifao -o "%C"` and keep the
running total under ~642 including the user's jobs. Any new path must have a distributed form
before it is tested on real data. Related: [[flex-nparts-spawns-slurm-workers]],
[[batch2-shared-account-cpu-cap]].
