---
name: numba-parallel-before-fork-deadlocks
description: "TRAP (batch2 Python, 2026-11-07): numba.set_num_threads() or any parallel kernel in the parent before a multiprocessing fork deadlocks the workers; set NUMBA_NUM_THREADS and run parallel kernels only after the fork"
metadata:
  node_type: memory
  type: project
  originSessionId: 1218bb2b-a21a-4d25-9e48-c40ddd89ceff
  modified: 2026-11-08T02:02:29.630Z
---

In `tree_latent/scripts/tl_tree5.py`, a module-level `numba.set_num_threads()` started numba's thread pool in the parent process.
The `multiprocessing` fork pool's workers then hung on their first `prange` kernel. The jobs used 15 CPU-seconds in 6 minutes.

**Why:** forked children inherit a thread pool whose threads do not exist in the child.

**How to apply:**
- Set the thread count with the `NUMBA_NUM_THREADS` environment variable in the job script.
- Create the fork pool before the parent calls any parallel kernel (the selfcheck ran after the fork).
- Check a new job's AveCPU with `sstat` after a few minutes: CPU time far below the wall time means a hang.

Related: [[tree-latent-verdict-2026-11-07]].
