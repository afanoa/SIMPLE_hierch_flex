---
name: box-oom-from-paired-pcg-solve
description: "2026-09-23: my paired even/odd flex PCG state solve (nested omp sections under the 80-thread master boost) took fwl-c142339 down (OOM, box unreachable); cap threads like production and measure RSS before any concurrent-solve test"
metadata:
  node_type: memory
  type: feedback
  originSessionId: 846f1d9d-8f6f-4d4a-b48d-b023b38c86e4
  modified: 2026-09-23T19:26:44.800Z
---

The flex_pca PCG state change that solves even/odd in `!$omp parallel sections` (2×nthr/2) was
built and run on the box without a thread cap; the machine went unreachable and the user
reported an OOM. Production's rec3D PCG master caps itself at PCG_MASTER_NTHR_CAP=32 for the
paired solve; flex's master runs at the boosted nthr_glob (80), and the solver's inner regions
size per-thread scratch from nthr_glob.

**Why:** the box is shared (elmlundho had jobs running); a memory blow-up there is not a
private failure.

**How to apply:** never test a new concurrency/memory change at scale first. Cap the threads
(production's 32), read the startup RAM report, measure RSS on a small case, and only then run
20k. Check dmesg for the OOM record and other users' processes after any such incident.
Related: [[box-shared-with-elmlundho]], [[flex-pca-staterec-memory-ceiling]].
