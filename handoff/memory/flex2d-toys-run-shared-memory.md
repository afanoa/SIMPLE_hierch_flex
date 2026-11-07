---
name: flex2d-toys-run-shared-memory
description: SUPERSEDED by use-the-full-cpu-cap-for-tests - flex2D toy runs are Slurm-latency bound when distributed; shared memory was faster but the user wants distributed 80-100-core arms
metadata:
  type: feedback
---
2026-09-30: 50k-particle hemoglobin toy runs with nparts=10 took 20-40 min while using ~40 cores of 1,933 idle: every
flex2D iteration submits 10 worker jobs that each run for seconds, so the run waits on Slurm scheduling. The user
("I should have enough cpus available") was right that CPUs were not the limit. Shared memory (NPARTS=1 NTHR=32,
sbatch -c 32, synth2D/scripts/dev_flex_smp.sh): 7 iterations in 72 s, a full run in ~10 min.

**Why:** per-iteration distributed dispatch has fixed latency far above the compute of a toy iteration.

**How to apply:** dev/toy runs (<= ~100k small particles) go single-node shared memory; keep nparts=10 x 8 threads for
real data sets. Related: [[simple-run-parallelism-preference]], [[flex2d-findings-2026-09-21]].

2026-09-30 follow-up: on real hemoglobin (nparts=10) each iteration's workers waited 21-28 s in the queue for ~30 s
of compute. The user declined SIMPLE's persistent-worker backend ("too complicated"): do not propose it; use one
shared-memory job on a large csbdevel node (nparts=1, many threads) when dispatch latency dominates.


SUPERSEDED 2026-09-30 evening: the user's standing rule is [[use-the-full-cpu-cap-for-tests]] (every test arm distributed,
nparts=10-12 x nthr=8, several in parallel) and [[keep-test-loops-minimal]] (no sbatch dependencies). Follow those;
the per-iteration Slurm wait (~20-30 s) is a known cost, not a reason to go shared-memory.
