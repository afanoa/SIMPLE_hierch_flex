---
name: keep-test-loops-minimal
description: "the 2026-09-30 cls_split test loop became \"a big mess\" through harness elaboration; scripts live on batch2 only (never staged through the Mac scratchpad, which vanishes at compaction), no sbatch dependencies, check qsys_name before distributed runs, one submit + one poll"
metadata:
  node_type: memory
  type: feedback
  originSessionId: 03e8b343-4213-4208-890c-56346e7de4b2
  modified: 2026-09-30T20:30:00.000Z
---

2026-09-30, cls_split pca_mode=flex on 10180. What went wrong, in order: sbatch --dependency on a
prep job whose last `ls` failed (arms held forever); arms on a project whose compenv is
qsys_name=local (workers ran on a 4-core node); after a context compaction, a chain that wrote its
job script into the Mac scratchpad died on `cd` (the scratchpad is wiped) and its waiter polled a
file that could never complete, so an hour passed with nothing running. The user: "absolutely
ridiculous ... you make this big mess" and "why is this taking so long, our test loop should be way
faster to actually develop".

**Why:** every staging layer (Mac scratch files, scp, dependency chains, waiters on waiters) is a
place to fail silently; each silent failure costs the user's waiting time, and the real tests are
fast (the flex unit suite runs in 2 s; a distributed 46k split in 2 min).

**How to apply:**
- Job scripts live on batch2 (heredoc over ssh); never stage through the Mac scratchpad.
- One sbatch, one poll on the job's own output; never --dependency, never a waiter on a waiter.
- Before any distributed run on a copied project: `update_project ... qsys_name=slurm`.
- cls_split worker stdout is NOT kept (qsys_cleanup); per-class diagnostic prints are only
  visible in a shared-memory `class=N` run or must be written to a part file.
- Give the user the one-line command that shows the result; scorers and montages after.
Related: [[flex2d-deliver-results-first]], [[use-the-full-cpu-cap-for-tests]],
[[head-node-runs-nothing-but-sbatch]].
- batch2 banner filters like `grep -v "^  "` also delete squeue's right-aligned rows, so a
  "while squeue shows the job" waiter exits at once; poll `squeue -h -o '%T'` and match the state.

2026-10-20 17:25 (user: "764 cores is a lot to use. And we aren't even testing fast"): I had put four 400k synthetic sets
through the FULL f692e chain (parent board, family refinement, region boards, zooms, consolidation; hours each, hundreds of
worker cores) to answer method questions. Fast loop instead (routes/toy55/fast_factor.sbatch): <= 100-120k particles, no
refinement (simulated poses are near-true), one node (-c 64, local boards, no worker jobs), whole-particle board -> regions
-> focused boards -> factor_sample.py + truth, ~20-30 min. Also count only MY jobs when reporting cores: other sessions' runs
(e.g. f692k board workers, 320 cores) share the account.
