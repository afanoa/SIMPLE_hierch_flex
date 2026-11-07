---
name: head-node-runs-nothing-but-sbatch
description: on batch2 the login node runs ONLY sbatch/squeue/file reads and builds; every SIMPLE call, scoring script, dump and montage goes inside an sbatch job
metadata:
  type: feedback
---
2026-09-22, the user (angry, "I have told you many many times"): nothing of ours may run on the batch2
head node except job submission, queue inspection, reading files, and compiling (which they allowed).
That includes the LIGHT things I had been doing there: `print_project_field` dumps, the Python
scoring/centroid/montage scripts, and long `sleep` wait loops inside ssh sessions.
**Why:** the head node is shared and policed; long or CPU-using processes there get killed and reflect on
the user. The rule is not about job size.
**How to apply:** put dumps, comparisons and montages at the END of the job script that produced the data
(the flex2D_tests run.sh already do this), or submit a tiny sbatch for analysis. To wait for results,
poll with short ssh calls (squeue + cat), never a sleep loop that lingers on the login node. Related:
[[flex-tests-run-on-cluster-only]], [[batch2-archopt-sigill-login-node]], [[announce-commands-before-running]].
