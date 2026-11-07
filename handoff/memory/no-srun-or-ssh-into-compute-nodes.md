---
name: no-srun-or-ssh-into-compute-nodes
description: "batch2: never run anything on a compute node directly (srun --jobid --overlap into a running job, ssh cnXXX); everything goes through sbatch; inspect running jobs only via logs, squeue, sacct, sstat"
metadata:
  node_type: memory
  type: feedback
  originSessionId: b12feb5b-ebbe-44c9-bfe5-fac8fee8014b
  modified: 2026-10-24T02:36:53.493Z
---

2026-10-23 (OLE tree, job 54181360 on cn124): I ran `srun --jobid=54181360 --overlap ... ps` to see the tree's processes
on its node. User: "you cant operate directly on a node you have to submit via slurm".

**Why:** compute nodes are shared and scheduled; attaching extra steps to a running allocation (or ssh to a node) bypasses
the submission rules the user set for batch2.

**How to apply:** on batch2 only sbatch submits work; to see what a running job does, read its output file and step dirs
on beegfs, or use `squeue`, `sacct`, `sstat -j <id>` (CPU/memory per step) from the login node. No `srun` into an
existing job, no `ssh <node>`. A diagnostic that needs a compute node is its own small sbatch job, announced first.
Related: [[head-node-runs-nothing-but-sbatch]], [[announce-commands-before-running]].
