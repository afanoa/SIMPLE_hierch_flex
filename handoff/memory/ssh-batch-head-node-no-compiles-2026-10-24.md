---
name: ssh-batch-head-node-no-compiles-2026-10-24
description: "USER RULE 2026-10-24 ~00:20 (via coordinator): ssh to batch.ncifcrf.gov (head node fsitgl-head01p, same beegfs + SLURM as batch2); NEVER compile on that head node - builds go through a SLURM job"
metadata:
  node_type: memory
  type: feedback
  originSessionId: 851a7b77-9a45-48b1-b2b6-396c1a7b5829
  modified: 2026-10-24T04:19:22.632Z
---

The batch2 login node rebooted around 00:00-00:15 on 2026-10-24. The user then made `ssh -o BatchMode=yes batch.ncifcrf.gov`
the way in: same /mnt/beegfs, same SLURM queues, jobs started from batch2 keep running.

**Why:** batch2's login node went down; the batch head node is shared and must not carry compile load.

**How to apply:** use batch.ncifcrf.gov for sbatch/squeue/sacct and file reads; any `make`/cmake build runs inside an sbatch
job (this overrides the "builds may run on the login node" allowance in [[build-on-batch2-login-node]] for this head node).
The banner noise goes to stderr (filter with 2>/dev/null). Related: [[head-node-runs-nothing-but-sbatch]].
