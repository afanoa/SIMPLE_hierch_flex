---
name: batch-login-node-when-batch2-down
description: "when batch2.ncifcrf.gov is down or refusing logins, ssh batch.ncifcrf.gov (fsitgl-head01p) reaches the same slurm and beegfs for checks/sbatch; never compile there"
metadata:
  node_type: memory
  type: reference
  originSessionId: b12feb5b-ebbe-44c9-bfe5-fac8fee8014b
  modified: 2026-10-24T04:19:20.980Z
---

2026-10-24 00:10: batch2 (alias fsitgl-xfer03p) refused logins ("System is going down"); the user: "you can ssh to
batch.ncifcrf.gov but you cant compile code there". `ssh -o BatchMode=yes batch.ncifcrf.gov` lands on fsitgl-head01p with
the same sacct/squeue and /mnt/beegfs paths; running jobs (cn124) were unaffected.
**How to apply:** fall back to batch.ncifcrf.gov for job checks, file reads and sbatch when batch2 is unreachable;
no builds (make/cmake) there. Coordinator relay of the user's rule (2026-10-24): batch.ncifcrf.gov is THE login for all
cluster work from now on; never compile or run anything on the head node; builds go through a SLURM job. The login-node rules still apply ([[head-node-runs-nothing-but-sbatch]],
[[no-srun-or-ssh-into-compute-nodes]]). Related: [[batch2-ssh-needs-no-vpn]].
REVERSE CASE (2026-10-27 ~19:00): batch.ncifcrf.gov (head01p) itself became unreachable (ssh timeout; module loads on batch2
print a metrics-DB connect error to batch.ncifcrf.gov:3306, harmless). batch2.ncifcrf.gov kept working (squeue/sbatch, beegfs).
Watchers that poll squeue through the dead host read "no jobs" and fire falsely -> re-arm every watcher on the live host and
make watchers fail closed on ssh errors.
