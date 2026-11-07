---
name: batch-unlimited-partition-shares-norm-nodes
description: "batch.ncifcrf.gov unlimited partition = 30 of norm's nodes with its own 432-CPU per-user cap; move pending workers with partition=unlimited qos=normal when the norm cap is full"
metadata:
  node_type: memory
  type: reference
  originSessionId: ea91e751-00ee-44f7-a73f-361aba072d8d
  modified: 2026-11-03T19:04:50.291Z
---

On batch.ncifcrf.gov the `unlimited` partition (cn084-109, cn164-167) is a subset of `norm`'s nodes: the same hardware, not a separate pool. Its per-user cap is separate: partition QOS `infinite`, cpu=432, against norm's QOS `norm`, cpu=600. Workers submitted by SIMPLE carry QOS `norm` explicitly. Moving them needs `scontrol update jobid=<id> partition=unlimited qos=normal`; changing only the partition leaves them under the norm cap.

**Why:** 2026-11-03. The IgM 3D test workers were stuck behind `QOSMaxCpuPerUserLimit` because the RNase P jobs used the 600-CPU norm cap, while 25 norm nodes sat idle. The user first doubted it was a different partition (it isn't, hardware-wise), then approved the move ("ok thats fine you can move to unlimited").

**How to apply:** ask before using it for a new campaign, and say it shares norm's nodes. Move only my own workers, by job ID, selected by WorkDir. SIMPLE re-submits each iteration's workers to the project's `qsys_partition`, so a run needs a mover loop until its masters end. Related: [[batch2-shared-account-cpu-cap]], [[scancel-only-by-job-id-never-by-name-pattern]].
