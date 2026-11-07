---
name: box600-simple-master-needs-64g
description: "TRAP: box-600 masters peak ~44-50 GB (request 64-96G); refine3D WORKERS peak ~100 GiB at a 4 A matching limit vs a 48 GB job_memory_per_task, so 3 per 256 GB node OOM — set job_memory_per_task=200000"
metadata:
  node_type: memory
  type: reference
  originSessionId: ea91e751-00ee-44f7-a73f-361aba072d8d
  modified: 2026-11-04T02:15:26.577Z
---

On batch.ncifcrf.gov the IgM 3D tests (box 600, 1.09 Å/px, nparts=10, nthr=16) ran their SIMPLE masters with `--mem=16G`. Peak RSS: T1 solve3D 40.6 GB (survived), T2 solve3D_addon 43.6 GB (killed), T3 solve3D nstates=3 43.6 GB (killed). The 16 GB request is not a hard cap, so jobs overshoot until the node runs short and the OOM killer fires, typically in the final full-box reconstruction.

**Why:** 2026-11-03. T3 died at stage 5 after 7 h. The T2 add-on died after its stage-8 iterations, before writing `solve3D_addon_report.txt`, the addon_diag map or the published project. A killed add-on leaves its working copy with frozen rows at state=0 and no frozen poses ([[solve3D-addon-working-copy-masks-frozen-rows]] if written).

**How to apply:** request 64G (`#SBATCH --mem=64G`) for any SIMPLE master at box ~600; scale from `sacct -o maxrss` of a comparable run ([[size-job-memory-to-measured-rss]], [[size-batch2-requests-to-the-job]]).

**Workers (2026-11-04):** refine3D refine=neigh, nspace=20000, box_crop=432, ~28-34k particles per part. Workers ran at about 100 GiB RSS once the matching limit reached 3.99 A, against `job_memory_per_task=48000`. Slurm does not hold them to that request; it packs 3 per 256 GB node, and the node OOM-killed A1 (part 7, exit 137, cn140). Set it with `update_project job_memory_per_task=200000`, one worker per node. Raising job_memory_per_task also enlarges the particle/record buffer, whose budget is half of it (`ptcl_buffer_budget`, src/main/strategies/search/simple_matcher_ptcl_io.f90): records at box_crop take box_crop^2*4 bytes per particle. So never pick a value that allows two workers per node and also switches the buffer on.
