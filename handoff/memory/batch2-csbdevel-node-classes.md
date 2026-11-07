---
name: batch2-csbdevel-node-classes
description: "csbdevel on batch2 has two hardware classes - cn001-007 are 96-core/2TB, the rest 48-core/251GB; pick cn001-007 with --constraint=x6740P for high-neigs flex_pca"
metadata:
  node_type: memory
  type: reference
---

`batch2.ncifcrf.gov` (login lands on fsitgl-xfer03p, 36 cores, Skylake - hence the ARCHOPT SIGILL trap).
Partition `csbdevel`: 32 nodes, 1872 CPUs, 5-day default AND max walltime, no GPUs. Two classes:

- **cn001-007** (7 nodes): 96 cores, RealMemory 2,064,022 MB (~2.0 TB), 1.4 TB local TmpDisk, feature `x6740P`
- **cn116-132, cn162-169** (25 nodes): 48 cores, ~251 GB, ~875-895 GB TmpDisk, feature `x6342`

Select the big ones with `#SBATCH --constraint=x6740P`. QoS caps a single user at **642 CPUs**, which
binds before the partition total. cn169 was drained for maintenance 2026-09-14.
`largemem` also exists (2 nodes, 96 cores, 3.1 TB, infinite time) but the user asked for csbdevel only.

Queue is usually empty on csbdevel; beegfs has ~4.4 PB free and **no quotas** (beegfs-ctl --getquota
returns empty tables). /mnt/beegfs is native BeeGFS on batch2 but an NFSv3 re-export on fwl-c142339.

**Why it matters:** flex_pca per-process RAM is ~2 x npairs-scaled kernel blocks (the PAIRED engine holds
both fits), so cnga1 at box_crop=128 needs 250 GB at neigs=32 - it OOMs on a 48-core node and fits on
cn001-007. See [[flex-pca-staterec-memory-ceiling]], [[batch2-archopt-sigill-login-node]].

2026-10-21 00:00: a `-c 64` job can only land on cn001-007 (the 96-core class); when those are busy it pends with
"ReqNodeNotAvail, UnavailableNodes:..." although ~25 48-core nodes sit idle. `scontrol update NumCPUs=48` is not enough
(MinCPUsNode stays 64); resubmit with `sbatch -c 48` and it starts at once (64 threads on 48 cores is fine for SIMPLE).
