---
name: batch2-memory-not-enforced
description: "batch2 slurm does not kill jobs over --mem (ConstrainRAMSpace=no, NoOverMemoryKill); an under-request oversubscribes a node instead of failing"
metadata:
  node_type: memory
  type: reference
  originSessionId: d63181a2-43e0-44aa-a87c-6f7c55606053
  modified: 2026-11-08T07:15:02.663Z
---

Checked 2026-11-08:
- `/etc/slurm/cgroup.conf` has `ConstrainRAMSpace=no`, and `scontrol show config` gives `JobAcctGatherParams = NoOverMemoryKill`.
- `SelectTypeParameters` includes `CR_CORE_MEMORY`, so `--mem` still drives scheduling.
- A job that uses more than it asked is NOT killed; it silently oversubscribes the node's memory.

Example: hetero_tree asks 8 GB + 0.025 MB per particle for the `site_split` worker (10.5 GB at 100k). On IgM with blocks, the
worker measured 13.1-13.7 GB MaxRSS and still completed.

**How to apply:** a job surviving is no proof that its request is right. Size requests from sacct MaxRSS
([[size-job-memory-to-measured-rss]]). When a node OOMs or swaps for no visible reason, look for co-located jobs that
under-requested. Related: [[size-batch2-requests-to-the-job]].
