---
name: ctffind-4114-segfaults-on-cluster
description: "TRAP: cluster module ctffind/4.1.14 segfaults silently under RELION (empty micrographs_ctf.star); use ctffind/4 (default) or 5.0.2"
metadata:
  node_type: memory
  type: reference
  originSessionId: 851a7b77-9a45-48b1-b2b6-396c1a7b5829
  modified: 2026-11-02T23:52:07.075Z
---

On the NCI cluster, `/mnt/nasapps/production/ctffind/4.1.14/bin/ctffind` segfaults right after reading its inputs (0.732 Å RNase P micrographs, 2026-11-02). relion_run_ctffind_mpi does not stop: it writes a micrographs_ctf.star with no micrograph rows, and extraction then fails with "no rlnMicrographName column". This cost an 8-hour stage-1 job.

`/mnt/nasapps/production/ctffind/4/bin/ctffind` (module `ctffind/4`, default) and `ctffind/5.0.2` both work. Pass the binary explicitly with `--ctffind_exe`. Test one micrograph before a full run.
