---
name: beegfs-transient-write-failure-2026-10-27
description: "TRAP: a ~6 min beegfs write failure (23:30-23:36 2026-10-27) killed trees cluster-wide: 'failed to write byte buffer' (simple_fileio 1115), FSYNC ERROR, 'file_size cannot be inquired' (binoris 103); not a code bug"
metadata:
  node_type: memory
  type: project
  originSessionId: 9643ce3d-fc7a-4dc1-860a-59db10867e7c
  modified: 2026-10-28T03:43:34.709Z
---

2026-10-27 ~23:30-23:36 EDT: /mnt/beegfs writes failed for a few minutes (no quota or space issue: 63 % used, quotas unlimited;
/home at 92 %). Symptoms in tree.log / step logs: `ERROR! failed to write byte buffer; simple_fileio.f90 line 1115` (simple_copy_file),
`FSYNC ERROR FILE UNIT: -10`, `file_size cannot be inquired; simple_binoris.f90 line 103`, and `step X FAILED (exit code 1), node
stopped`. A truncated 2,097,152-byte copy (one MAXBUFSZ chunk) is the signature. Killed: OLE recipe1 tree (0 leaves), qtoynull
recipe1, full692 route-A arm (8 nodes), olenull_sst, ole_fixm, data-lp olenull_on; the full692 base arm lost 2 nodes.
`dd ... conv=fsync` on beegfs worked again at 23:37.

UPDATE 23:45 (data-lp session): not only transient. `beegfs-ctl --listtargets --nodetype=storage --spaceinfo` showed targets 6202 and
6203 (bgfs06p) at 0.0 GiB / 0 inodes free (37 of 165 under 50 GiB). Files are striped over 4 consecutive target IDs, full ones
included; any file landing on 6202/6203 fails ("Remote I/O error" on append, "failed to write byte buffer", empty PCG artifact).
A rewritten file KEEPS its targets, so a resume fails again on the same file. Check with `beegfs-ctl --getentryinfo FILE`; move
the step dir aside so new files get fresh targets.

**Why:** an error right after a fallback-seeding log line looked like a code bug; the same error in 5 unrelated trees within
minutes showed it was storage.
**How to apply:** when a tree dies on an I/O error, grep every live tree.log for the same message before debugging code; test
with a small `dd conv=fsync`; move dead run dirs aside (mv, never rm) and relaunch, or resume in place by setting the failed
steps aside ([[never-delete-only-read-or-move]]). Related: [[batch-login-node-when-batch2-down]] (head01p was down the same evening).
