---
name: ptcl-oris-drop-unregistered-keys
description: "SIMPLE particle oris (ptcl2D/ptcl3D rows, is_ptcl) have a zero-capacity hash - o%set('newkey', v) is silently dropped; only the enumerated keys of simple_defs_ori survive; and a distributed master hangs forever when its slurm workers die"
metadata:
  node_type: memory
  type: reference
  originSessionId: 34b66a02-03d2-42df-91f8-915afbf80120
  modified: 2026-10-13T18:58:20.685Z
---

`new_ori(is_ptcl=.true.)` sets `htab = hash(0)` (`src/main/ori/simple_ori.f90:142`), so a key outside the
`simple_defs_ori` enumerators (e.g. `pop`) set on a ptcl3D row disappears without an error. cls2D/cls3D rows
(non-ptcl) do keep extra keys. Found 2026-10-13 when sigma_pop read `pop` from cavgs work-project ptcl3D rows:
every worker threw in 3 s, and the distributed master (refine3D_states_cavgs) sat in "ITERATION 1" for 30 min
because it does not notice FAILED slurm workers - check `sacct` of the worker job IDs in the run's out file.

**How to apply:** carry per-row metadata through a cls2D/cls3D field indexed by the row's `class`, or a registered
key; after any distributed submit, a run silent past ~5 min means check worker sacct, not wait. Related:
[[rnasep-cavg-states-findings-2026-10-13]].
