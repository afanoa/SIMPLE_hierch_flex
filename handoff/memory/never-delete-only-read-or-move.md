---
name: never-delete-only-read-or-move
description: "STANDING (user 2026-10-22): never delete anything on batch2 or the raw share - only read or move (mv aside); some data is delicate"
metadata:
  node_type: memory
  type: feedback
  originSessionId: 851a7b77-9a45-48b1-b2b6-396c1a7b5829
  modified: 2026-10-27T18:38:48.593Z
---

User 2026-10-22 (RNase P full-data planning): "We do not delete anything. We only read or move. Because I dont want
to delete something particularly delicate."

**Why:** the RNase P tree mixes delicate raw sessions, reprocessing rounds and other sessions' runs; a wrong rm is unrecoverable.

**How to apply:** no rm / rm -r / overwrite-in-place on cluster data; to retire a run, `mv` it into an `_old/` or dated folder
and say where it went. `tar xf` / `scp` uploads overwrite silently: give every new script a name not yet used in the target
dir and `ls` it first (2026-10-23: a tar upload replaced the Oct-22 stree_junkval.sbatch; restored from the session
transcript's Write call, the only copy: Slurm keeps no batch scripts here). Scripts that clean up (remove_chunks, tmp dirs) must be checked for deletes before use.
2026-10-27 (coordinator, after a stray `rm -f /dev/null` in an ad-hoc ssh command): NO rm of any kind on the cluster, not even
scratch or one-liners; grep scripts for rm/unlink/os.remove/rmtree before upload. Related:
[[never-overwrite-scripts-under-running-jobs]], [[say-upfront-when-i-cannot-do-the-ask]].
