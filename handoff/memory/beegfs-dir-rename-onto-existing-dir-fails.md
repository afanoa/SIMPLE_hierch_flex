---
name: beegfs-dir-rename-onto-existing-dir-fails
description: "TRAP: on batch2's BeeGFS, rename(2) of a folder onto an existing folder fails with EEXIST even when the target is EMPTY (POSIX/ext4 allow it); simple_rename only prints 'IOS# 17' and returns, so check the result"
metadata:
  node_type: memory
  type: feedback
  originSessionId: 6f195b76-225b-466c-a432-7fc8360148e8
  modified: 2026-11-08T23:33:20.990Z
---

2026-11-08, sparse-sites band ladder `eb1d99ede`: the node made an empty `sites/` at start. To take over the cap job's
results it then called `simple_rename(sites_cap/sites, sites)`. BeeGFS refused with errno 17 (EEXIST), confirmed with
`mv -T` on empty scratch folders. A rename to a fresh name one level up works.
- `simple_rename` (src/fileio/simple_syslib.f90) only prints "simple_rename failed ... ERROR: IOS# 17" and returns, so the run
  carried on and reported "0 site blocks in 0 sites" on IgM 10 %. The message's file name ends in a NUL, which makes
  `grep` treat flex.log as binary and stop printing; read such logs with `grep -a`.
- Fixed in `e4fbef54e`: the folder is created only where it is written, and the takeover happens only when no target
  exists, then is checked by `file_exists` (fallback: compute it in-process).

**Why:** a folder move that works on a laptop silently loses results on the cluster's filesystem, and box-size tests
that never reach the move path don't catch it ([[sparse-sites-band-ladder-verdict-2026-11-08]]).

**How to apply:** in SIMPLE code, never rename a folder onto one that may exist (even empty). Create output folders
lazily, test `dir_exists` before the move, and verify the move with `file_exists` on a file inside. Validate any new
takeover/move path on a fixture that actually reaches it. Related: [[simple-chdir-failure-does-not-abort]].
