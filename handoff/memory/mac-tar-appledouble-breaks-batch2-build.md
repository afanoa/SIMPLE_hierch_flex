---
name: mac-tar-appledouble-breaks-batch2-build
description: "TRAP: `tar cf - files | ssh batch2 tar xf -` from the Mac writes ._<file> AppleDouble twins; SIMPLE's build then fails (default_audit.py UnicodeDecodeError 0xa3; cmake GLOB also compiles ._*.f90) -> use COPYFILE_DISABLE=1 tar"
metadata:
  type: reference
---

2026-10-21 (state-tree port, ~/SIMPLE_stree build 10): syncing 25 changed files with a tar pipe from the Mac left 25
`._<name>` files next to them on batch2. `make install` died at the UI-defaults step: scripts/default_audit.py reads every
exec/commander/stream source as UTF-8 and the AppleDouble header byte 0xa3 raised UnicodeDecodeError; the configure-time
`GLOB_RECURSE main/*.f90` would also have picked up `._*.f90` as sources.

**Why:** macOS tar stores extended attributes / resource forks as `._` entries unless told not to.

**How to apply:** EVERY tar sync from the Mac starts with `COPYFILE_DISABLE=1` (recurred 2026-10-23, ~/SIMPLE_m8c build 17:
same 0xa3 failure, one wasted build): `COPYFILE_DISABLE=1 tar cf - <files> | ssh batch2.ncifcrf.gov 'cd <tree> && tar xf -'`
(or scp per file); if `._*` files appear, list them (`find src -name '._*'`), mv exactly those aside (no rm on batch2,
[[never-delete-only-read-or-move]]), rerun `cmake .` and make.
Related: [[batch2-new-source-needs-cmake-rerun]], [[build-on-batch2-login-node]].
