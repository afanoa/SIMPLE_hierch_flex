---
name: tar-upload-keeps-old-mtimes-make-skips
description: "TRAP: tar-piping worktree files to batch2 (tar xf) restores their LOCAL mtimes; if the build dir holds objects built later from other versions, make compiles 0 files and the snapshot silently carries the old code. Check the compile count (0 = wrong)."
metadata:
  node_type: memory
  type: feedback
  originSessionId: 58de2c8f-ad53-44fb-b8a2-05aadff16ec4
  modified: 2026-10-23T22:17:45.475Z
---

2026-10-23, m8c37:
- I uploaded the approved flex files with `tar cf - files | ssh ... tar xf -`.
- Builds 35/36 had meanwhile compiled the COMMITTED versions of those files on batch2, so their objects were newer than my local worktree copies.
- tar restored the older local mtimes. make built 0 Fortran files, and the gate passed on the OLD code.
- I caught it from "Building Fortran" = 0 and a `strings bin/simple_exec | grep <new string>` = 0, moved the snapshot aside, touched the files and rebuilt (23 files).

**Why:** a snapshot that doesn't contain the change looks green: the gate passes and the commit goes in, but every later run uses old code.

**How to apply:**
- Extract with `tar xmf -` (don't restore mtimes), or touch the uploaded files after extraction.
- After every build, check `grep -c "Building Fortran" makeN.log` is > 0 when sources changed.
- For important changes, grep the binary for a string only the new code has before snapshotting or committing.

Related: [[batch2-full-rebuild-trap-and-montage-first]], [[never-push-before-build-and-gate]].
