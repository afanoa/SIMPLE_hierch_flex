---
name: batch2-full-rebuild-trap-and-montage-first
description: "every `make -j install` on batch2 recompiled all 603 objects (10 min) during the cls_split work; the user: 'why are you building from scratch every time', 'stop wasting my time I want to see montages'. Diagnose the cascade once, never chain blind rebuilds, and deliver montages before analysis"
metadata:
  node_type: memory
  type: feedback
  originSessionId: 03e8b343-4213-4208-890c-56346e7de4b2
  modified: 2026-10-02T00:30:00.000Z
---

2026-10-01 night. Chasing a rotation-window defect in cluster2D's inpl mode I did five
edit -> `git reset --hard` -> `make -j 24 install` -> run cycles; each build recompiled every
object (603, starting at simple_defs), ~10 min each, and all five runs gave bit-identical results,
so the loop measured nothing. The user, waiting for montages: "why are you building from scratch
every time. stop being stupid", "ridiculous, this is an easy fix but you cant figure out how to
compile code. just sort this out, I am going to sleep".

**Why:** a full rebuild per cycle is a build-system defect to fix once, not a cost to pay every
time; and five identical results mean the change is not reaching the code path, which one
in-process (shared-memory) run would have shown before any rebuild.

**How to apply:**
- Before the second rebuild of a loop, count "Building Fortran object" in the make log; if it is
  the whole tree, stop and find the trigger (suspects: `make install` re-touching the generated
  lib64/simple_args.f90 and SimpleGitVersion.h, which simple_defs.f90 #includes).
- When a change has no effect, verify the code path runs (shared-memory run, a log line) before
  rebuilding again.
- When the user asks for montages, render and send what exists first; analysis after.
Related: [[keep-test-loops-minimal]], [[never-rebuild-under-running-jobs]],
[[compile-and-test-means-batch2-fast-gate]].

2026-10-09 again (user: "please stop building from scratch. you need to make -j install"): ~/SIMPLE_clsrender/build_clsrender.sh
reran `cmake -S . -B build` before every `make install`; cmake embeds `git describe` (SIMPLE_GIT_VERSION, CMakeLists.txt l.64-77)
so every new commit changed a definition all modules see -> 646 of 646 objects recompiled per build. Fixed: the script now
configures only when build/Makefile is missing; rebuilds are plain `make -j install`. Check `grep -c "Building Fortran" make.log`.
ROOT CAUSE (found 2026-10-09): the `simple_git` make target runs scripts/insert_git_commit_hash.pl on every make; it rewrites
lib/simple/SimpleGitVersion.h whenever `git describe` changed, and simple_defs.f90 depends on it -> EVERY NEW COMMIT = full
rebuild (646 objects), locally too. Iterate on UNCOMMITTED changes: copy changed sources into the batch2 checkout
(scp into ~/SIMPLE_clsrender/src/...) and `cd build && make -j 24 install` (2 objects, 28 s); snapshot under a new name
(snap/<label>) so running jobs keep theirs; commit in batches when done.
