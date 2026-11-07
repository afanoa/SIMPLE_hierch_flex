---
name: no-4m-synthetic-tests-too-expensive
description: "user 2026-10-19: the 4M-particle synthetic 55-conformer set is too expensive to test on and hogs the CPU cap; cancelled, do not relaunch; nparts=40 was called unacceptable (nparts=10 nthr=8 rule holds for everything)"
metadata:
  type: feedback
---

2026-10-19: the 4M synthetic scaling run (stack merge 410 GB, refine3D master > 130 GB, 40 parts) was
cancelled on the user's instruction: "Its too expensive to test on 4M set, and we are using too much
cpus for that (I could use those on more useful tests)". The same day nparts=40 on its refine3D was
called "unacceptable".

**Why:** the CPU cap is shared with the user's own tests; a 4M run holds hundreds of cores for days
for a scaling data point that 200k/700k already approximate.

**How to apply:** the scaling question (state count vs N) is answered with 200k and 700k only, as
half-tree pairs under the cross-tree rule ([[klem-cross-tree-merge-design-2026-10-19]]). Never
submit a 4M stage again without being asked. Every refine3D stays at nparts=10 nthr=8
([[simple-run-parallelism-preference]]); check the watcher scripts for stale 4M triggers before
re-arming them.
