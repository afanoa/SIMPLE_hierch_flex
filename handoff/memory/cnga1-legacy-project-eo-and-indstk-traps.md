---
name: cnga1-legacy-project-eo-and-indstk-traps
description: "TRAP: today's SIMPLE refuses the Sept CNGA1 project (10_refine3D_auto): ptcl2D/ptcl3D halves differ (49.4 % agree) and the 274,333-image stack has 245,424 rows with a legacy indstk. Its indstk is PROVEN right; use cnga1_tree/input/cnga1_ready.simple + SIMPLE_EO_ADOPT3D=yes"
metadata:
  type: project
---

Found 2026-10-30 launching CNGA1 root trees (hetero_tree_20261028/cnga1_tree). Both errors stop any reconstruct3D or refine3D in seconds.
- **eo:** "particle 8 has different halves in ptcl2D and ptcl3D". ptcl2D holds the 2D stage's halves and ptcl3D refine3D_auto's; every
  CNGA1 project's ptcl2D has the 2D labels. Test switch SIMPLE_EO_ADOPT3D=yes (cnga1_tree build, eo_partition adopt_legacy) gives ptcl2D the 3D halves.
- **stack index:** "nptcls_stk not present". fix_projfile refuses ("no recorded indstk contract") because legacy indstk is distrusted by policy.
  - The stored indstk (3..274,332) is right: each row's CTF equals import row indstk(i) exactly, and Thon rings correlate 0.60 / 0.34
    read by indstk vs -0.13 / -0.34 read by row.
  - Test switch SIMPLE_FIX_TRUST_INDSTK=yes keeps it.
  - Fixed copy: cnga1_tree/input/cnga1_ready.simple, all fields verified identical. The original CNGA1 folders are untouched.

**Why:** the coordinator's first two CNGA1 launches failed on these; the 2022 stack (testbench/cnga1/final_particles.mrcs) was never consolidated.
**How to apply:** run anything on CNGA1 from cnga1_ready.simple with SIMPLE_EO_ADOPT3D=yes on the cnga1_tree build (or a build carrying
both diffs). Related: [[cnga1-is-the-apo-closed-dataset]], [[hetero-tree-passes-pgrp-to-flex-pca]].
