---
name: no-new-builds-use-existing-hooks
description: "\"you are wasting time building\" (2026-10-14): before writing Fortran + building for a test, check memories/code for an existing hook or program that already does it (SIMPLE_COV_COMPOSE did the averages-basis -> particle embedding I rebuilt as SIMPLE_FLEX_BASIS_IMPORT)"
metadata:
  node_type: memory
  type: feedback
  originSessionId: 34b66a02-03d2-42df-91f8-915afbf80120
  modified: 2026-10-14T20:37:53.937Z
---

2026-10-14: I added SIMPLE_FLEX_BASIS_IMPORT to flex_pca and built reseed31 to embed particles in the basis learnt on
the expanded averages. The user: "you are wasting time building...". The flex2D/master builds already carry
SIMPLE_COV_COMPOSE, which composes another run's basis and embeds the particles in it. It had been tested that same day
on RNase P (see [[rnasep-averages-basis-propagation-2026-10-14]]).

**Why:** each build plus the code edit costs the user 10-30 min of waiting. A test that an existing hook or program
covers should run on the existing binary right away.

**How to apply:** for a new test idea, grep the memories (hub-simple-infrastructure, hub-flex-pca-method) and the code
for an existing env hook, program or script first. Run with what exists. Write and build new code only when nothing
covers it, and say so in one line before starting. Related: [[when-pointed-at-existing-code-implement-dont-probe]],
[[batch2-compose-hook-and-worker-names]].
