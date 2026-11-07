---
name: depth-cap-replay-closed-peaks-d5-2026-10-28
description: label-only depth-cap replay of the 692k RNase P trees: best closed group peaks at depth 5 (base F1 .153) / 3-4 (route A, closed-corner map); depths 6-7 only fragment it -> a cap at 5 loses nothing on closed
metadata:
  type: project
---

2026-10-28, user asked how deep the 692k tree should go (doubts RNA supports > ~128 states). Replay from node rows only
(fixtures/scripts/depth_replay.py, beegfs hetero_tree_20261028): cut at depth D, every particle to its ancestor at D.
- build-7 base best closed group F1 by depth: d3 .085, d4 .099, **d5 .153** (r12112, 26.6k, map .82/.63), d6 .146, d7 .126,
  actual leaves .093; best union flat .156-.158 from d6 on. Never reaches the closed corner (best map .87/.68 at d4).
- route A: d3 r111 F1 .140 map **.99/.57** (closed corner), d4 r1111 .158 map .96/.66.

**Why:** 57 of the base's 132 leaves end at the depth-7 cap, so the cap, not the split test, sets the depth; the deep splits cut
closed into ever smaller groups (no extra closed signal).

**How to apply:** for closed-state questions a maxdepth of 5 is enough; compare deeper caps on the union F1 and maps, not
leaf counts. Related: [[site-pole-closed-state-replay-2026-10-27]], [[tree-runtime-goal-5-7h-at-700k]], [[rnasep-closed-state-must-be-retained]].
