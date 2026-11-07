---
name: tree-speed-over-null-precision
description: "user 2026-10-27: offered 20-shuffle revert vs a sequential 50-shuffle null for the split re-test, said 'Just go for whatever is faster' -> pick the faster option for tree statistics' Monte-Carlo sizes"
metadata:
  node_type: memory
  type: feedback
  originSessionId: 0f67cd42-16a3-4392-a859-a59fdac654fc
  modified: 2026-10-27T14:57:48.222Z
---

2026-10-27, after I brought the choice for the 1.81x split slowdown (design-layer re-test on 50 null shuffles): option A revert
to 20 (always fastest; decisions near z_crit can differ, 5.8 % at 1 null-sd), option B sequential 20->50 (agrees with 50 to
0.24 %, usually 20). The user: "Just go for whatever is faster. and survey for other optimizations as well".

**Why:** tree wall time gates development ([[tree-runtime-goal-5-7h-at-700k]]); a Monte-Carlo null size is not a design
feature worth wall time.

**How to apply:** when a tree speed-up trades only Monte-Carlo precision (shuffle/permutation counts, null sizes), present
the trade-off in one line and default to the faster option. Scientific rules (minpop, split acceptance, small-state retention
[[rnasep-closed-state-must-be-retained]]) are not covered by this.
