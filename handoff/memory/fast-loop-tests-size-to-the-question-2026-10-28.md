---
name: fast-loop-tests-size-to-the-question-2026-10-28
description: user wants results in ~40 min, not 5 h subtrees: maxdepth 0 for a node's own split (35-40 min at 79k), maxdepth 1 doubles it; decision-step changes replay on saved stage-B outputs in minutes; report per node as it lands
metadata:
  type: feedback
---

2026-10-28 the user objected that a 5 h closed subtree "is not the fast loop" they asked for; standing rule since: no full 692k
tree per test, changes are tested on subtrees of interest (fixtures), and results come within ~40 min, node by node.
- A node's own split: maxdepth 0 (the fixture root only), about 35-40 min and ~19 core-h at 79k rows; maxdepth 1 doubles it
  (r112 fast runs: 74-78 min, ~38 core-h, 3.5 GB).
- Routing / split-test / merge / seeding changes: replay on saved stage-B outputs (as adaptive-split-merge did for choose), minutes.
- Report each node's closed readout the moment it finishes; don't wait for the run to end.

**Why:** the user develops by quick iterations; long runs that report only at the end stall the loop.

**How to apply:** size every test to the question before launching; say the expected wall up front; see
[[tree-runtime-goal-5-7h-at-700k]], [[depth-cap-replay-closed-peaks-d5-2026-10-28]], [[do-the-requested-job-not-diagnostic-sprawl]].
