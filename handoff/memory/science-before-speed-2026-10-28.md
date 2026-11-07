---
name: science-before-speed-2026-10-28
description: "STANDING (user 2026-10-28): 'we cannot sacrifice science for speed' -> a speed/CPU/disk lever counts only if the full-tree state gates hold; the 5-7 h goal never licenses dropping a state-finding step"
metadata:
  node_type: memory
  type: feedback
  originSessionId: 9643ce3d-fc7a-4dc1-860a-59db10867e7c
  modified: 2026-10-28T18:46:52.293Z
---

2026-10-28 ~15:00, right after the build-7 base 692k leaf scorecard failed the RNase P closed gate (F1 0.093) while the speed
work was centred on base-like recipe9/recipe11 trees: "we cannot sacrifice science for speed."

**Why:** a fast tree that loses the closed state (or other real states) is useless ([[rnasep-closed-state-must-be-retained]]);
the 5-7 h target ([[tree-runtime-goal-5-7h-at-700k]]) is a goal for the tree that HAS the science, not a reason to ship a lesser one.

**How to apply:**
- Split levers in two: pure engineering (same computation: cleanup, shared-memory small nodes, thread sizing, chaining steps,
  node-local stacks) vs science-touching (subsampled discovery, forced binary splits, fewer passes, stopping rules, dropping
  route A / site layer). Engineering levers need only a numerical-equivalence check; science-touching ones must not lose on any gate.
- The timing tree is the tree that passes the full-692k leaf scorecard (route A if it is what keeps closed), not the plain base.
- Report time and CPU next to the gate table, never a time win alone. [[tree-speed-over-null-precision]] (Monte-Carlo sizes) still holds.
