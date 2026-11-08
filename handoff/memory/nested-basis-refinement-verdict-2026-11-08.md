---
name: nested-basis-refinement-verdict-2026-11-08
description: "VERDICT: per-subtree refinement of the shared site basis (sppca pursuit refit per node, parent blocks fixed) adds no block on any fixture or null; qtoy2 1|2 and IgM a5 40/40 misses are proposal/test-limited, not basis-limited; frozen root basis stands"
metadata:
  node_type: memory
  type: project
  originSessionId: 1218bb2b-a21a-4d25-9e48-c40ddd89ceff
  modified: 2026-11-08T07:34:10.988Z
---

tree-latent, 2026-11-08 (batch2 `handoff/tree-latent/STATUS.md` addendum 4; `refine/summary.txt`). The test: at every node of the
frozen K=3 Sidak forests, sppca's pursuit was refitted on the node's y1 members (`tl_nested.py`):
- the parent blocks fixed as templates;
- the root's settings and stop;
- data y2, or y with the parents and nuisances profiled.

The result:
- No block entered at any of 30 fixture nodes (qtoyflex 4, qtoy2 15, IgM 10 % 11) or at the 4 null roots (IgM twin, qtoynull,
  olenull, micnull).
- The 2 SE gate passes many seeds, but the kappa 3 group penalty zeroes every grown block.

**Why:** the user wanted "refinements per subset". On these fixtures the misses are not basis-limited. qtoy2 r12's seed 36 gains
30.6 nats under the fitted parent prior but 2.1 with the parents profiled, so the gain lies inside the parent span. The frozen root
basis stands.

**How to apply:**
- Do not port per-subtree basis refinement. If it is revisited, use only the profiled-y form; on real data it needs the per-particle
  noise model first ([[real-data-fission-needs-misfit-model-2026-11-08]]).
- TRAP: a fitted-prior y pursuit on y1-selected members leaks (it added a block at qtoy2 r1 that y2 and profiled y did not).
- TRAP: micnull's y2 copy shows a 102-nat central seed that y does not.

Related: [[tree-latent-verdict-2026-11-07]], [[fission-select-members-on-y1-labels-may-read-y]].
