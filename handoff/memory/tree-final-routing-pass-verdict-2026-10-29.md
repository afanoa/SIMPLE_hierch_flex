---
name: tree-final-routing-pass-verdict-2026-10-29
description: "hetero_tree final routing pass (SIMPLE_TREE_ROUTE, default off, uncommitted): flat refit wins on qtoy2 (NMI .706 -> .827 at 4 its), no closed gain on RNase P (margins ~2 nats, corner 1.02 -> 0.97); no default change"
metadata:
  type: project
---

Built 2026-10-29 (worktree tree-final-routing, simple_hetero_tree.f90 only, uncommitted; cluster route/route_v5.diff, route/build_v5).
SIMPLE_TREE_ROUTE=only (offline on a finished tree via route/route_run.sbatch, mirror dir) or =yes (after STATES_DELIVERED); _NODES,
_LEAVES, _ROOT, _FLATITS scope it. Scoring = stage B's refine3D prob_neigh geom, one iteration, volrec=no, per half against
other-half maps; margin from prob_state_evidence.txt.
- qtoy2: route-down NMI +.011 only; flat .770/.789/.809/.827 (it1-4), still rising; reproducible >= 99.8 %.
- RNase P 692k closed lineage: r2's closed refs stay on r2 (11 % return, offset by r1 -> r2); r2322 keeps its refs; lineage recall
  .353 -> .33 under flat; corner stays a corner (0.97/0.59) but softens. Margin gate does not help either dataset.
**Why:** the user asked for "one final pass where we route each particle down the tree to see if they land in the correct one."
**How to apply:** do not make it a default; full 692k flat (~300 core-h) not recommended on this evidence; the user decides.
Notes: side-session-notes tree_final_routing_2026-10-29.md. Related: [[integration-closed-two-leaves-2026-10-29]],
[[rnasep-closed-state-must-be-retained]].
