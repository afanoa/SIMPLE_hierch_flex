---
name: integration-692k-r17d-closed-corner-2026-10-29
description: 692k RNase P integration tree on recipe17d delivered one closed-corner leaf in 5 h 31 min; base/route A comparison and where the files are
metadata:
  node_type: memory
  type: project
  originSessionId: 3d793c94-9a6d-4a38-bcdd-0083174b8bb5
  modified: 2026-10-29T13:22:42.472Z
---

The 692k RNase P integration tree on recipe17d is done. Recipe17d = route A to depth 2 (kway_maxdepth=2), choose, the locality guard, maxdepth 4, row-local projects, reuse/prune, the STATES_DELIVERED marker. The run used nodepar 12, nparts 12, the split env maps/0.03 and all 692,392 rows. It ran 2026-10-29 03:42 to 09:14.

- Cost: STATES_DELIVERED at 5 h 31 min, inside the [[tree-runtime-goal-5-7h-at-700k]] target. 873 core-h, 86 GB final disk (at least 176 GB peak), 84 nodes, 67 leaves.
- Closed lineage: r1 → r13 → r131 → r1311, with the closed arm at 0.45 → 0.67 → 0.87 → 0.90.
  - r1311_s2 is the one closed-corner leaf: 12,010 rows, 816 of the 8,720 reference rows, 1.00 / 0.58.
  - Sibling r1312: 21,870 rows, 0.84 / 0.56, best closed-leaf F1 0.119. The re-test refused to split it.
- Build-7 base: 132 leaves, F1 0.093, no corner leaf.
- Landscape at depth 4 sits closest to route A's cut: W1 0.017, against 0.021 to the base. The comparison used the integration tree's root embedding, because the build-7 embedding and leaf rows were lost with scratch. Route A's depth-4 cut was rebuilt from node rows.

**Why:** this is the first full-size tree that keeps the closed state ([[rnasep-closed-state-must-be-retained]]) at production cost.

**How to apply:** use it as the recipe17d baseline. The files are in hetero_tree_20261028/integration/ (final/, readouts/, logs/).

Related: [[kway-root-binary-below-2026-10-28]], [[choose-routing-fixes-closed-dilution-2026-10-28]].
