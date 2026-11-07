---
name: floor-off-guard-verdict-2026-10-29
description: "ADOPTED 2026-10-29 (user: drop the floor, guard decides; shrink shelved). recipe17d floor-off test: with SIMPLE_SPLIT_RHOMIN=0 + locality guard every null still merges (L/L_ref 0.20-0.24); only real change = integration r24 splits (open-side, guard 0.705, no closed signal); ancestral shrink fails r131 via the floor"
metadata:
  node_type: memory
  type: project
  originSessionId: 9643ce3d-fc7a-4dc1-860a-59db10867e7c
  modified: 2026-10-29T19:56:42.776Z
---

2026-10-29, factorial/integration session, hetero_tree_20261028/floor_off/ (runs, replay, readouts; ~59 core-h).
- The floor is env-only: simple_state_mergeup.f90 rhomin defaults 0, read from SIMPLE_SPLIT_RHOMIN; production exports 0.03.
- Floor off, guard on: qtoynull s1/s3, olenull s1/s3/s4 all pass the re-test and are merged by the guard at 0.20-0.24.
- m33 r2121212 was never a floor case: in-tree re-test fails (z 0.98 vs 3.71); the fixture (m 3) passes, then the guard merges at 0.370.
- Integration r24 (33,385): floor refused it; floor off it splits 17,309/16,076, guard 0.705, both open-side (closed arm 0.40/0.43),
  closed refs 154/173 (no enrichment), corr AUC 0.616. Integration r1321 (also floor-refused) untested: stage-B halves pruned.
- Fixtures use their own-root Bonferroni (m 3-52), looser than the 692k m 460; replays use the tree's m.
- Ancestral shrink prototype (resolution session, 2026-10-29): r131 root z 14.48 passes but rho 0.0298 < 0.03 -> one leaf (off arm
  splits, guard 0.93). Shrink removes the 12-10 A band the held-out re-test reads. OLE r32 node replay degenerate (seed-invariant,
  no split on recipe17d even off).

**Why:** the user asked (2026-10-29) whether the fixed 0.03 floor can go now that the guard exists (open question #3).
**How to apply:** USER 2026-10-29: "drop the 0.03 floor ... the guard decides. Shelve ancestral shrinkage for now." Production stops
exporting SIMPLE_SPLIT_RHOMIN (unset, not 0); no source change (default 0). Shrink diff/install kept, no further runs. Between-band guard cases
(0.37, 0.705) make 0.5 the deciding line. Related: [[null-split-diffuse-real-split-localized-2026-10-28]],
[[recipe9-lineage-split-test-needs-env-2026-10-28]], [[tree-data-derived-lp-verdict-2026-10-29]].
