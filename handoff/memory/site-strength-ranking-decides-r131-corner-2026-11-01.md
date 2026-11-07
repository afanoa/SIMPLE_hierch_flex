---
name: site-strength-ranking-decides-r131-corner-2026-11-01
description: TRAP - at the r131 root the closed corner at md1 hinges on which site strongest_pole ranks first (site 1, 2981 seeds = the closed population), not on the pole tests
metadata:
  type: project
---
2026-11-01 math reworks: the subspace site test (rework 1) found the same two root poles as the halves test on the r131 fixture
(site 1: 2981 seeds; site 3: 9076 seeds), but its seed-set Fisher z ranked site 3 strongest (4.66 vs 3.29) where the halves margin
ranked site 1 (12.2 vs 7.0). Seeding from site 3 shifted the root partition, r1's re-test dropped from z 4.31 to 2.45 (< 3.71) and
the 1.02/0.49 closed-corner stage-B state was not delivered as a leaf: HARD GATE FAIL by the ranking alone. Ranking by the weakest
endpoint's hold-out z (5.8 vs 5.1) restores site 1 (variant B, sub2 arm).
**How to apply:** any change to site_rank's margins or to strongest_pole must be read at the r131 root as "which site seeds" first;
report the site order with the corner verdict. Related: [[hetero-tree-math-reworks-setup-2026-11-01]], [[rnasep-closed-state-must-be-retained]].
