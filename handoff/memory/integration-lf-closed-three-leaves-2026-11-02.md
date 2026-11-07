---
name: integration-lf-closed-three-leaves-2026-11-02
description: "USER VERDICT 2026-11-02 on the combined-code 692k (f692_lf_s1): leaves r1312_g1_2 (1.04 closed / 0.63 open, 954 closed rows, purity .077), r2112_g1_2 (0.89 / 0.66, 1,259 rows, .095) and r1321_g1_2 (0.89 / 0.61, 739 rows, .061) are ALL closed conformers; the strict corner rule (>= 0.95 / <= 0.6) admits none of them and undercounts; closed state RETAINED on c5176bddc, in three lineages"
metadata:
  type: project
---

User, 2026-11-02 18:10, judging the 10 A unsharpened maps (final_maps/raw units 30, 50, 28 of f692_lf_s1, delivered beside the r17d corner
r1311_s2 = unit 52): "yes they are all closed".
- The combined code (claude/tree-lever-floor c5176bddc: cleanup, subspace site test, SVD discovery, draw-free split test, took-off pole rule,
  effect floor) keeps the RNase P closed state ([[rnasep-closed-state-must-be-retained]]) in three leaves from three root groups (r1 -> r13 ->
  r131; r2 -> r21 -> r211; r1 -> r13 -> r132), 2,952 closed-reference rows in them (r17d: one corner leaf with 816 + siblings), best 3-leaf
  union F1 0.137 vs the r17d corner's ~0.08.
- The scorecard's closed-corner rule fails all three on the open-arm residual (0.61-0.66 vs 0.58 for the old corner), so it must not be used
  alone as the hard gate: count a leaf as closed when its closed-arm occupancy is >= 0.95 OR the user's eye says so; report both arms.
**How to apply:** when gating the closed state, read the closed-arm occupancy and the closed-reference rows per leaf, not only the corner
flag; never fuse r1312 / r2112 / r1321 in a merge rule. Related: [[integration-692k-lf-launch-2026-11-02]],
[[integration-closed-two-leaves-2026-10-29]], [[integration-692k-r17d-closed-corner-2026-10-29]].
