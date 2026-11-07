---
name: depth-wiener-no-truth-map-benefit-2026-11-07
description: "depth-tied tree_wiener (depth_wiener=yes, port piece 1) has no truth-map benefit on qtoy2/qtoyflex; emp2 is equal or a hair better per band; its only claimed gain is nominal resolution on deep trees"
metadata:
  node_type: memory
  type: project
  originSessionId: 7da7e826-c848-444f-a66e-ad78379bb739
  modified: 2026-11-07T06:27:06.590Z
---

depth_wiener (piece 1 of the 2026-11-06 ports, commit d142ea8af on claude/tree-ports, gate-matched on labels) was tested
against truth on 2026-11-07 with `map_truth_error.py` (handoff/ports/scripts/hetero_tree_mest on batch2): every final-map
estimator of a gate tree (final_maps/{raw,own,tree2,emp2,depth2}/unit_NN.mrc) vs the leaf's composition-weighted truth
conformer maps (routes/toy55/maps160, tables from qtoy2.csv / qtoyflex.csv), FSC and relative error per band.

- qtoy2 (8 leaves): emp2 error .264 / .543 / .820 (40-20 / 20-14 / 14-10 A) vs depth2 .264 / .553 / .836; raw .966 at 14-10 A.
- qtoyflex: emp2 .278 / .582 / .863 vs depth2 .277 / .584 / .867 (tree2 the same).
- The depth-tied covariance R^2 is .04 on these depth-1 trees (.53 on r131's three levels).

**Why:** the ship decision on depth_wiener should rest on map accuracy against truth, not on the nominal half-map resolution
(4.7 -> 4.3 A on 692k) that motivated it; on what truth can measure, emp2 >= depth2.
**How to apply:** do not claim a map-quality gain for depth_wiener=yes; if it ships, it ships for the deep-tree resolution
readout only, and the default is the maintainer's call. Related: [[mest-formulations-verdict-2026-11-06]],
[[hub-hetero-tree-current]].
