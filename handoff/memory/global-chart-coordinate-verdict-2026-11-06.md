---
name: global-chart-coordinate-verdict-2026-11-06
description: "VERDICT (2026-11-06): the composed global chart coordinate (child -> parent affine maps on fflex charts) is a piecewise-affine atlas with SEAMS; its split-axis valleys and cross-node leaf distances are artifacts, it does not isolate the closed minority on real data (d' <= 1.3; 692k closed pair among the closest sibling pairs) and reproduces worse than the root chart (.45-.96 vs 1.000 across recipes); nothing to port"
metadata:
  node_type: memory
  type: project
  originSessionId: 5482e66b-b4b2-4b69-b534-a7b318bf3360
  modified: 2026-11-07T03:42:37.965Z
---

Tests (scripts/hetero_tree_proto_latent/f2r_global.py on mest_charts.py's composition; runs math_reworks/latent/global/) on the mm2 /
pb15 gate trees and the 692k trees rnasep_klem/stree_port/full692_{m33,base92,base80}_s1; tables in the side-session note
f2r_flex_to_refinement_status.md.
- Minority: closed AUROC global / root chart: r131 .67-.69 / .62 (the root-latent ceiling), qtoy2 .994 / .982, 692k .75-.82 / .77.
  Closed-enriched leaves are no more isolated than other sibling pairs (r131 d' 1.35-1.51; 692k closed pair r1111221_s1/_s2 d' 1.12,
  lower quartile of 112 pairs: a latent merge rule would fuse it, against [[rnasep-closed-state-must-be-retained]]).
- Reproducibility: the root chart reproduces at .998-1.000 (same input -> same flex fit); the composed coordinate at .89 across r131
  seeds (leaf NMI .19) and .45-.96 across 692k recipes (NMI .24-.45); extrema overlap .13-.80 (chance .05).
- SEAMS: particles of different children reach the root frame through different maps (R^2 .1-.6). Valleys of .8-.9 on r131 and the
  qtoyflex continuum vanish in the root chart (one frame): they are chart seams, as are cross-node d' of 20-54 on 692k. Only qtoy2's
  discrete conformers show a real, half-reproducible root-chart valley (.6).

**Why:** the user asked whether a global embedding can say when to stop splitting / merge up and find a minority (closed state).
**How to apply:** do not use the composed coordinate for between-node decisions or density tests; within-node (sibling) statistics
and the root chart are seam-free. A continuous atlas would need overlapping particle sets between neighbouring charts (LLC-style),
which the tree's disjoint nodes do not provide. Related: [[tree-charts-coordination-gap-verdict-2026-11-06]],
[[flex-to-refinement-coupling-verdicts-2026-11-06]], [[landscape-must-be-reproducible-and-splits-data-derived]].
