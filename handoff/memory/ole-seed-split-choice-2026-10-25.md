---
name: ole-seed-split-choice-2026-10-25
description: "OLE seeds 1/2 agree only at the root (unit ARI 0.18, no unit at Jaccard >= 0.5); r1 is decided by site ORDER (first kept pole), r2 is one axis cut at opposite ends (a swing); site_pole_rule row index is NOT a site_seeds state"
metadata:
  node_type: memory
  type: project
  originSessionId: 58de2c8f-ad53-44fb-b8a2-05aadff16ec4
  modified: 2026-10-25T16:10:30.045Z
---

2026-10-25, for the coordinator "2D to 3D heterogeneous state translation" (worktree sharp-proskuriakova-b68280).

- OLE (ole_12707) 6_hetero_tree (seed 1) vs 7_hetero_tree (seed 2), stree_port/ole_seed_overlap_s1s2.txt (seed_overlap.py):
  - 14 vs 12 units, unit ARI 0.182, NMI 0.341. No unit or leaf of either seed has a partner at Jaccard >= 0.5.
  - The root split reproduces (ARI 0.945); r1 and r2 split along different axes (ARI about 0).
  - So OLE fails the standing whole-landscape check below depth 1 ([[landscape-must-be-reproducible-and-splits-data-derived]]).
- Why, from each node's kept poles (all sites ranked by the factorial layer's fsite_ dirs):
  - r1: both seeds split on their site 2, a seed-PRIVATE axis (|r| 0.003 between the two, <= 0.012 with anything in the other seed).
    Sites 3 and 4 recur in both seeds (|r| 0.92 / 0.99) and were never used.
    The tree takes the FIRST site in order that keeps a pole, not the strongest: seed 2's winner had the weakest margins (5.4 / 6.4 vs 28 / 32, 40 / 22).
  - r2: same axis in both seeds (|r| 0.885), but the poles sit at OPPOSITE ends (z3hi vs z3lo). Each seed cut a different end off one continuous axis, so ARI -0.015: a swing sliced twice.
- Trap: site_pole_rule.txt's row index is NOT a site_seeds state. After a rejection, pole_rule re-seeds the kept poles at posterior >= 0.5 and relabels them.
  - Example: cdisc M r1 site 3, rule row 1 has 32,263 seeds, and site_seeds pole 1 has 182.
  - The kept pole = the non-rest site_poles.txt line with the most site_seeds rows (factorial's null_score3.py).
- Response, all default off and uncommitted in the working tree, snapshot m8c67t:
  - node_motion=yes: the factorial switch / swing call per kept pole, with multiway product cells;
  - strongest_pole=yes: the largest min margin seeds the split;
  - re-entry once per lineage;
  - site_pole_rule gains the post-rule state and its seeds (two columns before the kept word).
  - Replays queued: r1 strongest (both seeds), r2 node_motion (both seeds), then r1 node_motion.

- Replays, both seeds (stree_port/motion_ole, arm_score.py, two_pole_table.py):
  - Cross-seed ARI under each arm:

    | arm | r1 | r2 |
    |---|---|---|
    | site order (as run) | 0.006 | -0.015 |
    | strongest_pole | 0.393 | 0.235 |
    | + two_poles (the opposite pole by endpoint name, else the same site's opposite end at posterior >= 0.5) | 0.538 | 0.418 |
    | m8c75t from fflex: component sort + design (node_motion, motion_multi=no, gate, inherit, strongest, two_poles) | 0.668 | 0.772 |

  - The last row is stree_port/flexsort/ole_s{1,2}_r1r2, scored with node_children_ari.py and tree_pair_verify.py. Both seeds split r1 on the same axis (r +0.963) and r2 on the same axis, same orientation (+0.855). The gain combines the sort, the fflex refit and the design.

  - node_motion calls every OLE candidate a switch (nb_SF 0.001), matching factorial. minpop folding turns the multiway back into pole-vs-rest, so node_motion leaves r2 unchanged.
  - Two-pole aligns the children's mean axis positions across seeds. The residual is the boundary rows of a graded axis: the ends reproduce, a binary cut does not.

- 7_tree_MF (m8c63t, split contrast maps, flex_min_neff 1000, seed 1) vs 6_hetero_tree, same seed: unit ARI 0.169, NMI 0.309, 0 of 14
  units at Jaccard >= 0.5, best-partner unit-map deviation corr median 0.77 (min 0.66). Seed 1 vs 2 (6_ vs 7_): ARI 0.182, NMI 0.341,
  deviation corr median 0.79 (min 0.42). So the partitions differ like a reseed while the maps stay within the seed spread
  (stree_port/unit_overlap.py; unit maps = leaf-weighted means of final_maps/emp2).

**Why:** the user's rule is that landscapes must reproduce across seeds; split choice by site order makes deep splits arbitrary.
**How to apply:**
- Judge any split-choice change by the cross-seed ARI of the same node.
- Read kept poles via site_poles + site_seeds, never the rule's index.
