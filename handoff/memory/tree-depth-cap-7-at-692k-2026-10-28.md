---
name: tree-depth-cap-7-at-692k-2026-10-28
description: "user 2026-10-28: RNA at 692k can't support ~200 states; never go past depth 7 (<= ~128 leaves). 692k trees already run maxdepth=7 (code default 4); 57 of build-7's 132 leaves sit AT the cap"
metadata:
  type: feedback
---

2026-10-28 ~17:00: "we cant really even support 200 states with 692k particles (on the RNA atleast). we probably shouldnt even go
past depth 7 (2^7 = 128 leaves)".

**Why:** ~5k particles per state is already the floor for useful RNA maps; more states only slice the same motions thinner
([[discrete-splits-capture-the-span-states-that-matter]]). The split test accepts nearly every split at this N, so the depth cap,
not the test, ends most lineages (build-7 base: leaves at node depth 4/5/6/7 = 3/21/51/57; 213 nodes incl. internal).
**How to apply:** keep maxdepth <= 7 on 692k RNA trees; "~200" in reports was NODES, not states - say leaves. Whether 5 or 6 loses
anything is a labels-only truncation replay (closed F1/arm per depth, qtoy2 truth per depth), not a new tree.
- VERDICT 2026-10-28 ~17:00 (labels-only replays, 17 s): RNase P base best closed group peaks at node depth 5 (F1 .153; d6/d7/leaves
  .146/.126/.093, union flat .156-.158); route A closed corner at d3 (r111 0.99/0.57), d4 0.96/0.66. qtoy2 (trees9): all 8 states
  by depth 4, NMI peaks at 4 (.659 order / .649 strongest); OLE seed ARI d2 .77 -> d3 .37 -> leaves .14. ADOPTED maxdepth=4 (code
  default; maxdepth=k delivers depth-(k+1) groups) for recipe12 gates + integration 692k. Base build-7 d0-d4 = 53 of 213 nodes.
  Outputs: beegfs hetero_tree_20261028/depth_cap/, fixtures/readouts/depth_replay_20261028_1647.out.
