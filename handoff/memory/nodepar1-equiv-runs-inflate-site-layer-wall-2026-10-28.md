---
name: nodepar1-equiv-runs-inflate-site-layer-wall-2026-10-28
description: "TRAP: hetero_tree site-layer (fs) jobs share one global nodepar cap; at nodepar=1 a node's home site ranks serialize and hold seeding, so nodepar-1 equivalence runs overstate fs-lever wall gains (13 min on r111 md1, ~0-2 min/level at 692k)"
metadata:
  node_type: memory
  type: project
  originSessionId: c8c6f555-6cba-44a4-9ec4-48197b959f2f
  modified: 2026-10-29T02:55:29.218Z
---

In the hetero_tree controller every factorial-site job (home ranks, climbs, refines, maps) counts against a single nodepar cap
(fs_nrun), and fs_home_ranked holds a node's seeding until its own k>1 home ranks finish. At nodepar=1 those ranks run one after
another: 5-9 ranks of ~1 min each. So recipe15b's rank reuse cut r111 md1 wall from 76 to 63 min. Build 7 ran about 8 fs jobs at a
time with 4-5 home ranks per node, so at 692k nodepar 12-16 the same lever saves mostly CPU (52 of 101 site-layer core-h).

**Why:** a wall-time gain measured at nodepar 1 can be a queue artifact, not something that holds at scale.
**How to apply:** take wall gains of site-layer levers from runs at production nodepar, or model the fs queue. Check the CPU and
same-computation verdict at nodepar 1. Also: cmp between runs in different folders fails on absolute paths in .simple and text
outputs; normalize the run path first. Related: [[hetero-tree-numerics-sensitivity-of-engineering-levers-2026-10-28]],
[[m33-tree-wall-profile-2026-10-25]].
