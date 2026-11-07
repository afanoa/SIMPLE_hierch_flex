---
name: test-on-subtrees-not-full-692k-2026-10-28
description: "STANDING (user 2026-10-28): 'we cant really run the whole tree each time' -> test changes on subtrees of interest (closed-holding node to leaves, r1211, small fixtures); full 692k only for the integration build"
metadata:
  node_type: memory
  type: feedback
  originSessionId: 9643ce3d-fc7a-4dc1-860a-59db10867e7c
  modified: 2026-10-28T18:48:43.365Z
---

2026-10-28 ~14:50, while recipe11 asked for a 692k timing tree and three 692k trees had run 5-14 h: "we cant really run the
whole tree each time. So maybe limit testing to specific sub trees of interest."

**Why:** a 692k tree costs 8-14 h and hundreds of CPUs per answer; most changes act below depth 2, where a subtree from a saved
node's rows gives the same verdict in a fraction of the time ([[science-before-speed-2026-10-28]] still applies to the verdict).

**How to apply:**
- Standard set: (1) the closed subtree = the node that first isolates closed (r11, ~80k rows) run to leaves + leaf_scorecard.py on
  its rows; (2) r1211 fixture; (3) qtoy2 + qtoynull + olenull; (4) OLE 80k for reproducibility; (5) one small-node subtree
  (~20-30k node) to time engineering levers.
- A change that acts at the root/depth 1 (discovery, route A seeding, strongest_pole) needs a root-level test, not a subtree below it.
- Full 692k only for the one integration build after the set passes ([[integration-plan-2026-10-28-builds]]).
- FEEDBACK 2026-10-28 ~17:00: "why do I have to wait 7 hours for results ... This is disappointing." A maxdepth-4 subtree (23 nodes,
  nodepar 3-4, ~5 h) is NOT the loop. Fast loop: (1) minutes - replay decision steps (split test, merge, seeding/pole choice) on saved
  node outputs; (2) ~40 min - maxdepth 1 runs on the closed-lineage fixtures (r112, r1211) + qtoy2/nulls md1, all in parallel;
  (3) report each node's closed readout as it finishes, never only at the end; (4) a subtree to leaves only to confirm a winner,
  with nodepar >= nodes per level. Keep closed-lineage node dirs (r11/r112/r1121..., r12/r121/r1211...) for replays.
