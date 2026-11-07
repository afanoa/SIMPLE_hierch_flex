---
name: final-routing-pass-verdict-2026-10-29
description: "Final routing pass (SIMPLE_TREE_ROUTE, route_v5.diff, uncommitted, worktree tree-final-routing): flat refit over leaves lifts qtoy2 NMI .706 -> .827 (4 it); on RNase P closed lineage it churns low-margin rows, no closed gain, corner 1.00 -> 0.97; keep default off"
metadata:
  node_type: memory
  type: project
  originSessionId: 9643ce3d-fc7a-4dc1-860a-59db10867e7c
  modified: 2026-10-29T23:05:16.494Z
---

2026-10-29, routing session. Held-out scoring (each half vs maps from the other half), refine3D prob_neigh, volrec=no; margin = best minus second
log-evidence. Build route/build_v5 (gate 15/15).
- qtoy2 (13 leaves): route-down own-leaf .83, NMI .717 (near neutral); flat refit NMI .770/.789/.809/.827 (it1-4, still rising), every combo
  F1 up, moved rows right .50 -> .81; two runs agree >= 99.76 %. 4 flat iterations 10.7 core-h.
- RNase P closed lineage (route/f692_lineage_v3, 34.7 core-h): low separability (median margin ~2 nats vs qtoy2 root 10.5; 45 % of rows
  < 1 nat). Flat over a 9-leaf pool swaps ~40 % of the corner's rows per iteration with unchanged closed-ref counts; lineage recall
  .353 -> .327, purity .069 -> .065; corner arm 1.00/0.58 -> 0.97/0.59 (still a corner). r2's closed refs mostly stay in r2; r2322 stays
  a separate sub-state. A margin gate removes the qtoy2 gain.

**Why:** the user asked (2026-10-29) for "one final pass where we route each particle down the tree to see if they land in the correct one".
**How to apply:** route-down = diagnostic (where disagreements start, per-particle margin), not a refit. Flat refit helps well-separated
discrete states, not RNase P's low-margin closed lineage; no full 692k flat (~300 core-h) without new evidence. Related:
[[integration-closed-two-leaves-2026-10-29]], [[rnasep-closed-state-must-be-retained]].
