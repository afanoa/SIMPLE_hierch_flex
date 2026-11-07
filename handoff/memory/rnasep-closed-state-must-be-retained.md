---
name: rnasep-closed-state-must-be-retained
description: "STANDING HARD GATE (user 2026-10-25): any tree change that loses the RNase P closed state is useless, whatever else it improves; the scorecard checks closed-conformer recall first"
metadata:
  node_type: memory
  type: feedback
  originSessionId: 25d47be4-7d1d-48d8-abdf-97b1c3c946ab
  modified: 2026-10-25T22:48:26.058Z
---

2026-10-25 about 18:50. The base80 scorecard (the improved baseline on m8c83t) recalled 34 of 55 deposited conformers against
27 for yesterday's 84-leaf tree, but 0 of 3 CLOSED conformers (yesterday's tree: 2 of 3). User: "well we need to retain the closed
state. If we are losing that then its uselleess."

**Why:** the closed state (P10.1-P12 region, the paper's key rare state; see [[rnasep-degenhardt-paper-closed-state-and-pipeline]])
is the biologically important result. Gains in recall, resolution or reproducibility don't make up for losing it.

**How to apply:**
- Every RNase P tree, and every change proposed for the tree, is gated FIRST on closed-conformer recall (paired arm / core test,
  scorecard section A) and on the closed share of particles (paper 11.6 %).
- A variant that drops closed recall below yesterday's tree (2 of 3, closed share 9 %) is rejected outright, even if other metrics
  improve.
- On toys (qtoy), the rare closed state's recall and purity are the first row.
- Related: [[tree-runtime-goal-5-7h-at-700k]], [[full692-tree-vs-deposited-benchmark-2026-10-22]].
