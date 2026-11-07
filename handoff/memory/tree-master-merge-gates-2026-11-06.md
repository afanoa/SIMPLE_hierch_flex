---
name: tree-master-merge-gates-2026-11-06
description: "master merged into the tree (claude/tree-master-merge @ ee26d9bbf) passes the gates at master's default prob_beta; prob_beta=1 is the only arm on record that loses the RNase P closed corner"
metadata:
  node_type: memory
  type: project
  originSessionId: cc23038c-fea9-4b5f-b321-d58eb2942350
  modified: 2026-11-06T21:34:56.547Z
---

Merge of master (posterior machinery, memory/speed work) into the hetero tree: branch `claude/tree-master-merge`, worktree
`.claude/worktrees/tree-master-merge`, built on batch2 at `math_reworks/merge_master2` (fast gate 15/15). Gate arms 2026-11-06,
seed 2, keys as the bc baseline (`math_reworks/gates/<arm>_*`, readouts `gates/<arm>_readout.txt`):

- mm2 (master's default prob_beta): qtoynull, olenull, micnull 1 leaf each; qtoy2 8 leaves, NMI 0.780 (bc 0.812, 25 earlier arms
  0.71-0.87); r131 closed corner kept (r2 SEPARATES, map 1.03/0.52); OLE 8 leaves. PASS.
- mb1 (prob_beta=1): qtoy2 NMI 0.812, but r131 never reaches the closed corner (best 0.94/0.46), the only arm on record
  without it. Fails [[rnasep-closed-state-must-be-retained]].

**Why:** the merge is gated; the next step is the user's push to the fork ([[never-push-to-origin-without-permission]]).
**How to apply:** do not run the tree with prob_beta=1; judge qtoy2 NMI against the 0.71-0.87 arm spread, not one baseline run.
