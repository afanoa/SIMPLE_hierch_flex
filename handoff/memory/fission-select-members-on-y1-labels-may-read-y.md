---
name: fission-select-members-on-y1-labels-may-read-y
description: "RULE (tree-latent 2026-11-07): in a fission-tested tree, node membership and child weights must be selected on fission copy X1; split proposals and sat5's input labels may read the full data y"
metadata:
  node_type: memory
  type: project
  originSessionId: 1218bb2b-a21a-4d25-9e48-c40ddd89ceff
  modified: 2026-11-08T02:02:19.187Z
---

In a tree whose splits are accepted by sat5's fission saturation test (igm-refs Part 5), two kinds of input differ:
- **Node membership** (which particles' y2 data enter the test's maps) is a selection step. If children are selected on y, their y2
  maps carry the parent's selection, because y2's noise correlates .71 with y's. That compresses the near-boundary step and makes a
  re-carved continuum read DISCRETE. So membership and child weights must come from y1.
- **Split labels** only set sat5's score direction D, a global map. Each particle's group comes from its own y1 score, so labels
  from y are valid; sat5's random-split control and Part 5's round-1 labels already relied on this.

**Why:** the strict variant (everything on y1) lost qtoyflex's closed leaf (F1 .337 vs .558) and IgM a4. The rule above (v3)
recovered both, with nulls still one leaf.

**How to apply:** when porting the tree or a fission split test, draw the per-particle fission once and select every subset on X1.
Proposals and labels may use full-data latents, and final labels may be re-assigned on y under the same deconvolved mixture.
Related: [[tree-latent-verdict-2026-11-07]], [[data-fission-merge-test-verdict-2026-11-07]].
