---
name: replicate-or-redesign-positives-only-2026-11-08
description: "USER 2026-11-08: 'If this is not replicable, the methods must be redesigned. Negatives are not publishable in high impact journals' -> judge methods by replicability across seeds; a non-replicable result triggers a redesign, not a write-up; papers carry positive, replicated results"
metadata:
  type: feedback
---

2026-11-08 ~23:00. I proposed writing up the day's negative results (what limits closed-state recovery) as publishable findings.
User: "well test reliability. If this is not replicable, the methods must be redesigned. Negatives are not publishable in high
impact journals. but beyond reliability, and considering the model itself how can we make it more robust".

**Why:** the paper needs a method that reproducibly recovers RNase P's closed state, not a catalogue of what failed.

**How to apply:**
- Every RNase P claim gets a replicate check: several seeds at the node or subtree where it is decided, with a hit rule set before
  the runs. A result that doesn't replicate means a method redesign.
- Negative results go into the testing guide as design constraints, never into a paper outline as findings.
- Robustness first: make the tree's hard decisions (seed site, K, accept) uncertainty-aware, so near-ties can't flip lineages
  (guide 7.10 thread 14).
Related: [[shared-tree-model-rnasep-direction-2026-11-07]], [[hetero-tree-paper-plan-2026-11-06]], [[rnasep-closed-state-must-be-retained]].
