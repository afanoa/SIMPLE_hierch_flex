---
name: mergeup-fuses-siblings-shown-distinct-2026-10-27
description: "state_mergeup fuses sibling leaves their own node had shown distinct (OLE s1 10/11, s2 10/10, RNase P m33 16/20 incl. the closed conf18/conf19 pair); merge-up units lose 6-25 % of OLE span; proposal = merge only on demonstrated equivalence (doc planned/hetero_tree_adaptive_split_merge.md)"
metadata:
  node_type: memory
  type: project
  originSessionId: 01619764-7985-49aa-b2e9-539d9cd5a16f
  modified: 2026-10-27T15:01:02.177Z
---

2026-10-27, an offline replay from saved split and merge-up logs. Worktree adaptive-split-merge-b384e9; doc
`doc/implementation_notes/planned/hetero_tree_adaptive_split_merge.md`; scripts in
`scripts/hetero_tree_math_experiments/split_merge_replay/`. All uncommitted.

- The merge-up re-tests sibling pairs that the node test already called distinct (beyond-EM re-test, own mask). It then fuses
  most of them, because the repeated test is weaker:
  - "not shown different" is treated as "same";
  - pooled 3-35-unit tests run inside the union of the units' masks (closed node 2.9 % vs level-7 union 17.7 %);
  - the threshold is the max over 60-11,900 null statistics;
  - the pooled EM drains small groups to 256-344 per half.
- OLE r12111: identical inputs gave distinct at the node (T 0.32 vs M 0.235) and a merge in the merge-up (0.23 vs 0.339). At
  about 500 per group per half the call is a coin toss.
- The consistency margin (the tree's smallest demonstrated effect r_hat) shows 0 of 17 / 0 of 14 OLE merges equivalent, and
  3 of 49 on RNase P (all from pooled tests).
- OLE K* = 1 stage-B pairs: splitting them is NOT supported (7 of 10 have the same counterpart in the other seed).
- m33's node splits had the em_min 1000 floor: 13 of its 23 K* = 1 nodes are drains with no re-test.

**Why:** the user wants states that represent the landscape, and small states are fine
([[discrete-splits-capture-the-span-states-that-matter]]). The merge-up costs span and reproduced states, and on RNase P it
cost the closed pair ([[rnasep-closed-state-must-be-retained]]).

**How to apply:**
- Use the leaves, not the merge-up units, as the state set until the merge change is built.
- Any merge rule must not overrule a node's own distinct verdict.
- The `minpop` replacement (power rule) is untested and needs a native A/B ([[tree-design-feature-freeze-2026-10-25]]).
