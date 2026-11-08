---
name: fission-split-test-verdict-2026-11-07
description: "VERDICT 2026-11-07: fission + locality guard in state_split passes every gate null (tau 0.6 and 1; the effect floor never acts) and keeps qtoy2 (NMI .717 vs .811), but FAILS the r131 hard gate at both tau: the r131 root is a near-tie unrelated to closed, too weak for an honest test, so the depth-1 corner is never reached"
metadata:
  node_type: memory
  type: project
  originSessionId: cc23038c-fea9-4b5f-b321-d58eb2942350
  modified: 2026-11-08T01:37:16.611Z
---

fission-split2 session; `hetero_tree_20261028/handoff/fission-split/{VERDICT,STATUS}.md`, patch v3 against cc6d9cbb6
(prototype env vars `SIMPLE_SPLIT_FISSION=<X2 project>`, `SIMPLE_SPLIT_GUARDS=no|locality`); runs `fission_split/runs/`.
- Design: the whole tree runs on an X1 project (all selection on X1); only the split decision reads X2.
- Gate nulls with fission + locality guard:
  - qtoynull: z -0.19 at tau 0.6, 0.42 at tau 1. Self-selection; X1 maps at the same labels read z 4.10.
  - micnull: z 0.41 / 0.80. Self-selection with a weak pose tilt; X1 maps read z 61.
  - olenull at tau 0.6: fission keeps a real contrast split (z 7.78; the small EM group is brighter in truth) and the locality
    guard merges it. At tau 1, fission alone gives one leaf.
  - The effect floor never acts.
- qtoy2 at tau 0.6: 8 leaves, NMI .717 (pp_on .811); the 1/2 and 22 pairs merge.
- r131 fails the hard gate at both tau: the root never splits, so the depth-1 corner is never reached.
  - X2 z at tau 0.6: 2.81 / 1.97; at tau 1: 0.88 / 1.67 (z_crit 3.71). X1 replays at the same labels read z 7.6-13.8.
  - The r131 root is not the closed split in any tree: closed enrichment of root states x0.87-1.13. It is a near-tie: root NMI
    between seeds .43 on X, .04 under fission.
  - Its real signal falls 2-4x at tau 1, because the noisier X1 labels change the partition itself.
  - Honest X2 rebuilds of the root maps read .83-.87 closed: no corner.
- TRAP: groups selected on X and re-read on X2 keep their imprint, because X2 = S + N - Z/tau shares N with X. Only X1
  selection makes an X2 read honest.

**Why:** the user asked to replace state_split's hand-set guards with a principled test (decision 6) without losing the RNase P
closed state.
**How to apply:** don't adopt fission in the split test as it stands. The guards do two jobs: fission can replace the
self-selection job, while the locality guard's "states differ locally" job stays. The r131 corner sits under a weak root
near-tie that no honest test confirms at 44k particles. Routes that keep full-data labels (band split, gold-standard halves) or
cross-fitting are untested, and the closed lineage could be judged at 692k instead. Related:
[[extremes-fission-verdict-2026-11-07]], [[cross-half-test-needs-band-above-assignment-2026-11-07]],
[[data-fission-merge-test-verdict-2026-11-07]], [[rnasep-closed-state-must-be-retained]].
