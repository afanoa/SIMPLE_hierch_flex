---
name: two-run-consensus-verdict-2026-11-09
description: "VERDICT: 2-run consensus of replicate splits is a confidence filter: purer only when both runs made the same choice (fixtures +.002-.005 purity; RNase P root x3.28, r1 x4.29) at recall -.07-.21, loses states when choices differ; null-safe only after the split test (stage-B consensus beats random on all 4 nulls); IgM 10 % null twin splits its root at SIMPLE_SEED 2 (b11)"
metadata:
  node_type: memory
  type: project
  originSessionId: 6f195b76-225b-466c-a432-7fc8360148e8
  modified: 2026-11-09T09:29:54.722Z
---

Overnight panel read (coordinator, 2026-11-09). Notes: `hetero_tree_20261028/handoff/ole-panel/CONSENSUS.md`; `ole_panel/scripts/consensus.py`. Cost 110 core-h, against ~100 approved (the IgM replicates ran long).

- **Same choice** (same K / same site):
  - qtoy2 s1/s3 and the 14 same-K qtoyflex pairs: purity +.002-.005, NMI +.01-.02, recall -.005-.02.
  - RNase P root replays 2/3: x3.28 vs x2.91/x2.71, recall .52 vs .62/.59.
  - r1 two hits: x4.29 vs x2.00/x3.59, recall .62 vs .70/.83.
- **Different choice:**
  - qtoy2 (5 vs 4 children): purity .71 vs .61-.67, but one state is lost (recall .002).
  - qtoyflex (3 vs 2): purity .098 vs .14-.18.
  - r1 hit with miss: the closed group is lost (x0.22).
  - IgM fix (2 vs 3 children, ARI .08): best-arm purity .607 -> .639 with half the rows dropped.
- **Third-run agreement** rises on the consensus rows (OLE .74-.92), but that is selection: on those rows each run equals the consensus.
- **Nulls:** delivered splits are one leaf, so the consensus is one group. Stage-B states agree beyond random pairs (consensus .62 / .71 / .81 / .60 vs floors .50-.51, for qtoynull / olenull / micnull / IgM null).
- **FINDING:** the IgM 10 % null twin on b11 keeps one leaf at seed 1 but splits its root at SIMPLE_SEED 2 (50,501 / 49,499, NMI 0 on every arm): stock null-safety is seed-dependent there.

**Why:** consensus was proposed as a purity lever ([[methods-must-generalize-across-datasets-2026-11-08]]).
**How to apply:** use consensus only on the split test's children, and only once the discrete choice is reproducible (see [[near-tie-rule-verdict-2026-11-09]]: K 6/6/6 at the RNase P root with four XD starts). Check null gates at more than one seed ([[null-acceptance-majority-leaf-2026-10-27]]).
