---
name: tree-distribution-build6-2026-11-01
description: "hetero_tree distribution builds 6-12 (pushed as fork/hetero-tree 9005a8c77): NOT r2 fixture 6.5 h -> 1.9 h, peak 544 -> 256, core-h 353 -> 121, md5-identical"
metadata:
  node_type: memory
  type: project
  originSessionId: 15979298-7fb4-455a-9eb0-b32cccd74b25
  modified: 2026-11-01T23:34:14.647Z
---

2026-11-01 (session funny-pike, branch claude/tree-dist, PUSHED as fork/hetero-tree 21fac34db..9005a8c77; builds in hetero_tree_20261028/dist_sdisc):
- build6: site discovery as a worker job (S_SDISC) + 6-10x faster discover_sites; site layer gets nodepar x nparts x nthr cores.
- build7: node-level site batch shares that 256-core budget (sb_fill/busy_update).
- build8: 1 rank thread per 5,000 ptcls (8 is the sweet spot), site layer reuses node ranks: 2.53 h, peak 264, 165 core-h.
- build9: site_em sums all-particle maps once, class 1 = all - class 2 (exact): 2.26 h, 143 core-h, but peak 376.
- build11: pre_rec/refine_rec class maps gated on the core budget (they were never counted against nodepar).
- build12 (= 9005a8c77): junk_screen identity/smear/member_gap cut (report-only): 1.92 h, peak 256, 121 core-h.
Fixture labels/merged/sites/unplaced md5-identical build6..12; RNase P closed_r11 on build9 md5-identical to production (1h07 vs 1h51).
Report: side-session-notes/local_5bd42cba.../tree_distribution_report.md.

**Why:** user wants 5-6 h NOT trees on ~160 cores without changing results; said 592 peak cores is too many.
**How to apply:** build from 9005a8c77 or later (dist_sdisc/build12/install). TRAP: pad2() names sites >= 100 site_**.mrc (pre-existing, unfixed).
Related: [[cluster-use-cut-nodepar4-2026-10-27]], [[tree-runtime-goal-5-7h-at-700k]].
