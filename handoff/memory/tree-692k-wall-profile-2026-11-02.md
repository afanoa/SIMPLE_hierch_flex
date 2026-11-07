---
name: tree-692k-wall-profile-2026-11-02
description: "where the 692k's 5 h 31 min go (recipe17d accounting, 873 core-h): root 1.91 h serial at ~60 of 642 cores (nparts 12 x nthr 8 cap), body 2.3 h at 230-310 cores, end phase 1.3 h = a 73-min serial merge-up; site layer 38 % of core-h; levers without science change: budget-aware nparts (~-1.5 h), draw-free merge-up, the site layer's per-site refine_rec"
metadata:
  type: project
---

2026-11-02, this session (audit note section 9, script scripts/hetero_tree_profile/profile_692k.py on claude/tree-dist; data
integration/readouts/r17d_sacct.txt on batch2):
- root phase 0-1.91 h (35 % of wall): a serial chain, every distributed step at nparts 12 x nthr 8 = 96 cores over 692k particles;
  wall minutes mean 3, flex 16, cheap 3, fflex 12, stA0 5, recB0 3, site layer 19, seeds 5, stage A 13, recB 2, stage B 24, split 6.
- body 1.91-4.20 h (41 %): 230-310 cores busy (264 while only depth 1's four nodes ran); end phase 4.20-5.52 h (24 %): plain maps
  3 min, landscape 2 min, merge-up 73 min as ONE 48-core job (58 core-h).
- core-h by step: stage B 176 (20 %), stage A 103, flex+fflex 88, stA0/recB0/recB/seeds 93, site layer 332 (38 %: node site_rank 41,
  per-site site_rank 169, per-site refine_rec 122), split 14, merge-up 58, plain 7.
**How to apply:** the wall levers that keep the science are scheduling: (1) a core budget instead of fixed nparts (root's five
nparts-bound steps = 70 of its 113 min at 96 cores -> ~1/4 at 60-80 parts; expected ~4 h total), (2) the merge-up on the draw-free
statistic (reworks note section 11), (3) what the subspace site test leaves of the site layer (read f692_lf_s1). nodepar alone
cannot speed the root. Related: [[tree-runtime-goal-5-7h-at-700k]], [[integration-692k-lf-launch-2026-11-02]], [[hetero-tree-algorithmic-audit-2026-11-01]].
