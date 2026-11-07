---
name: cls-expansion-class-path-defect
description: cls_expansion class=N runs calc_pspec and the greedy refine2D round on the WHOLE project (805k on NOT, 15+ min for one class) and the subclass ids 1..nsub collide with non-member class ids; restrict to members before using class= again
metadata:
  type: project
---

Found 2026-10-07 on NOT (805k, box 420), `cls_expansion class=269 ncls=100 nthr=8` in memory: the split of
one class finished in 3.5 min, then the commander's sigma2 estimate (calc_pspec) and the greedy
`refine2D` round ran over the entire project, 15+ min and still going when killed. Worse, the split
writes subclass ids 1..nsub into the members' class field while non-members keep their original class
ids 1..300, so the round treats non-members with class <= nsub as subclass members.

**How to apply:** before any `class=` use, make exec_cls_expansion deactivate non-members (state=0 where
cluster /= class) in the run copy before calc_pspec/refine2D/restoration, or skip steps 2-3 for class=.
Also: 100 subclasses per class degenerates the divisive bisection (a few big leaves + dozens of 2-15
member leaves); a target-size rule is needed for large ncls. Master's Release 4 needs `fix_projfile`
on old projects (`nptcls_stk`), and old-format sigma2_state.bin files are rejected by refine2D.
Slim model (no cross-fit/nuisance/residual) = local branch cls-expansion-slim e39834e91, gate 14/14,
10180 0.889/0.666, 10076 0.825/0.793; NOT pushed. See [[flex-cls-split-state-2026-09-30]].
