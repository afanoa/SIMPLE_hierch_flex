---
name: tree-core-budget-tree-ncores-2026-11-02
description: "tree_ncores core budget built (claude/objective-meitner-e64b05): parts = budget/(nthr x nodes in flight) - 1, floor 15k particles per part; r131 A/B shows the masters grow 11.5 s per extra part (partial-file reads) so the 692k root's optimum is 36-58 parts; projected 692k gain only ~20 min (5.52 -> ~5.15 h), not the audit's 1.5 h: the root's nparts-bound steps are half master-serial at 12 parts already"
metadata:
  type: project
---

2026-11-02, session objective-meitner (worktree objective-meitner-e64b05, branch claude/objective-meitner-e64b05 off
claude/tree-lever-floor c5176bddc). Key `tree_ncores` (default 0 = today's fixed sizing): before each distributed step a node
gets budget / (nthr x nodes in flight) - 1 parts when more than one-per-20k (>= 4, <= nparts), capped at 128 and one part per
15,000 particles (PTCLS_PER_PART_MIN); core_room uses the budget; log line "stageB submitted on 46 parts (640-core budget, 1 in
flight)". Builds math_reworks/build_cb (floor 500, measurement) and build_cb2 (floor 15k, committed), fast gate 15/15; gates
cb_r131 / cb48_r131 / cb48_qtoy2 / cbf_r131 on norm; analysis scripts/hetero_tree_profile/step_rounds.py; note = audit 9.1.
- recipe17d root (114 min, 12 parts): part rounds 41.6 min, the masters' own serial time 41.7 min (flex 11.5, fflex 8.4,
  seeds 4.3), single jobs 31. Parts time linear in particles/part (stage B 39 + 16.7 s per 1k, stage A 51 + 8.1, fflex 51 + 3.3).
- r131 4 vs 48 parts: parts 4-8 s each at 924 ptcls (fixed ~5 s), 48 sbatch in 2 s, BUT the masters grow 11.5 s per extra part
  (stA0 +131 s, stage B +155, stage A +70: one partial file per part x state x half x iteration at 0.15-0.3 s). Fixture chain
  44 -> 50 min, 18 -> 52 core-h. Leaves r1 within 0.5 %, corner kept; r2's marginal split flipped (z 4.7 -> 3.1 vs 3.71) =
  the production split test's single-draw variance, not the budget.
- 692k projection: root 114 -> 94-98 min at 36-58 parts (the 15k floor gives 46 at 640 cores; 79 parts without the floor is
  worse, 96-105); depth 1 -1 to -2 min per node; deeper nodes keep 4-5 parts; wall 5.52 -> ~5.15 h, core-h +2-4 %.
- qtoy2 md1 at the budget (48 root parts, 11 per child): NMI .804 (bc .812, head .773); root split identical to 3 particles,
  r1 split differs (707 + 21,227 vs 11,783 + 10,154: a pole accepted at z 18-20 where bc saw z 11-12). TRAP: the part count
  reaches flex_pca beyond rounding: same-nparts arms are identical (head vs bc: 0.0000), r131's root latent identical at 4 vs
  48 parts, but qtoy2's root latent differs 3-10 % rms (corr .995-.9995) and children inherit more (r1 corr .93 -> <.3).
**How to apply:** launch the 692k with tree_ncores=640 (keys otherwise as f692_lf) for ~20 min; the root's remaining levers
are the masters' serial time (flex_pca master 11.5 min) and the per-part partial-file reads; the body's is nodes in flight
within the budget. Related: [[tree-692k-wall-profile-2026-11-02]], [[split-test-single-draw-variance-2026-11-02]],
[[tree-runtime-goal-5-7h-at-700k]], [[hetero-tree-math-reworks-setup-2026-11-01]].
