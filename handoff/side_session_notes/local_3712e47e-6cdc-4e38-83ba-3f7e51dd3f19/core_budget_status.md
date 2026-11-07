# Core budget for the tree's distributed steps (session objective-meitner, worktree objective-meitner-e64b05, branch claude/objective-meitner-e64b05 off claude/tree-lever-floor c5176bddc)

Work dir on batch2: math_reworks/ builds build_cb (floor 500, the measurement build) and build_cb2 (floor 15,000, the committed code),
sources src_cb / src_cb2, gates cb_r131 (fixed sizing), cb48_r131 and cb48_qtoy2 (tree_ncores=392 on build_cb), cbf_r131 (build_cb2,
tree_ncores=392), readouts gates/{cb,cb48}_readout.txt, latent comparison gates/cb48_latent.txt (compare_latent.py), accounting
dumps gates/cb*_r131_sacct.txt. Note: audit note section 9.1 on the branch; script scripts/hetero_tree_profile/step_rounds.py.

## 2026-11-02 (evening, batch2 time)
- Key `tree_ncores` (default 0 = today's sizing) committed dc72c27fc: parts = budget / (nthr x nodes in flight) - 1 when more than
  the fixed one-per-20k (>= 4, <= nparts), at most 128 and one per 15,000 particles; core_room shares the budget; each
  submission logs "stageB submitted on 46 parts (640-core budget, 1 in flight)". Fast gate 15/15 on both builds.
- Root anatomy (recipe17d accounting, step_rounds.py): of the root's 114 min, part rounds 41.6, the masters' own serial time
  41.7 (flex 11.5, fflex 8.4, seeds 4.3), single jobs 31. The parts' time is linear in particles/part (stage B 39 + 16.7 s per
  1k, stage A 51 + 8.1, fflex 51 + 3.3, ...).
- r131 A/B 4 vs 48 parts (norm, both arms waited on the cap at times; per-step table separates waits): parts 4-8 s at 924
  particles (fixed ~5 s per part, 48 sbatch in 2 s) but the masters grow 11.5 s per extra part (stA0 +131 s, stage B +155,
  stage A +70: a partial file per part x state x half x iteration at 0.15-0.3 s). Chain 44 -> 50 min, 18 -> 52 core-h.
  Leaves r1 within 0.5 %, corner 1.00/0.48 kept; r2's marginal split flipped (z 4.7 -> 3.1 vs 3.71).
- qtoy2 md1 at the budget (48 root parts, 11 per child, 50 core-h): NMI .804 (bc .812, head .773); root split identical to 3
  particles; r1 splits 707 + 21,227 vs bc's 11,783 + 10,154 (site pole accepted at z 18-20 vs z 11-12).
- Latent comparison: same-nparts arms identical (head vs bc 0.0000); r131 root identical at 4 vs 48 parts; qtoy2 root differs
  3-10 % rms (corr .995-.9995), children more (r1 .93 -> < .3). The part count reaches flex_pca beyond rounding (iterative
  stopping point or per-part sampling); marginal site/split verdicts flip downstream; the aggregate holds.
- Projection with the per-part master growth: 692k root optimum 36-58 parts (the 15k floor gives 46 at 640 cores): root 114 ->
  94-98 min; depth 1 -1 to -2 min per node; deeper nodes unchanged (12 in flight share the budget); 5.52 h -> ~5.15 h, core-h
  +2-4 %. The audit's 1.5 h hope was wrong: the root's nparts-bound steps are half master-serial already and the masters' per-part
  cost caps the useful part count. Remaining root levers: the masters' serial time (flex_pca master 11.5 min), the per-part
  partial-file reads; the body's: more nodes in flight within the budget.
- Pending at writing: cbf_r131 (committed code, 4 parts expected = today's tree) and the final note/script commit. Nothing pushed.
- Cost so far: builds ~4 core-h, cb 18.3, cb48 51.9, qtoy2 50.2, cbf ~18 (pending), readouts < 1.
