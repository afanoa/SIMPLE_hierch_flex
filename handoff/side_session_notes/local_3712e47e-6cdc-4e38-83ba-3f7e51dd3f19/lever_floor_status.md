# Rare-state lever, seed anomaly, effect floor, merge-up inputs, pose statistics (session tree-lever-floor, worktree tree-lever-floor-bccc29, branch claude/tree-lever-floor-bccc29 off claude/tree-dist c7ac3b0bd)

Work dir on batch2: /mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/lever/ (sites/, splits*/, runs/bc_s{1,2}/, traces/,
site_score.py, launch_*.sh); builds math_reworks/build_lever (climb + complement, rejected), build_head (c7ac3b0bd baseline),
build_bc ((B)+(C)), build_de / build_lf ((D)-(i)+(E); build_lf = branch HEAD code, fast gate 15/15); gate arms gates/head_* and
gates/bc_* (+ _r131_s1), readouts gates/{head,bc}_readout.txt, gates/ole_lf.txt. Note on the branch: hetero_tree_reworks_2026_11.md
sections 9-12. Nothing pushed.

## 2026-11-02 (afternoon, batch2 time)
- (B) committed cb3196f0a: an endpoint whose site class grew past half the node is no pole (claims no seeds); a rejected majority
  rest stays the complement; the pole rule ignores a rest without endpoints. Cause: seed_poles' posterior >= 0.5 maxloc let the
  took-off class (share 0.83, posterior ~0.9 everywhere) take the real pole's seeds (717 / 59,145 at the 8 % root). Gates: nulls 1
  leaf, r131 both seeds identical trees to HEAD (corner kept), qtoy2 root moves <= 200 particles; 8 % root pole now 1,863 / 1,780
  seeds at purity ~0.75 but its own seed z 2.2-2.3 < 3 -> fallback as before; 6.5 % lost both seeds -> LIMIT UNCHANGED (16-17).
- (A) NO-GO, not committed: the climb (hard E-step from the top-N*) lifts no 6.5 % endpoint (best z 2.2 vs 3 over 1-3 passes, both
  seeds), turns the pure 8-10 % endpoints into anti-poles (AUC .89 -> .15) and gives the null roots kept poles (27k / 26k seeds).
  The complement form of (B) flips the r131 seed-2 root pole and was rejected too. Table: note 9.2; memory site-climb-lever-no-go.
- (C) committed 577ee25fc: effect floor after the guard, P/N < 2 merges (EFF_MIN 2; pair_maps = pair_locality's body). Fortran =
  Python port to 1-10 % on the same nodes, same verdicts. bc arm: fires only at qtoy2 r3 and r5 (both within-conformer, P/N 0.22):
  10 -> 8 leaves, NMI .773 -> .812; r131 untouched both seeds (corner 2.6-2.7, the 0.93 r2 substate 2.08 kept); nulls untouched;
  OLE depth-1 kept at 5.3-7.9 (landscape comparison pending at writing: gates/ole_lf.txt).
- (D) fe02cffb8: equiv_tree_inputs reads split_stats.txt only, node_effect deleted. Full deletion of test/decide/eqn_stats needs
  the merge-up's ugrp from decide_once (scope in note section 11, ~1,100 lines); not done.
- (E) d062ba9f2: report-only '# pose ...' lines in split_stats.txt (node flip30, per pair dflip / z within the split test's strata,
  flag dflip >= .04 & z >= 5). Reproduces the study: qtoy2 .0009, OLE root .021, micnull .273; only micnull's refused root pair flags
  (.071, z 20.6). Note section 12.
- Cost: builds ~12 core-h, site replays ~3, split replays ~2, roots 32, head arm 115, bc arm 112 (OLE included).
