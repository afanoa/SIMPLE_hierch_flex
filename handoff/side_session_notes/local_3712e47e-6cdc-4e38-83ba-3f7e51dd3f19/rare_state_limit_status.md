# Rare-state detection limit under the subspace site test (session tree-rare-state-limit, worktree tree-rare-state-limit-d3f0e0, branch claude/tree-rare-state-limit-d3f0e0 off claude/tree-dist 26c350abb)

Work dir on batch2: /mnt/beegfs/elmlund/afan/hetero_tree_20261028/detlim2 (scripts/, runs/<arm>_s<seed>/<fixture>, traces/, logs/, build_vd).
Fixtures: detlim/fix/P_L03 P_L04 P_L05 P_L065 P_L08 P_L10 (conf18 closed arm at 3 / 4 / 5 / 6.5 / 8 / 10 % of 60k from pool2) + P_N00 (pool2 null).
Arms: old = dist_sdisc/build12 (9005a8c77, halves test); new = detlim2/build_vd (claude/tree-dist 26c350abb, variant D subspace test; built on x6342
like build12, job 54827647). Keys: recipe18d gate keys, minpop 260, maxdepth 0, nodepar 1, SIMPLE_SPLIT_CONTRAST=maps, SIMPLE_TREE_KEEPMAPS=yes.

## 2026-11-02 ~01:10 (batch2 time)
- Launched seed 1 of both arms, 4 trees at a time (sbatch singleton slots dl2_slotA-D): jobs 54827648-54827661. Seed 2 at the two
  shares nearest each arm's limit once seed 1 reads out. No Fortran change; I do not touch state_mergeup, site_rank, discover_sites.
- Second part (coordinator's addition): oversplit cost (forced split one level down, both arms) and a needed/unneeded characterization of
  every split on the saved gate trees (math_reworks/gates/base_* and sub4_*, maps kept), read-only on those directories.

## 2026-11-02 ~02:00 (batch2 time): first four roots (5 % and 6.5 % of 60k, both arms)
- Both arms IDENTICAL below the site layer (same sites, fallback seeding, stage A/B, re-test PASS z 5.9-6.2, guard merges L/L_ref 0.27-0.33):
  one leaf at both shares in both arms. The subspace test accepts endpoints the halves test rejected (6.5 %: site 3 poles at z 7.6 and
  15.0/4.1) but the pole rule rejects them and none is the rare state.
- Losing stage, with numbers (6.5 %, arm site 5 = 37 % of the rare-difference energy in 0.8 % of the box): the EM endpoint tracks the rare
  state (posterior AUC 0.75-0.78, ~5,000 equivalent pure members vs 3,900 true) but its cross-half deviation correlates at r 0.13-0.15
  while 3 Fisher SE at the site's n_eff 116.5 is r 0.275 (z 1.4-1.6). Old test: margin 2.1 SE. Relevant to the site_rank owner: the
  acceptance bar at small-site n_eff, and starts (fflex z1/z3 tails rank the rare particles at AUC ~0.6) are the levers, not the test form.
- Oversplit characterization written (note section 3.2): diffuse nuisance splits carry the largest power (micnull 48 % of the signal power,
  44x the noise) and only localization catches them; among delivered splits the one nuisance split (qtoy2 r4, P/N 0.8) sits where the
  rare-state detections at the limit sit (P/N 1.0-1.7). Proposal: guard + P/N >= 2 floor, margin 1.35x on the r131 closed corner.

## 2026-11-02 ~02:30: 8 % and 10 % roots
- Both arms FIND the rare state at 8 % (19.6) and the new arm at 10 % (24.5); lost at 4 / 5 / 6.5 % in both. Limit unchanged at
  n/sqrt(N) 16-20 (= 1.9-2.4 % at 692k). Seed 2 at 6.5 / 8 % and the forced-split halves (P_L08_h*, P_L10_h*, 30k, minpop 130) launched.
- FOR THE site_rank OWNER (variant D, detlim2/runs/new_s1/P_L08/1_hetero_tree/node_r/site): at 8 % the subspace test accepts the arm
  endpoints z1lo / z3lo / z3hi (r 0.27-0.29, z 4.0-4.3 at n_eff 217) but groups {z1lo, z3hi} (both N* 40,000) vs {z3lo} (N* 10,000) and
  the seed assignment gives pole 1 717 seeds (seed z 1.36 -> rejected) and pole 2 59,145 seeds (share 0.99 -> rejected): "no pole kept".
  The halves build on the same starts keeps one pole {z3lo, z3hi} of 18,600 seeds (rare P 0.21 R 0.81, margin 3.42) and seeds from it.
  A 10,000-N* endpoint class drawing 59,145 seeds at posterior >= 0.5 looks like a sign / complement issue in the pole seed assignment.
  The state is still delivered via the fflex fallback (stage A 0.30/0.94), so the delivered-leaf limit is the same in both arms.

## 2026-11-02 ~03:50 (batch2 time): DONE. Note doc/implementation_notes/tree_rare_state_limit_2026_11.md in the worktree tree-rare-state-limit-d3f0e0 (branch claude/tree-rare-state-limit-d3f0e0, UNCOMMITTED: the maintainer commits)
- 26 roots, 218 core-h, detlim2/runs 235 GB kept maps (nothing deleted). Memory rare-state-detection-limit-2026-10-29 updated.
- VERDICT: the subspace site test does not move the root limit. Both arms identical as delivered leaves on all 26 runs: lost at
  n/sqrt(N) <= 15.9, found at >= 17.1 -> 1.9-2.1 % at 692k (the 2026-10-29 bracket 11-15.5 closes to 16-17).
- Losing stage: the arm site's endpoint acceptance (r 0.13-0.15 vs 0.275 at n_eff 116; both tests), then the guard merges the fallback
  split (0.33-0.48 < 0.5). Lever proposed (not implemented): climb-by-scoring at the site (audit 6.5) to double endpoint purity.
- Oversplit: sqrt(2) per unnecessary level verified (8 % found at 19.6, both halves lost at 13.7/14.0; 10 % found at 24.5, halves found at
  17.1/17.5). Split characterization: 145 child pairs scored (detlim2/splits/*.tsv, analysis_saved_v03.txt); the guard catches the diffuse
  nuisance splits (largest power of all); proposed guard + P/N >= 2 floor; the OLE root split is real by the maps.
- For the site_rank owner: new_s1/P_L08 site 1 pole-seed assignment (10,000-N* class -> 59,145 seeds) and the old test seeding from
  non-rare poles at 8 % s2 / 10 % s1 (detlim2/traces/old_s*/P_L08.txt, P_L10.txt).
