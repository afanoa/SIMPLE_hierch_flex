---
name: site-pole-closed-state-replay-2026-10-27
description: "route A (site_states, branch claude/site-states): r1 replay extra closed seed survives (K*=3, 0.81 at node level); fixtures: qtoy closed rescued but diluted, 0 extra states on both nulls, r1211 unchanged; full692 3+3 pair 54435946/47 running"
metadata:
  type: project
---

2026-10-27, factorial-layer session (branch claude/site-states on 8c1b66f31, uncommitted; key site_states; build
site_states/build_sst2, fast gate 15/15).
- Test 0 (full692 V3 r1 on its own build s47d, closed pole as seed state 3): stage A 39k -> 64k, stage B 65k, state_split K*=3,
  child 69,700 rows (79 % from the extra state). Closed-arm / open-arm at 10 A: extra state 0.81 / 0.50; r1's own stage-B states
  0.62/0.76 and 0.57/0.78; pole top-2.5k 0.84/0.52; the tree's 8.7k closed leaf r12112_g1_2 0.96/0.65. ~32 core-h, 41 min.
- Test 0b (seed = unrefined home ranking, 25k): same, K*=3, 0.79/0.45.
- Gate (>= 0.9 at >= 5k) not met at NODE level: a 65k state is a closed-enriched mixture; the original tree reached 0.96 three
  levels below r1. Next: maxdepth-3 subtrees from the extra-state child vs r1's original closed-holding child.
- Route B prototype (cdisc toy, leaves at the 0.08 closed base rate): local E-step seeded by the closed pole -> closed recall
  0.31 -> 0.85 at precision 0.38 (4.8x base); leaves-only control 0.08. 14 core-h.
Related: [[port-tree-site-layer-recovers-closed-2026-10-22]], [[focus-three-threads-2026-10-27]].
- ROUTE A (site_states v2, build_sst3) on qtoy2 seed 1 (where the plain tree loses closed): at the ROOT it adds 3 extra states,
  all REAL (truth): r_s2 5,040 -> child 6,288 rows 93 % A-domain/holo (base 11 %); r_s3 15,305 = open end of the swing ->
  children s1 -4.1 / -7.2 (node -1.4); r_s4 8,887 closed share 0.196 -> child 17,260 rows holding 72 % of closed (share 0.149,
  base 0.060). K*=5. Their higher snrq is a top-of-posterior effect, not quality. v2 guard: extra states < 1/2 node, driving
  seeds >= minpop (v1 bug: 5 poles drained the driving seeds). Leaves + nulls pending.
- Build-vs-chance (qtoy2, seeds 1-3): 8c1b66f31 keeps closed on 2/3, seeds 2/3 identical to d3513ff70; seed 1 differs at r1's
  fragile pole ranking (244-row pole) -> chance, not a regression.
- USER 2026-10-27: "Wait for the r1211 iteration results (1+1 / 2+2 / 3+3), then launch the full 692k route-A tree with whichever
  iteration count wins." (One-pass root splits the closed rows ~60/40 in the replay.) 1+1 = amp_mstep_v1/trees/r1211_base_s1;
  2+2 = r1211_it22x_s1 (job 54409782); 3+3 = start_point_v1 job 54409886.
- RESULTS 2026-10-27 evening (plan note doc/implementation_notes/site_layer_tree_integration_plan_2026_10_27.md in /Users/atifao/SIMPLE):
  - qtoy route A leaves: closed leaf 17,260 at P .149 / R .715 (F1 .25; same-seed base .16, base seeds 2/3 .50/.42); structure MI
    .96 vs .59; the closed child cannot split again at maxdepth 1. r1211 route A = base (no non-driving pole anywhere).
  - NULLS: qtoynull route A adds 0 extra states at all 7 nodes, leaves identical to base (largest 55 %: the BASE fails the 80 %
    bar via quality splits). OLE null route A (maxdepth 2): root 0 (sole candidate = the seeding site itself, rho 1.0), r2 0.
  - r1-child subtrees (maxdepth 3, 1+1): control (orig child 2) best leaf 0.93 closed-arm at 11.6k; route-A t0b child gives
    0.97 at 4.9k and 0.94 at 6.8k but holds few reference closed rows (purity .10-.16 vs .26) -> different particles.
  - full692 3+3 pair (build_sst3_it3; route A 54435946, base 54435947; site_states/full692/): root route A K*=4 (both 223k/122k
    extras became children; not quality axes by corr AUC .52/.46), base K*=2; root wall 105 vs 93 min; depth-1 extras absorbed.
- 2026-10-28 MOVED: beegfs rnasep_klem/site_states -> /mnt/beegfs/elmlund/afan/_to_delete_20261028/rnasep_klem__site_states (user will
  delete it; beegfs Default pool full). Result text copied to /scratch/cluster_scratch/atifao/site_states/beegfs_results_20261027/,
  scripts to /scratch/cluster_scratch/atifao/tools_backup_20261028/. New work: /scratch/cluster_scratch/atifao/site_states/ (build_sst7_it3crop_emb =
  route A v2 + validated_crop.diff + embedding fix, 3+3; full692 pair 54457751 route A / 54457752 base, nodepar 8). Never write to beegfs (SUPERSEDED 2026-10-28 16:3x: user said go back to beegfs, see cluster-scratch-user-quota-exceeded-2026-10-28).
- 2026-10-28 CLOSED NODE (sst7 pair, scripts/closed_nodes.py; ref = old 8,720-row closed leaf; maps at 10 A): both trees first separate
  closed at r11 (depth 2). Base r11: group 79,392 holds 95 % of ref, best map 0.82/0.74 (r112 0.87/0.68). ROUTE A r11 (K*=4): group
  45,796 holds 68 % at purity .083 and stage-B state 1 reads 0.99/0.57 = closed corner at depth 2 (a driving-site state; the root's
  extras concentrated closed into r1). Root extras r_s4/r_s5 -> children r5/r6 with corr AUC .28/.25 = quality-like: P4 needed on RNase P.
  Route A disk 2x base (root-sized node projects with site_layer=yes).
- 2026-10-28 BUILD-7 BASE DONE: 213 nodes, 132 leaves, nodes 791 min (13.2 h), labels 822 min; per-depth level walls 134/53/56/96/
  150/227/240/202 min at nodepar 8 (deep node COUNT sets the wall). recipe9 (54510093) deep nodes 7-12 min vs ~20 -> ~8-9 h projected.
  TRAP: the end-of-tree plain/ and landscape/ steps write hetero_step.* into every node folder, so span readers must skip them
  (scripts/level_stats3.sh); reading all step files gave 361-494 min "spans" and re-fired the node-finish monitor.
- 2026-10-28 BASE LEAF SCORECARD (scripts/leaf_scorecard.py, scratch): build-7 base has NO closed leaf: best F1 0.093, 0 closed-corner
  leaves (max 0.86/0.63), ref rows in 128/132 leaves, best 6-leaf union F1 0.158 -> the K*=2 cascade below r11 shreds closed (fails
  the hard gate). Route A leaf scorecard (same script) decides whether r11 state 1 (0.99/0.57 at node level) survives to a leaf.
- 2026-10-28 BUILD-7 BASE LEAF SCORECARD (job 54540166, site_states/leaf_scorecard_base_20261028_1438): FAILS the closed gate.
  Best closed-leaf F1 0.093 (r1211211_g1_2, 9,662 rows, arm 0.82/0.70); 0 closed-corner leaves (best closed-arm 0.86); the 8,720 ref
  closed rows sit in 128 of 132 leaves; best purity-ranked union F1 0.158. The K*=2 cascade below r11 scatters closed. Route A leaf
  scorecard decides; base-like recipe9/recipe11 trees must be scored with the same script before any timing win counts.
- 2026-10-28 ~17:30 FAST-LOOP FIXTURES (beegfs hetero_tree_20261028/fixtures): r112 (base lineage, 79k) has NO closed corner in any
  build (build-7 K*3, recipe9 strongest no K*2 / yes K*3, best map 0.80/0.63, groups at node purity) -> strongest_pole does NOT fix
  closed; d3513ff70 doesn't split r112 (junk screen fires despite junk_screen=no). Route A r111 fixture: stage-B s02 25,270 rows reads
  1.00/0.59 (closed corner) but state_split delivers it mixed with 15k of s01 (40,257 / 5,539) -> delivery fix being replayed
  (adaptive). Route A ported to recipe12b as recipe13 (keys site_layer=yes site_states=yes). Fixture roots run the ROOT pipeline,
  not the oneflex child path: compare fixtures among themselves only.
- 2026-10-28 ~18:40 SUBTREE REFERENCES (build-7, fixtures on beegfs): base r112 (site_layer=no): 13 leaves, 0 closed-corner, best map
  0.88/0.63, F1 .105 (full tree .121: fixture reproduces). ROUTE A r111 (site_layer=yes site_states=yes): 9 leaves, 2 CLOSED-CORNER
  LEAVES (r1111_s1 5,373 0.97/0.60; r1121_s2 3,976 0.95/0.55) + r1111_s2 6,436 0.99/0.61 -> ~15.8k closed rows over 3 leaves (split by
  the s02+s01 fusion; recipe14 choose should keep it whole). Disk 129 GB vs base 44 GB (root-sized node projects) -> row-local site
  layer needed before 692k (asked of optimization session, recipe15).
