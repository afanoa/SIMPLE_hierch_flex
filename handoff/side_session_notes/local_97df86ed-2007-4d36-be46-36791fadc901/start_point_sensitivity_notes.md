# Start-point sensitivity (session: start-point-sensitivity-cb0d4a)

Write-up: .claude/worktrees/start-point-sensitivity-cb0d4a/doc/implementation_notes/planned/hetero_tree_start_point_sensitivity_2026_10_27.md
Cluster: /mnt/beegfs/elmlund/afan/hetero_tree_regression_20261026/start_point_v1 (install/, runs/, starts/, scores/, tools/)
No edits to simple_hetero_tree.f90 / simple_state_mergeup.f90 / simple_flex_probe_fit_estep.f90.

## Native change (uncommitted)
refine3D key `state_anneal=n` (default 0 = unchanged): prob_neigh iterations j=1..n draw states from the table's
state posterior at beta=j/n. Built from d3513ff70 + patch, fast gate 0, unit_pftc_align2D3D 0; key off = byte-identical
labels to the clean build.

## Findings (A3+B3 node recipe, standalone refine3D fits)
- No existing knob softens STATE under geom search (prob_athres only shapes the projection draw; snhc is global).
- Annealing fails: r1211 latent start stays .54/.59 closed P/R; site-pole start drops .984 -> .92/.86/.80.
  r12111 cross-start ARI 0.002 -> 0.12-0.31; qtoy2/CAHRA no gain.
- Not non-convergence: B9 runs stop at 4 its (converged) with unchanged partitions. Distinct fixed points.
- Majority-vote consensus follows the restart mix (r1211: latent basin) and averages within basin (.846 vs .984).
- Candidate: pick the restart with the highest mean max state posterior (prob_state_evidence.txt) of its STAGE-B
  partition. Rank corr with truth r1211 0.88, CAHRA 0.93, qtoy2 0.89; right basin in every pool. Fails at stage A on
  r1211 (0.12) and on annealed restarts. qtoy2 quality-start adversary: no confidence advantage, refines to structure.
- qtoy2 and CAHRA barely start-sensitive; RNase P (r1211, r12111) and OLE r1111 are.
- Halves deviation-map gate saturated (0.93-0.998) on real fixtures; null floor pending.
- qtoy NULL (one conformer): 7/8 random K2 restarts give the same quality split (snrq 1.25 vs 0.82, NMI .19-.21),
  confidence 0.86-0.90 (above real data), halves gate 0.97; published state_split ACCEPTS it (pair stat 0.66 vs
  threshold 0.096, K*=2). state_split accepted a split on every node tested. Sent to adaptive-split-merge.
- OLE null not run (8 restarts cancelled: ~1.8 h each, over their 1 h 40 min limit).
- CLOSED 2026-10-27 (wind-down). Cost: r1211 3.1, qtoy2 2.6, CAHRA ~3, r12111 1.8, OLE r111 3.6 core-h per A3+B3 fit;
  study total 545 core-h (406 completed, 139 cancelled).

## Reopened 2026-10-27 evening: iteration budget (ITA = ITB = 3, uncommitted, simple_hetero_tree.f90 line 35)
Tree A/B, matched d3513ff70 builds (start_point_v1/install = 1+1, it3/install = 3+3), seed 1, trees in start_point_v1/trees.
Root node, interim (leaves pending; watcher on batch2):
- r1211 subtree: closed group P/R .345/.958 (1+1) -> .474/.914 (3+3); closed stage-B state .337/.948 -> .500/.937.
- qtoy2: unchanged (stage-B NMI conformer .422 vs .426).
- qtoy null 3+3: root split 54/46 % by quality (snrq .91 vs 1.24), already below the 80 % bar; stage-B states not by quality.
r1211 leaves (site_layer=no oneflex=yes minpop=3000 maxdepth=2; no site splits the root, fflex groups seed):
- 1+1 (5 h 01 min wall): closed scattered: r12_s1 4215 P .691 R .334; r12_s2 6147 .472/.333; r11_s2 8156 .213/.200.
- 3+3: closed kept in ONE leaf r1_g1_2 16813 P .474 R .914; node r1's stage-B states split another axis (.354/.591 closed),
  state_split K*=1 (board pair stat 0.98 > 0.193, but cross-half EM re-test 0.08 < 0.144 in band 16..23).
- Reference full692 tree made r1211 a 3-state node from 3 site-rank poles (3117/25004/900 seeds, 3+3 passes); state 1 = closed.
  Neither test tree isolates closed: the start decides; 3 passes necessary, not sufficient. sacct down (head01p) -> core-h pending.
Parent cancelled r1211_it3_s1, qtoy2_it1_s1, qtoynull_it3_s1 (user: cut queued work); qtoy2_it3_s1 still running (watcher on batch2).
qtoy2 3+3 vs qtoy_base_s1..3 (maxdepth=1: leaves = depth-1 stage-B states; our 1+1 tree reproduces base_s1 exactly):
- base: sizes 6.2k/9.3-12.1k/18.5k/23.2-25.9k, NMI conformer .513-.514, A .297-.338, closed best P/R .272-.353/.914-.916
- 3+3 depth-1 states: 6.1k/10.6k/19.6k/23.7k, conformer .508, A .309, closed .305/.897 -> within seed spread; no regression.
- Side finding: below maxdepth children get the split step's regrouped groups (move 13-20 % of a node), which follow junk
  (depth-2 groups NMI conformer .369/.412, junk .18 vs states .51/.04). Owner: split test (adaptive-split-merge).
CLOSED again 2026-10-27 night: parent cancelled qtoy2_it3_s1 too (user: cut cluster use). No watchers, no runs left.
Write-up sections 5.6 and 7 hold the final iteration-budget numbers; tree core-h not counted (sacct down with head01p).
