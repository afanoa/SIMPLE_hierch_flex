---
name: integration-plan-2026-10-28-builds
description: "Current production-candidate build = recipe7 on /scratch/cluster_scratch/atifao/recipe7 (fixm split test + 3+3 + validated crops + embedding fix + map cleanup v3 keeping node masks); recipe1-6 lineage and what each test showed"
metadata:
  type: project
---

Builds (all on /scratch/cluster_scratch/atifao, x6342, fast gate 15/15; source = ~/SIMPLE_recipe1_src frozen fixm split-test code):
- recipe1: + ITA=ITB=3. recipe3: + intopt.diff (pre-fix crop patch: sets box_crop before prg -> distributed parts die "prg= or
  test= must come first" at node fflex; fast gate does NOT cover it). recipe4: + validated_crop.diff (intB->intBopt2) + embedding fix.
- recipe5: + prune_maps.diff (partials after each step; non-leaf .mrc under node_* at tree end). recipe6: + reprojection models,
  dist_part_neigh. recipe7: + keeps node mask.mrc (needed to re-run merge-up). qtoy2 maxdepth-1 test: leaves identical across
  recipe5/6/7, 12/12 leaf maps, tree ends 2.9 GB (5.6 GB recipe5). Env SIMPLE_TREE_KEEPMAPS=yes disables cleanup.
- recipe4 nulls at 3+3: OLE null 1 leaf; qtoy null FAILS (76 %, rho 0.044 > 0.03 floor; floor call flips with shuffle seed).
**How to apply:** launch new trees from recipe7/{tree,null_tree,prune_test}.sbatch (scratch run dirs, SIMPLE_SPLIT_CONTRAST=maps
SIMPLE_SPLIT_RHOMIN=0.03); the null guard and the Bonferroni count (node-local recovers deep RNase P splits) are still open
([[split-test-accepts-quality-split-on-null-2026-10-27]], [[test-runs-moved-to-delete-folder-2026-10-28]]).
- 2026-10-29 INTEGRATION RUN approved by the user: 692k RNase P on recipe17d (recipe13 route A + choose + guard + kway_maxdepth=2 +
  row-local site layer + rank reuse/prune + STATES_DELIVERED marker; maxdepth 4; maps + 0.03; strongest_pole + equiv merge-up
  defaults), nodepar 12-16 <= ~300 CPUs, beegfs hetero_tree_20261028/integration/f692_r17d_s1. Launched by the factorial session
  (owner of the recipe17 line). Scoring: closed lineage per node, leaf scorecard + landscape at STATES_DELIVERED.
- Integration root (job 54634275, 114 min, 133.6 core-h): K*=4 (252k/193k/147k/100k); choose delta 61,320 >> 3 SE 712; guard
  1.01-1.06; extras from r_s4/r_s5 poles. Closed refs: r1 56 % (purity 0.0194 = 1.5x node), r2 29 %, extras r3/r4 15 % at low purity.
- Integration r1 (44 min, 50 core-h): K*=3; closed group 80,775 rows, 3,317 refs (38 % of 8,720), purity 0.041 = route A r11's.
  r3/r4 binary, no closed signal. Watch r2 (29 % of refs): if it also isolates closed, leaf F1 is capped by the root's 4-way split.
- Integration r2 (48 min, 45 core-h): K*=6 (3 extras), closed refs spread over all six (best 982 refs at 0.43 arm) = minor side
  branch. Main lineage r1 g3: F1 .074 at depth 1 vs route A r11 .077 (P .041 both; R .38 vs .64), so the 4-way root costs no F1 yet.
- recipe17f junk-leaf rule on IgM: fired once (root junk group 27,743 rows, 69 genuine), -11 % core-h, -6 GB; NMI .429 vs 17b .447,
  span .45/.57 vs .53/.70 (likely the r1 seeding near-tie, not the rule). Mixed -> NOT in the default; RNase P runs junk_screen=no.
- Integration: side branch ends as leaf r2322 (7,927 rows, 743 refs, 0.78/0.68). Main closed r131 -> r1311 (0.90/0.57),
  r1312 (0.85/0.56); their depth-4 split is the last chance for a closed-corner leaf.
- INTEGRATION CLOSED GATE: r1311 -> leaf r1311_s2 12,010 rows, 816 refs, stage-B map 1.00/0.58 = CLOSED CORNER (P .068, F1 ~.08).
  r1312 (21,870, 1,824 refs, 0.85/0.56) refused (z 0.64 < 3.71) though its try had a 0.94/0.56 state -> one near-corner leaf.
  Closed lineage ends in 33.9k rows / 2,640 refs (30 %); route A r111 corner was 25.5k / 2,800 refs. Leaf scorecard pending.
- INTEGRATION FINAL (recipe17d, 2026-10-29): states in 5 h 31 min, 873 core-h, 10,130 jobs, 84 nodes / 67 leaves, disk 86 GB
  final (peak >= 176 GB). Leaf scorecard: 1 closed-corner leaf (base 0); best leaf F1 .119 r1312 (base .093); best union .134 over
  3 leaves (base .158 over 6). Landscape d4: vs route A W1 .017, coverage .94/.61; extrema closest to route A. Files under
  integration/final/. Next on the user's ask: consolidate recipe13-17d into hetero-tree commits.
- Integration embedding figures (user asked 2026-10-29): integration/final/plots/integ_umap_v2.py (+ umap_200k.txt), copies in
  ~/Downloads/rnasep_integ692k_r17d_{umap,embedding_z1z2}.png. Root embedding = ONE continuous cloud; closed refs only ~2-3x
  enriched toward the high-z2 edge; the tree reaches closed by per-node re-embedding, not the root axes. Use hexbin, not overplot.
- CLOSED-IN-TWO-LEAVES DIAGNOSIS (split-merge session, route9/integ_closed): corner r1311_s2 vs r1312 SNR_d 14.7/3.1/0.3
  (40-20/20-12/12-), L/L_ref 0.87, 58 % of the difference ~22 A from the closed arm -> r1312 = distinct closed sub-state, not
  dilution (open arm unchanged). r2322 differs from both mostly ~30 A away (more open arm): the closed refs under r2 are a
  different sub-state, not misrouted corner particles. Merge-up never equivalence-tested them (re-test always above threshold).
  User's eye decides whether closed sub-states are worth states ([[close-closed-substates-not-worth-a-split-2026-10-28]]).
