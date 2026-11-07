# Equivalence merge-up port and quality-matched null: hand-over (2026-10-28)

Full write-up: worktree adaptive-split-merge-b384e9, doc/implementation_notes/planned/hetero_tree_adaptive_split_merge.md
sections 10.10, 10.11 and 11. Diffs (copies here and on the cluster):
- equiv_port.diff: /scratch/cluster_scratch/atifao/equiv_port/equiv_port.diff (md5 ab7d7471...), one file
  simple_state_mergeup.f90, SIMPLE_MERGE_RULE=equiv, default off. Against recipe4 source; recipe7's state_mergeup is
  byte-identical (md5 7475ab88...), so it applies unchanged. Gate 15/15 on recipe4 (build1). recipe7 build: build7.
- qnull_v4.diff: /scratch/cluster_scratch/atifao/qnull_port/qnull_v4.diff (md5 af08c005...), stacked on equiv_port.diff,
  SIMPLE_SPLIT_QNULL=yes|core, default off. Gate passed (qnull_port/build1).

## Equivalence merge-up: native = offline on every tree
| tree | leaves | delta | distinct pairs | units native/offline | units today |
|---|---|---|---|---|---|
| qtoy recipe5 s1 | 4 | 0.00246 | 2 | 4/4 | 4 |
| qtoy2_bp base | 20 | 0.00044 | 0 | 19/19 | 16 |
| CAHRA sub_a | 8 | 0.00103 | 2 | 8/8 | 5 (today fuses both node-distinct pairs) |
| IgM m69M | 58 | 0.00009 | 17 | 53/53 | 33 |
| qtoy recipe7 s1 (masks, build7 = recipe7 + equiv, gate 15/15) | 4 | 0.00246 | 2 | 4/4 | 4 |
| qtoy null recipe4 | 2 | 0.00002 | 0 | 2/2 | 1 |
Gap to know: the cannot-link list holds stage-B leaves only (_s1/_s2, _pA/_pB); sibling child NODES that end K*=1
(r1_g1_2 / r2_g1_2) are not on it (offline rule has the same convention).

## Quality null (core score) 2x2, 3 seeds, maps column
closed pair kept 3/3 (z_q 2.1-2.25); qtoy2 kept 3/3; qtoy null refused 3/3; OLE null refused only 2/3 (z_q 2.57/0.78/1.46).
unset + core keeps the OLE null 3/3 -> maps still needed. Alone it cannot replace the 0.03 floor; with the floor it changes
no decision at the tree count. Node-local count (2 N_node/minpop - 1) keeps r2111112 2/3 and the closed pair even unset;
nulls unchanged. No default changed.

## 692k route-A root (build 7) re-test: no guard drops children 5/6
Diff: /scratch/cluster_scratch/atifao/qnull_port/qnull_v5_on_v4.diff (SIMPLE_SPLIT_GUARD=pairs, default off, gate passed): for
>= 3 EM groups the floor and the core quality null act per nearest-group pair; a failing pair merges (the old floor reads the
weakest STARTING pair, even one the re-test already merged, and collapses the whole node to K*=1; old qnull is 2-group only).
maps + floor 0.03 + QNULL=core, ntests 460, 3 seeds, with and without pairs (battery5, 12 runs, 7-10 min at 16 thr):
- route A: K* 6 in every run (205462/176770/125554/101725/43239/39642); weakest pair 0.87-0.88, z 23-24; nearest pairs
  0.87-0.93 vs core pseudo-split 0.26-0.48, z_q 9-13.5, rho 0.24-0.54. Child 5 vs child 2: 0.87 vs 0.34, rho 0.24;
  5 vs 6: 0.92 vs 0.26, rho 0.40. Child 3 kept.
- base: K* 3, rho 0.17-0.19, z_q 9-12.
- Under maps the EM keeps 5/6 at 70 % / 76 % of their tree members (the rest move mostly to children 1, 2); seeds change nothing.

## 692k deep K*=2 splits under the guards (battery6): no arm stops them
The build-7 trees ran NO guard: build_sst7_it3crop_emb = route-A branch on 8c1b66f31 without the split-test port (contrast
unset, threshold = max of 20 shuffles, no z test, no floor). Replay on qnull_port/build1 v5, seed 1, 48 runs, 8.3 core-h
(/scratch/cluster_scratch/atifao/qnull_port/battery6, qnull_task6.sh / qnull6.tsv / qnull_table6.py):
a = maps + rho 0.03 @ 460; b = a + QNULL=core + GUARD=pairs; c = b @ 2N/3000-1.
- a keeps all 12 nodes: z 4.1-20 vs 3.71, rho 0.036-0.20 (deep 7-13k nodes 0.07-0.15).
- b/c refuse only base r2321 (12,781: real 0.213 vs core pseudo-split 0.135, z_q 1.46/1.52). Others z_q 2.1-12.4.
  c = b (the local count only lowers z_crit).
- Closed pair m33 r1111221 kept in a/b/c (rho 0.047, z_q 2.18/1.90).
- Build-7 closed controls (factorial scorer): base r11 (z 16.97, rho 0.164, z_q 9.69), base r112 (12.69, 0.148, 7.88),
  route A r11 K*=4 (weakest 19.22; pairs: rho 0.15-0.19, z_q 11.3-13.0, K* 4 -> 4): kept in every arm.
- Gap: build-7 closed splits z_q 7.9-13 vs deep nodes 1.5-4.7, so a z_q ~5 bar would split them, but it refuses the m33
  pair, and the closed purification chain below r112 is untested.
- The deep splits have 1.5-3x the closed pair's rho and z_q at or above it -> any null guard that stops them drops closed
  first. Not null-like (qtoy null rho <= 0.029). Stopping needs a non-null rule (depth/minpop/budget).
Doc: hetero_tree_adaptive_split_merge.md section 10.12.

## Top-down z_q >= 5 stop on the build-7 692k base tree (battery7, 42 replays, 11.8 core-h)
Arm b settings; walk = qnull_walk7.py; state at 8 h tree wall (145 nodes listed, 96 finished).
- 27 pass (lowest r321 5.59, r1112 5.27); 18 refused -> leaves 11-52k (median 23k), 415,624 particles (60 %):
  r3111 4.94, r3121 4.70, r3112 4.55, r322 4.42, r122 4.34, r2112 4.16, r3122 4.14, r232 3.85, r3212 3.68, r231 3.51,
  r22112 3.42, r222 3.41, r3211 2.74, r2212 2.46; replay K*=1: r212, r1122, r1212, r22111.
- 12 unfinished under passing parents: r11111/2, r11121/2, r11211/2, r21111/2, r121111/2, r121121/2.
- Saved: 88 of 145 listed nodes lie under refused nodes (51 run, 37 still to run).
- Closed lineages (scorer rerun in battery7/closed_score_*): A r11 9.69 -> r112 7.88 -> r1121 8.36; B r12 7.71 -> r121 7.66
  -> r1211 11.70 -> r12112 6.23. All pass (lowest 6.23 vs highest refused 4.94). No finished node at the closed corner yet
  (best 0.88/0.68 r1121, 0.86/0.72 r12112 child); the corner-reaching splits are among the 12 unfinished.
- m33 lineage: r1 7.96, r11 11.66, r111 10.98, r1111 10.91, r11112 6.04, r111122 2.41 with z 3.45 FAILING at 460, pair 2.18.
  The pair is a weak within-closed split (conf19/conf18) the operating point already misses at its parent, not an artefact.

## z_q >= 5 on trees with known answers (battery8, 41 replays, 4.7 core-h): FAILS
qnull_walk8.py; each tree at its own count; bar 1.67 (= arm b's own quality-null level) shown to isolate the bar.
- qtoy2 recipe7 (maxdepth 1, 4 leaves): nothing removed (z_q 17-23); all conformers unchanged.
- qtoy2_bp (m8c90t base, 20 leaves): bar stops r12 (z_q 2.08) above r121 (7.98) -> r1211 (12.82, closed-arm split):
  c18 closed best unit 0.70/0.93 (r12111) -> 0.78/0.27 (r12*); c22 0.27/0.62 -> 0.55/0.25. At 1.67 both kept.
  Holo c101 r22211 0.28/0.83 -> r22* 0.67/0.74 at BOTH bars (arm b keeps r22 at K*=1).
- OLE recipe4 vs native equiv (17 units): bar frontier 8 units; 9 equiv units vanish (r21* swallows 7 small states incl.
  r2112_s1/_s2; r122_s1/_s2; r11221/r1122; r1112 pair already at 1.67). All 4 node-distinct sibling pairs folded.
- Nulls: qtoy null root refused by arm b's own qnull (z_q 0.51); OLE null root K*=1. Bar opens nothing.
- Cause: z_q is not monotone down a lineage (weak parent over strong children). Top-down stops can't see them.
- Factorial base tree: r21111 7.01 pass, r21112 3.99 refused; closed-lineage nodes still running at 8.5 h.

## Default-on equivalence merge-up (recipe12) and the scratch move
- equiv_port_default_on.diff (md5 0869573b; one switch: l_equiv = .not. (SIMPLE_MERGE_RULE == 'off')), against recipe11
  source. Built + gated on beegfs (job 54543757, /mnt/beegfs/elmlund/afan/hetero_tree_20261028/equiv). OLE recipe4 tree:
  unset 17 units == offline == recipe4+equiv build; off == the tree's delivered merge-up (7 units). Verdict sent to
  tree-optim-budget, which applied the diff in recipe12.
- Scratch went over the user's quota at ~15:09. What this work still needs is copied and verified to
  /mnt/beegfs/elmlund/afan/hetero_tree_20261028/from_scratch/ (job 54543775, 25,252 files, 4.35 GB): qnull_port
  (build1, logs, scripts/diffs, battery logs + text), equiv_port (diffs, scripts, src_r11_defon, build1, build7),
  mergeup_inputs.py code. Battery s.simple copies were not kept (re-copy from the trees if a replay is needed).

## Split membership: EM groups dilute the closed stage-B state (route9, 95 replays, 28 core-h, beegfs route9/)
Why: the children are the cross-half EM's groups (shells 2..12, <= 20 A). On the route-A r111 fixture the EM moves 33 % of
particles off the stage-B start and lands orthogonal to the closed axis (a 5k slice of s01 vs the rest). The existing
SIMPLE_SPLIT_ROUTE modes map each stage-B state to the EM group holding most of it: both map to one group, so EM stands.
Fix tested (route9_on_v5.diff): when K* == number of stage-B states, the states themselves are the K* groups; modes
choose (data pick in the re-test band, stage-B beyond 3 SE), stageb (whole), spost (one cross-half E-step from stage B).
- s02 intact with choose = stageb on route-A r111 root (25,528), base fixture r1 (15,139), strongest-pole r1 (15,814);
  ctl fuses 15.5k / 7.1k / 8.2k of s01; spost leaves 3.3-5.8k; band_direct = EM. Base fixture r2: K*=1 at the operating
  point in every arm (z 2.93 < 3.71).
- Guards: nulls K*=1 all arms (routing never changes K*); qtoy2 r1211 c18 0.97/0.91 -> 0.90/0.98, junk split r2211 same
  partition; m33 pair K*=2 as its stage-B states. Closed reference rows: s02 child holds 73-76 % at higher purity (ctl's
  fused child 91-93 %). OLE cross-seed ARI unmeasurable here (seeds split r1/r2 on different axes; ~0 in every arm).
- Choice: choose with the identity mapping. Doc section 10.15.
- ADOPTED (coordinator): choose_default_on.diff (md5 a7fdd9cc, against recipe13 merge-up md5 9c6acef5; SIMPLE_SPLIT_ROUTE unset =
  choose, =em = EM groups). Verified job 54549369 (gate 0; r111 root unset == route9 choose, em == ctl). recipe14 =
  recipe13 + this diff (factorial session, single writer, build 54549505).

## Null vs real split: localization of the state-map difference (map_sep.py, beegfs route9/sep, seconds per node)
L = top-5 % voxel share of the envelope-matched, noise-subtracted residual between two stage-B state maps (40-10 A, half
maps, soft sphere); L_ref = the same for the consensus density.
- Nulls: qtoy2-matched root (recipe14 = recipe13 stage B) L 0.19 (ratio 0.3); OLE-matched root L 0.13 (0.2).
- Real: route-A r111 0.45 (0.7, weakest), m33 closed pair 0.71 (1.1), r1211 closed/open child 0.85-0.87, qtoy2 root pairs
  0.73-0.84, qtoy2_bp closed arm 0.91, junk split 0.82, OLE r1 0.87, all 34 build-7 closed-lineage pairs 0.74-0.96.
- Holds for top 2/5/10 % and sigma 2/3/4. D (attenuation-corrected correlation) and dB do NOT separate.
- Guard candidate: refuse/merge a split with L/L_ref < 0.5. Only 2 independent nulls; confirm on more null seeds first.
- Null seeds 2-4 (recipe14, production env, root only): all six roots L/L_ref 0.2; none split in the tree. Guard built:
  locality_guard_on_r14.diff (md5 1653fc1b; state_split, default on, SIMPLE_SPLIT_LOCALITY=off disables; reads ../stageB
  halves). Replays: null 0.27 merged, r111 0.66 / qtoy2 0.97-1.08 / r1211 1.25 / closed arm 1.26 / OLE r1 1.48 kept.
  Handed to the factorial session for recipe15.

## Consolidation replays (route9/agg, /repro, /stab): none consistent
- Bottom-up map equivalence (vs the tree's smallest kept effect): 0-1 merges per tree; +localisation 2-5 merges.
- Reproducibility merge (overlap over the other run's units vs re-split null): OLE ref 31->29 / 26->22, strong 0, IgM 1; null too tight.
- Stability cut: scale-dependent (OLE ref6 14 / 2 / 31 units for rho / SNR_d / sqrt rho); sqrt rho best on qtoy2 F1 (0.616 vs
  0.574), SNR_d loses IgM extremes; RNase P refs keep all leaves. Doc section 10.18.

## Held-out split test (route9/hold, 45 replays, ~25 core-h incl. build): the no-EM variant meets every bar
holdout_on_r17b.diff (md5 f44d1607, recipe17b mergeup 3d42a759 -> 463b54c0, fast gate 0): SIMPLE_SPLIT_HOLDOUT=R replaces the
re-test verdict. Each half's K group maps + consensus come from that half only (SIMPLE_SPLIT_HOLDOUT_ITERS=n own-half E-steps
first; 0 = stage-B board labels as they are), are docked to the other half's consensus, and score the other half's per-particle
log-likelihood gain (mixture - consensus) over klo..kb (40-20 A); swap; z vs R=10 stratified shuffled-label refits.
- Guard off: every null seed K*=1 in both holdout arms. Production splits qtoynull s1 and s3. Nulls' stage-B labels score
  BELOW shuffled labels (T < 0, z -24 to -376): own-half selected groups do not generalise.
- No EM (ITERS=0): all 7 positives pass.
  - m33 closed pair z 21.0 (production 3.66), K*=2; r2121212 z 6.42 (0.98) but the floor refuses it (rho 0.0296).
  - r1211 30.5, qtoy2 1718, closed arm 34.4, IgM 776, OLE r1 243.
- 3 EM its: the shuffled refits rediscover real structure (qtoy2 null 2.96, closed arm 8.8-12.1), and the closed arm is LOST (z 0.56).
- Per-particle z: m33 pair 7.6, r2121212 0.67, nulls -6 to -42. The shuffle sd is tiny, so the weakest positive relies on it.
- Cost/split step (16 thr): no EM 4.1-5.0x cross_em, state_split wall +7-63 %; 3 EM its 14-16x, 1.5-3.4x.
- Production guard on qtoynull_s1 skipped: recipe14 gates_env qtoynull_s1 node_r/stageB has no half maps.
- Handed to the factorial session (Assess and integrate the factorial site layer). Doc section 10.20.
- route9 is 71 GB (nulls/ 59 GB = six null-seed trees), over the ~50 GB budget; nothing moved. hold/_to_delete_20261029 holds
  a cancelled build's partial copy.

## Holdout v3: v1's per-map docking caused the recipe17e false negatives; r1211 r11 still fails
- **Cause:** v1 docked each group map (translation up to 3 px, ~15 Å at box 48) to the other half's consensus. Inside
  1.4-1.7 % focus masks this spoiled one map, so T = log π: qtoy2 r3 -0.54 = log 0.61; r11's refits bimodal 0 / log 0.41.
  It also inflated the closed arm, T 12.0 vs 0.83 without docking.
- **v3:** no docking. holdout_on_r17b_v3.diff (md5 fb230799, against 3d42a759); v1_to_v3.diff (04195910).
  Fast gate 0; patched mergeup a1bbe68e.
- **v3, no EM, 15 nodes + two 17e nodes:**
  - all 8 nulls K*=1 without the guard;
  - qtoy2 r3 z 342 (recovered), m33 pair 21.3, r2121212 10.6 (floor refuses), closed arm 86.7, other roots all pass.
- **r1211 r11 FAILS** in both arms: no EM z 1.47 vs 2.79, labels beat shuffled by 0.0035 nats/pp, pp z -0.21; EM z 0.50.
  - dsnr_band: SNR_d 5.7 at 40-20 Å, 86 % of the information. The maps differ, but particles do not prefer their state.
  - Production EM moves 27 % of its labels.
- **EM arm on v3** rejected: shuffled refits rediscover strong structure (qtoy2 r3 null 0.87, z 2.13 -> K*=1).
- The answer to the coordinator's Q1: the score is the mixture gain vs consensus, not under the particle's own label.
- Doc section 10.20 rewritten. Nothing new built beyond the in-place v3 test build (route9/hold/build).

## Holdout: two more 17e nodes and variants v4/v5. No variant meets every bar; production + floor + guard stays.
- **Truth (child_truth.py):**
  - qtoy2 r11 = stage-B states: real but impure. The 5,559 side is c1 88 %; the 7,223 side is c1 62 / c2 26 / c22 10 %.
  - qtoy2 r4, as v1/v3 split it: a false positive. Both sides are 98 % c101; junk 6.7 vs 13.9 %.
- **v3:** qtoy2 r11 z -23.7 (FN); qtoy2 r4 z 10.6 (FP; per-particle z 2.9).
- **v4** (noise fit per particle's best group): selection on noise; nulls s3/s4/olenull s2 pass. REJECTED.
- **v5** (one consensus noise fit shared by both scores; holdout_on_r17b_v5.diff, md5 62ae44e9):
  - qtoy2 r11 z -8.9 (shuffled partitions gain 0.083 vs the real 0.042);
  - qtoy2 r4 z 13.5 (FP); RNA r11 3.22 (pass);
  - qtoy null s3 z 10.9 on zero absolute gain; the guard merges it.
- **Best combination, v3 with both z and per-particle z >= z_crit:**
  - keeps the m33 pair; refuses r4 and all nulls;
  - loses qtoy2 r11 and RNA r11, which production keeps.
- route9/hold/build now holds v5. Scripts and tables in the worktree: hold5/6_table.txt, child_truth.py.

## 692k integration (recipe17d): why the closed state ends in two leaves (route9/integ_closed, read-only, ~0.3 core-h)
- **Merge-up:** equivalence rule, mdelta 0.000075.
  - At every node where the four leaves met (levels 1-7), their cross-half re-test stats sat above the shuffle threshold
    (0.115-0.139). So none of the pairs ever reached the equivalence check: shown distinct, not undecided.
  - At node r131 (level 1) the EM-group stats were 0.24 / 0.36 / 0.62 vs 0.126. r1311_s1/s2 is also a cannot-link pair.
- **Halves** (pair_maps.py): SNR_d at 40-20 / 20-12 Å; L/L_ref:
  - corner/r1312: 14.7 / 3.1; L/L_ref 0.87;
  - corner/sibling: 10.5 / 0.8; 1.17;
  - corner/r2322: 30.4 / 5.3; 0.96;
  - r1312/r2322: 24.6 / 4.8; 0.93.
  - All are localised.
- **10 Å views** (view_cmp.py, scorecard arm masks; leafmap = exact size-weighted mix of its stage-B states, so no pipeline
  confound; all top voxels lie inside both node masks):
  - Within-group differences are closed-arm occupancy only: corner vs sibling, and r1312 s02 vs s01 (2 Å from the arm).
  - corner vs r1312 (and vs r1312's 0.94-closed state s02, closed-arm mean diff -0.06): 58-70 % of the top difference sits
    off the arm, a median 22 Å from it.
  - r2322 differs from both about 30 Å from the arm, with more open arm.
  - The closed-open axis explains only 3-23 % of each pair's reproducible difference.
- **Read:** r1312 is a distinct closed-lineage sub-state, not the corner diluted by open particles. Its open-arm occupancy is
  unchanged (0.56 vs 0.58); its own internal variation is closed-arm occupancy.
- **USER VERDICT (2026-10-29):** "those 4 maps are distinct closed conformers" (r1311_s2, r1312_g1_2, r1311_s1, r2322_g1_2). The tree is right to keep them apart. The scorecard closed corner (>= 0.95 / <= 0.6) counts only r1311_s2.

## Locality-guard blind spot fixed (route9/guardfix; all gates passed; ~138 core-h)
- **Diff:** locality_guard_children_r17d.diff (md5 d6e45d9f) against recipe17d / 5709bc083. Fast gate 15/15.
  - When a delivered group pair has no spanning stage-B state pair, or the stage-B halves are missing, state_split
    reconstructs the children's half maps (Fourier-padded to the stage-B box 88) and judges the pair with the same
    pair_locality and 0.5 bar. It logs the path.
- **692k integration tree:** the blind spot never fired (59 split nodes, every group pair judged, all >= 1.01).
- **Gate 1, known instance qtoy2_inh r4:** the replay reproduces z 3.61 (5,074 / 1,647). The fix judges it: 0.639, kept.
- **Calibration** (15 nodes, stage-B vs children path):
  - nulls 0.20-0.24 vs 0.22-0.24; children-only nulls 0.26 / 0.36, all merged;
  - real splits 1.09-1.48 vs 0.87-1.01, kept; children-only reals about 1.03.
- **Gates 2-4:**
  - qtoy2 md2 identical to the recipe17d baseline (11 leaves, NMI 0.686);
  - qtoynull / olenull md0: 1 leaf each (guard merged at 0.205 / 0.224);
  - r131 md1: 4 leaves identical to the off arm; corner 0.97 / 0.54.
- Handed to the integration session (Assess and integrate) for recipe18c.

## Gold-standard node fix: SIMPLE_TREE_HALVES + SIMPLE_TREE_HALFSEED (2026-10-30, this session)
- Base fork hetero-tree 21fac34db (= route9/guardfix/src, md5 check in the build). Diff halves_halfseed_r17d.diff (md5 12bbfa0a, 2 comment lines),
  beegfs route9/halves/{src,build/install}; build job 54724894 on norm, fast gate 15/15. Both switches unset = 21fac34db behaviour.
- SIMPLE_TREE_HALVES=yes (matcher, read_mask_filter_refvols): with lp set, each half is matched against its own e/o half map (logs
  '>>> TREE HALVES: own e/o half maps for state s'); no fixposes restriction (works with prob_neigh re-alignment); no low-res e/o blend
  (that would share < 20 A between halves, the EM band). The tree's stA0/stage-A/stage-B inputs all have true _even/_odd beside volN.
- SIMPLE_TREE_HALFSEED=yes, chosen: seed each half separately and match labels (not cross-half assignment). Reason: a label made
  only from its own half's data stays independent through every later step (stage A/B with halves, site-layer chains); a cross-half
  label is re-mixed as soon as a step refits on own-half data, and any direction shared by both halves (flex on all data, a shared
  start tail) re-couples them through the selection effect (splitting along d makes both halves' group means differ along d).
  - per-half flex: fflex_e / fflex_o (flex_pca on one eo half, its rows re-split in two alternately) with the node's fflex; S_FFLEX waits for all three.
  - site_rank: starts from each half's own latent tails (merged file fflex_h_coordinates.txt), every endpoint EM fitted inside each
    eo half (state_mergeup set_eo: the half's own alternating split), odd endpoints greedily matched to even ones by site-map
    correlation, posteriors own-half; the halves accept rule then compares truly independent half versions.
  - two-group fallback: each half's groups from its own fflex states; odd labels swapped if the even/odd group-difference maps
    (40-10 A) anti-correlate.
- Gates launched 2026-10-30 ~02:40 (norm, sst_gate5, CONTRAST=maps, no floor): micnull_{off,h,s,hs} (0.47 deg / 0.31 px, re-align),
  qtoynull_md0_hs, olenull_md0_hs (jobs 54724934-39). Price ~290 core-h for the whole set.
- v2 (halves_halfseed_r17d_v2.diff md5 4f7724ed, build route9/halves/v2, fast gate 15/15): the half re-split for the per-half flex by
  a row hash. v1 alternated, and flex_pca's paired fits split by row mod 4: on index-parity eo (qtoynull, olenull) fit A got no odd
  particles -> "zero right-hand side" in its PCG. micnull (random eo) ran on v1 unaffected.
- NULLS (all one leaf):
  - micnull 0.47 deg / 0.31 px, re-align: off z 61.7 PASS K*2 (guard merged 0.221); h z -0.32 K*1; s z 52.3 K*2 (guard merged 0.214); hs z -0.91 K*1.
  - CONFOUND: the tree's stA0 re-alignment takes the fixture from 0.47 deg to 7.3-7.8 deg median (p90 ~133 deg, 25-34 % > 30 deg)
    in every arm; stage-B groups track micelle radius (AUC 0.32-0.36), never state. The good-pose regime (where the micelle session
    saw z 4.8 with halves) needs fixposes, not in this base -> per-half seeding is untested at good poses.
  - mic3 hs: K*1 (z 0.88), stage-B NMI 0.000 -> no positive control under re-alignment (states lost to pose scatter).
  - qtoynull_md0_hs2: K*1 (cross-half EM one group); olenull_md0_hs2: K*1 (every site's poles rejected by the halves rule).
- Positives launched 03:10: qtoy2_md2_hs, r131_md1_hs, ole_md2_hs (54727709-11); mic3_off/h for the control comparison.
- mic3 arms (re-align): off z 54.9 K*2 -> guard merged (stage-B groups = pose error AUC 0.60 + micelle radius 0.35, state NMI 0.001);
  h z 2.82 < 3.01 K*1 (NMI 0.000); hs z 0.88 K*1 (NMI 0.000). No state signal survives the re-alignment scatter in any arm.
- POSITIVES (v2, both switches on; score job 54733094, route9/halves/score; logs/smr_halves_posscore_54733094.out):
  - qtoy2 md2: 9 leaves (base 11), NMI 0.712 (base 0.686), quality NMI 0.000; best-leaf F1 c18 .89 (.87), c1 .79 (.54), c2 .71 (.65),
    c22 .69 (.74), c28 .64 (.74), c43 .56 (.51), c49 .73 (.80), c101 .94 (.94). PASS.
  - r131 md1: root splits as in base (22336 / 21998, z 9.43 vs base 19.86, guard 1.175); BOTH children fail the split test
    (r1 z 1.03, r2 z 1.18 vs z_crit 3.71; base r1 z 6.07, r2 z 5.25) -> 2 leaves, closed .86/.60 and .88/.52 vs the off/base corner
    leaf r2_s1 .97/.54. FAIL (closed corner lost). Mechanism: the beyond-EM band (shells 16..23) pair statistic falls 0.25 -> 0.075
    (r2: group-difference direction still replicates, coarse corr .94) while the null sd rises; r1's per-half groups also disagree
    (board corr .37, cross-half EM .17 vs base .99/.91). Halves vs halfseed not yet separated.
  - OLE md2: 6 leaves (r14 pair 21/25 at full depth); depth-2 cut 5 states (r14 10). r2 child fails (z 2.40 < 3.12), r12 cross-half EM
    collapses to one group (prints z Infinity, correct FAIL). Landscape: r14 covers hs .88-.91 but hs covers r14 .11; ends a2 .40..52
    vs r14 .35..56, a3 .46..57 vs .40..60; W1 vs r14 .065-.070 (r14 pair .044). Narrower landscape. WEAKER.
  - Cost: 314.8 core-h for the whole halves gate set incl. builds and scoring (r131 run 25, qtoy2 48, OLE 112).
- VERDICT: both switches kill every null false pass and help qtoy2, but fail the hard r131 closed-corner gate and shrink OLE: true
  halves remove the shared-reference correlation the split test's beyond-EM band relied on for depth-1 splits. Not a default.
  Next (needs OK, ~50 core-h): r131_md1_h and r131_md1_s to attribute the corner loss to halves vs per-half seeding.
