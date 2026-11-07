# hetero_tree on a continuum (session tree-continuum, 2026-10-29), final

## Verdict
- Baseline (recipe17d): a continuum gets ONE real cut (root, at the mass median); every child loses the motion axis (no pole -> fflex
  two-group fallback on an unrelated axis -> re-test fails). 2-3 broad leaves; maps near the poles (pose-pushed), none in t 0.15-0.7;
  far end not reached (0.76-0.78), dwell not represented.
- Remedy SIMPLE_TREE_INHERIT=yes (default off, build_inh; inherit_on_r17d.diff in this folder, md5 01764c0f): 4 leaves on both
  fixtures (the resolution limit at this SNR: every depth-2 re-test fails); closed end reached (maps 0.00-0.08 in all 4 runs), far end
  0.81-0.88 (just outside 0.1), ctpk s1 has a map at the dwell (0.68); W1 on maps -22 to -28 %. Reproducibility: ctun s1 = s2 (maps
  identical to 0.01, ARI 0.83); ctpk extremes reproduce (closed 0.076 / 0.081, far 0.875 / 0.810) but interior maps do not (map W1
  between seeds 0.112: seed 2's root took a route-A extra seed and split 3-way; same in baseline code).
- NOT fixed: occupancy. Leaf counts are near-equal and attached to maps the poses pushed outward; dwell underweighted 2.3x. Soft
  mixture weights over the leaf maps do not fix it on real data (pose-biased maps). Needs conformation-insensitive poses.
- Controls: qtoynull 1 leaf 100 % (unchanged by construction). qtoy2 changes at its 2 pole-less nodes (r1, r4): NMI 0.686 -> 0.685,
  closed c18 F1 0.87 = 0.87, c1 0.54 -> 0.77, c2 0.65 -> 0.56, holo c101 0.94 -> 0.85 (holo node r4 splits at z 3.61 vs 3.38 into
  two 92 %-holo leaves 5,074 / 1,312 + a 335-row 76 %-junk leaf; baseline fflex seeds z 2.65 -> one leaf).
- Guard blind spot (pre-existing, state_mergeup locality_guard): when choose routes to EM groups and every stage-B state's majority
  lands in the same EM group, no state pair spans the groups -> the guard skips silently and keeps K* (qtoy2_inh r4 printed only
  "K* 2 -> 2"). Not fixed here.
- The fflex fallback seeds 71 of 84 nodes of the 692k integration tree: the switch would change most of the RNase P tree; it needs the
  RNase P closed-state gates (r1211, closed_r111) before any default change. No default changed.
- Cost: ~340 core-h in total (sims 60, baselines 41, remedy 104 + seed 2 ~60, qtoy2 base/inh 50+41, diagnostics ~25, build 16).

Work dir (batch2): /mnt/beegfs/elmlund/afan/hetero_tree_20261028/continuum/ (scripts/, sim_<tag>/, runs/, diag/, logs/, RUN_ENV.txt).
Nothing committed. Baselines on recipe17d; the remedy on continuum/build_inh (recipe17d + inherit_on_r17d.diff, fast gate 15/15).

## Fixture (reusable)
- scripts/run_ctoy_sim.sbatch (env TAG, DENS=peak|uniform, PEAK/PSD/FLOOR, TRAVEL 30, NSTEPS 31, N 60000, JF 0.0002, EXTRA=conf,share,kind,map):
  qtoynull recipe (real RNase P poses + CTF via make_oris55.py, noise-free simulate_particles, qtoy_merge.py quality model B ~ U(0,200),
  q lognormal 0.35, SNR0 0.0012; refine3D prob lp 8 4 its from 10 deg / 3 px perturbed truth). Conformers = a rigid swing of the closed
  conformer's mobile domain (make_swing_maps.py, docked_all/conf18 + variants/domain_mask, 30 A tip travel in 31 x 1 A steps, ids 500-530).
- sim_ctpk: dwell peak t 0.65 sd 0.12 + floor 0.12: 564/step at the closed end (t 0), 5,244 at step 19, 632 at t 1 (59,978 rows).
  sim_ctun: uniform, 1,934-1,956 per step (60,000). Pose error 4.2 deg median (qtoynull 3.8).
- Truth: sim_<tag>/truth_qtoy.txt, column closed_arm = t (step / 30), s1 = tip travel A. Scoring only.
- Trees: sst_gate4, qtoy keys (nparts 4, mskdiam 232, minpop 900, site_layer/site_states yes, kway_maxdepth 2), CONTRAST=maps only.
- Scorer: scripts/ctoy_score.py leaves RUN [RUN2] truth=.. maps=SWING_MAP_DIR | nodes RUN. Leaf / unit member t (mean, sd, p10/p90),
  map content along t (NNLS fit of the 31 truth step maps, 12 A, top-1 % variance voxels; exact on member-averaged truth maps; the
  single best-step CC is biased for broad mixtures), end reach vs the same-size oracle, W1 of leaf-weighted positions vs truth (with
  best-K and evenly spaced K references), dwell share; nodes: seeds / split groups per t decile, re-test, choose, guard lines.

## Baseline (recipe17d) diagnosis
- Root splits 2-way at the mass median (ctpk z 52.2, guard 1.15; ctun z 24.9, guard 1.23). Route A adds nothing (sites' Spearman with
  the seeding site 0.998-1.000). Groups overlap: per-particle position error sd 0.33 of the span (root stage-B ev1-ev2: slope 17 nats
  per span, residual sd 5.8; best-quality third 0.27); evidence flattens at both ends; evidence tails of 900 reach only t 0.17 (ctun) /
  0.30 (ctpk).
- Depth 1: the children's sites give no pole; the fflex two-group fallback seeds groups with the same t distribution (30/70 in every
  decile); re-test fails (ctpk r1 1.54, r2 2.63; ctun r2 3.30) -> leaves. ctun r1 passes on a weak axis (t 0.305 vs 0.342).
- ctpk final: 2 leaves, members t 0.51 / 0.69; their MAPS sit at 0.22 (p10-p90 0.00-0.30) / 0.78 (0.70-1.00): the maps are pushed
  apart beyond their members (not quality weighting: < 0.01). No leaf at the dwell.
- ctun final: 3 leaves / 3 units: r2 members 0.69 map 0.76; r11 / r12 near-duplicates (members 0.31 / 0.34, maps 0.09 / 0.10;
  depth-2 re-tests 2.50 / 0.53 FAIL). No baseline map between t 0.15 and 0.7.
- Why the maps sit at the poles (diag/recpose_ctpk_s1, scripts/recon_poses.sbatch; NB import_particles resets state -> selection):
  the same two leaves reconstructed with true poses 0.475 / 0.720 (= members within 0.04), input refine poses 0.218 / 0.696,
  root stage-B 0.247 / 0.782, own stage-B 0.228 / 0.786 (= the tree's maps). Pose error grows with distance from the map's own end
  (leaf 1 own stage-B 3.4 deg at t < 0.3, 7.6 deg at t 0.9; mirrored in leaf 2), so far members blur instead of adding density.
  A broad leaf's map shows its well-aligned members; its count is not where its map is.
- The fflex fallback seeds 71 of 84 nodes of the 692k integration tree: any change there needs the RNase P gate.

## Remedy prototype (default off): SIMPLE_TREE_INHERIT=yes
- A node below the root whose sites give no pole seeds its binary split by the median of its rows' leading parent stage-B evidence
  contrast (PC1 of the parent's per-particle state evidence, prob_state_evidence.txt, mapped through node_rows.txt / rows.txt); falls
  back to the fflex two groups when the parent wrote no evidence. Re-test, choose and guard unchanged. Never at the root (nulls whose
  root merges are unchanged by construction).
- Parent contrast inside the baseline children: median halves t 0.236 / 0.409 (ctun r1), 0.605 / 0.778 (ctun r2), 0.645 / 0.726
  (ctpk r1), 0.444 / 0.579 (ctpk r2); the fflex seeds were equal.
- Seed 1 (md2), all depth-1 children now split along t (ctpk r1 z 5.72 / r2 10.67, ctun r1 12.76 / r2 17.20; guard 1.29-1.33);
  every depth-2 re-test fails (z -0.56 to 2.28): ~4 resolvable segments at this SNR. Cost 29-30 core-h per tree (baseline 17-23).
  - ctpk_inh: 4 leaves = 4 units; members 0.41 / 0.56 / 0.65 / 0.72; maps 0.08 / 0.30 / 0.68 / 0.87 (baseline 0.22 / 0.78).
    W1 on maps 0.130 (baseline 0.166), on members 0.076 (0.102; best-4 0.051, even 0.063). A map at the dwell (0.68, 0.63-0.67).
  - ctun_inh: 4 leaves = 4 units; members 0.25 / 0.40 / 0.61 / 0.77; maps 0.00 / 0.21 / 0.62 / 0.88 (baseline 0.09 / 0.10 / 0.76).
    W1 on maps 0.109 (0.152), on members 0.097 (0.133; best-4 0.065).
  - Ends (map within 0.1): closed end yes in both (0.08 / 0.00); far end just short (0.875 / 0.878).
  - qtoynull_inh: root re-test PASS 6.16, guard merges -> 1 leaf 100 % (unchanged by construction).
  - qtoy2_inh: root K* 4; r2 / r3 have poles (unchanged), r1 (22,811) and r4 (6,721) use the inherited seeding -> leaves change;
    recipe17d baseline qtoy2_md2_s1 running for the truth comparison.
- Occupancy is NOT fixed: near-equal leaf counts (median-like cuts) attach to maps the poses pushed outward. ctpk share per map vs
  the true mass in its Voronoi cell: 0.08 0.15/0.06, 0.30 0.31/0.14, 0.68 0.28/0.64, 0.88 0.26/0.16 (dwell 2.3x under).
  - 1-D simulation said soft mixture weights would recover it (0.65 vs 0.64) where hard routing would not (0.29-0.34).
  - Real data (diag/occ_<tag>_inh_s1, scripts/occ_estep.sbatch + occ_fit.py: one refine3D E-step of mergeup/final.simple vs
    the leaf maps, EM on pi): soft pi 0.12 / 0.32 / 0.33 / 0.22, W1 0.116 vs counts 0.130. Particles fit the map their poses
    were refined against, so likelihoods don't know where a map sits on t. The error is in the map positions (pose bias),
    not the counts: a final routing pass or mixture weights over these maps cannot fix it.
