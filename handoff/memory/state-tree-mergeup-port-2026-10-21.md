---
name: state-tree-mergeup-port-2026-10-21
description: "SIMPLE port of state_merge.py (prg=state_tree) and splitbench9+bu_merge (prg=state_mergeup), site_rank, landscape_compare harness and pole-seeded hetero_tree: uncommitted in worktree sharp-proskuriakova-b68280, batch2 ~/SIMPLE_stree snap/stree24; needs CONSENSUS poses and UNMASKED EM maps; swing30hs 16 vs 14 = SIMPLE per-particle noise normalization on unnormalized sim stacks"
metadata:
  node_type: memory
  type: project
  originSessionId: 58de2c8f-ad53-44fb-b8a2-05aadff16ec4
  modified: 2026-10-22T23:19:44.909Z
---

2026-10-21 (task: port the general heterogeneity-organization prototypes, see [[klem-bottom-up-organization-2026-10-20]]).
Code (uncommitted, Mac worktree .claude/worktrees/sharp-proskuriakova-b68280 = batch2 ~/SIMPLE_stree at 5703662e7, identical diff):
src/main/volume/simple_state_tree.f90, src/main/strategies/search/simple_state_mergeup.f90, commanders in
simple_commanders_refine3D.f90, UI simple_ui_heterogeneity.f90 (+tester count 7), param `linkage`, reconstructor
add_shell_rho/get_rho_shell_medians, gridding thread-safety fix. Snapshot for jobs: ~/SIMPLE_stree/snap/stree5.
Tests: rnasep_klem/stree_port (stree_accept3.sbatch = consensus-pose project + selection of the leaf states, then state_mergeup).
Results vs Python: state_tree identical merge order + d2/z on f692r 218 and the 3 panels; state_mergeup null_50k 3 = Python
(NMI 1.000), two_50k 4 = Python (NMI 1.000); swing30hs 16 vs 14 (NMI 0.954, swing 6.1/14.6/20.5 vs 6.1/14.9/21.7). The
Python gives 14 also on the port's eo partition, so the gap is the port: level-1 pairs 16+19 and 2+11 stay apart because
SIMPLE's KB-gridded maps resolve the leaves' view/defocus differences more sharply (stratified null 0.3-0.5 vs 0.05-0.13;
random-label control 0.03, so no coupling bug); EM moves 12-22 % vs 39-50 %.
Variant tried (user, 2026-10-21): EM insertion band-limited to the EM samples (2 <= r <= kb, as the prototype) = snap/stree6:
identical outcomes on all three controls (null 3, two 4, swing 16 with the same labels), 1.7-2.2x faster, kept in the tree.
Node 16+19 shows why: SIMPLE groups differ at ~58 A (deviation FSC shell 4/8 0.93/0.23, EM 1 % change) where the prototype's
are weaker (0.80 -> 0.47, EM 27-31 %/iteration).
DIAGNOSED (job stree_node6, rnasep_klem/stree_port/node6_diag): node 6 = leaves 16+19 = the WHOLE "blob" junk class split by defocus
(median 1.11 vs 1.51 um), so merging is right. simulate_particles draws each image's noise sd from its own CTF'd projection
(sd = sqrt(var/snr)) and its stacks are NOT normalized; SIMPLE's rec path divides every particle by its background sd, so the
low-defocus leaf gets 1.45x the signal at shells 3-6 (SIMPLE maps amp 16/19 1.38, 87 % of the difference = scaled consensus).
splitbench9n.py --simplenorm 1 (prototype + SIMPLE normalization) alone reproduces the port: null 0.13 -> 0.42, EM 50 % -> 18 %
moved, re-test 0.04 -> 0.47, K* 1 -> 2. Leaves 2+11 (other missed merge) = same conformers 48/49, defocus 1.56 vs 1.09 (not run).
Real SIMPLE-extracted stacks are normalized at extraction (simple_particle_extractor norm_noise), so the gap is sim-specific.
FIX applied (user, snap stree7): state_mergeup%new scales every plane by the background sd norm_noise returns (no-op on
normalized stacks): null 3 and two 4 unchanged (NMI 1.000), swing30hs 16 -> 15 (NMI 0.961): level 1 now = Python
(20 -> 17); left: level-2 node (units 2 5 6 8 9) null max 0.52 (Python 0.32) merges the 6.1 A swing end (port K* 4, Python 5).
**Why/How to apply:**
- state_mergeup must get CONSENSUS (pre-state) poses: on refine3D_states' state-refined poses it under-merges (null 6, two 8,
  swing 19 vs 3/4/14) - state-specific pose refinement makes leaves reproducible across halves. Put the leaf labels into the
  consensus project with `prg=selection infile=states.txt`; state_mergeup makes eo halves itself (ensure_field_halves).
- EM maps must stay UNMASKED like the prototype; masking them sharpened references and kept same-conformer units apart (null 6).
- TRAP (Python): state_merge.py's "15 A" envelope low-pass is exp(-2(15*smpd*k)^2), k cycles/px (~32 A at 1.464); ported as coded.
- TRAP (SIMPLE): `real, pointer :: x => null()` in a routine = implicit SAVE; deapodize3D_inplace crashed when called from threads.

PORT TASK 2 (2026-10-21, relayed by the 2D-to-3D session; spec /Users/atifao/SIMPLE/doc/implementation_notes/state_tree_port_spec_2026_10_21.md):
the worktree now also carries the uncommitted pcg_focus diff (copied from worktree vigilant-cartwright-dcd646, same base 5703662e7);
batch2 build tree ~/SIMPLE_stree (snaps stree8 = splitbench16 test, stree9/10 = + hetero_tree).
- state_mergeup's test = splitbench16 rules (SSNR junk screen + re-reference, Newton docking to the OTHER half's consensus on the
  tested band / EM band, per-pair equal N, re-test beyond the EM band with its own null, EM argmax over alive groups, em_min 1000,
  focus via pcg_mskfile = msk3 x read_and_crop'd mask, rb margin 4 with focus); prg=state_split tests one node.
  Validation vs splitbench16 (stree_port/stree_split1.sbatch): null50 / two50 / r11-foc identical partitions + junk; r22-foc K* 6 vs 5,
  f692st r1-foc K* 7 vs 8 (merges 4+6 at 0.15 < null 0.189; Python 0.18 kept); closed state 5 alone in both. 40-125 s per node.
- prg=hetero_tree (src/main/strategies/search/simple_hetero_tree.f90): master = node table + queue; every heavy step a qsys async
  job (exit file hetero_step.exit) in node_<id>/<step>/ with mkdir=no; in-process: pruned node projects (prune_particles, rows in
  master space), cheap children, falloff mask, decisions, leaf maps, leaf tree. Params kfloor/maxdepth/nodepar (+ minpop).
- TRAP (gfortran): a deferred-length character FUNCTION RESULT passed to a class(*) dummy (simple_mkdir(self%ndir(1))) segfaults;
  pass string(...) to simple_mkdir / file_exists / simple_chdir / simple_rename.
- E2E null_50k (stree_port/tree_null50, job 54051779, snap stree11, NP 4, mskdiam 220): root ran mean -> flex -> flex-state
  falloff mask -> fflex -> stage A/B -> state_split -> decision in ~40 min: junk leaf (unit 7, 4975), partial leaf (state 8, 4081),
  family r1 (40944) = splitbench16's "1 + partial + junk"; r1 cheap children 31 bins x 10 children ~4090 each, mask 95 % support.
  Steps must run simple_exec (qsys default exec is simple_private_exec -> reconstruct3D sigma2 error); set SIMPLE_EMAIL=nil.
  f692 r1 test: stree_tree2.sbatch TAG=f692r1 MAXD=0 (root project rflex_f692st/node_r/fstageB/2_refine3D/b.simple, rows_r1.txt).
- null_50k tree DONE (52 min, stree11): 3 leaves (1 junk) over 50000, NMI 0.71 purity 0.98; r1 K* 1 -> one leaf (no splits beyond
  junk = spec target). Merge-up resume (stree12): leaf tree + final project (leaf poses from each node's stage B) + state_mergeup.
- TRAP: sp_project%prune_particles THROWS "Too many stacks! get_mic2stk_inds" on the real f692 project (more mics with particles than
  stacks); node projects drop os_mic before pruning (refinement/flex/split need only stk + ptcl segments).
- f692 r1-equivalent node (hetero_tree, stree14, infile rows_r1, maxdepth 0; 57 min): K* 8, best closed leaf r_s8 11.1k at
  closed-arm 0.90 / open-arm 0.52 (script r1 state 5: 14.8k at 1.00 / 0.47). Node mask vs the script's for the same node:
  167k vs 167k voxels > 0.5, identical arm coverage (0.90 / 0.74), Dice 0.96, value corr 1.00 -> the mask port is exact; the gap
  is the stochastic flex / refinement. Node projects now keep ALL rows (prune corrupts f692: 149 particles' stkind vs stk ranges).
  Label check (stree_closedcmp, 54059768): the split made 0 merges (not the cause); port vs script r1 stage-B states NMI 0.24 on the
  same 92,554; script r1_s5 = 33 % port r_s8 + 33 % r_s3 + 17 % r_s4; port r_s8 = 44 % r1_s5 + 32 % r1_s6. Only a repeat port run
  on r1 (flex_pca / refine3D seed from the clock) tells run-to-run variation from a port difference.
- 2026-10-21 ~22:40: build 19 = snap stree15 (gated): state_mergeup infile2 = per-state focus masks (each group test inside the
  union of its leaves' node masks; hetero_tree's merge-up writes mergeup/leaf_masks.txt from node_*/mask.mrc). Runs: r1 seed rerun
  tree_f692r1_seed2 (54060016, stree14) -> stree_seedcmp.sbatch; f692st 64-leaf focused merge-up via hetero_tree RESUME
  (stree_mup_f692st.sbatch, 54060315: leaf table from rflex_f692st/leaves.txt + final_labels.txt, node dirs symlink the script's
  fstageB/2_refine3D/b.simple and rsmask.mrc).
- PUSHED 2026-10-21 (user: "push here: github.com/afanoa/SIMPLE_hierch_flex"): fork branch hetero-tree = 899189bf5 (pcg_focus,
  17 files) + b3e2bac72 (port, 11 files), base origin/master 5703662e7, 26 files +4459/-16, tree == gated build 19 (md5 of all
  26 files vs ~/SIMPLE_stree). Fork master abeeebfba is an ancestor (branch 17 ahead). Uncommitted after: pair_jobs warning fix.
- SEED TEST (stree_seedcmp 54062662; per-step via a copy of toy55/repro_r1.py = repro_r1_seed2.py, labels stree_port/seedcmp_labels):
  port run 1 vs rerun NMI whole flex 0.65, fflex 0.61, stage A 0.55, stage B 0.53 (closed group Jaccard 0.60) vs port vs script
  0.66-0.82 / 0.35 / 0.25 / 0.24 -> NOT the same range (the peer's audit predicted it would be); the gap opens at the focused flex.
  Same flex/focus code (SIMPLE_focus vs SIMPLE_stree differ only in port files), same fflex keys (port mkdir=no), script flex labels
  all 92,554 (no random fill). Port's closed leaf reproduces at 0.87-0.90 / 0.49-0.52 vs script r1_s5 1.00 / 0.47 (focused readouts).
  Masks don't explain it (stree_maskcmp3): Dice run1-rerun 0.976, rerun-script 0.971, run1-script 0.956 (value corr >= 0.995).
  Inputs identical (stree_inputcmp/2): every ptcl3D key over the r1 rows equal script vs port and run 1 vs rerun (poses, eo, dfx,
  history fields); rflex13's focused-flex selection (flex_rows.txt) keeps exactly the node rows. flex writes state/cfar/res/res05.
  -> nothing systematic found; third port replicate offered to the user to tell chance in the one script run.
- 2026-10-22 ~00:20 (user via peer): third run tree_f692r1_rep3 (54063860, stree14) -> stree_rep3cmp.sbatch (four-run per step);
  greedy default ATHR prob_athres=0.01 in refine_keys = build 20 = snap stree16 (gated, uncommitted, + pair_jobs warning fix).
- SPLIT-THEN-MERGE (stree_mup_r1.sbatch, dirs stree_port/mup_r1_<run>, report mup_r1_report.txt): the focused merge-up of each r1
  run's 8 leaves merges NOTHING (script, port1, port2: 8 -> 8, every group beats the shuffle null in both bands) -> the map-level
  cross-half test cannot reject arbitrary cuts of a near-continuum; not a width rule. Needs a particle-level criterion (best-vs-second
  cross-half gap >= 3 nats share) or restart co-association - user's call.
- USER 2026-10-22 01:05 (via peer): "the number of supported splits needs to be data derived. the conformational landscape needs to
  be reproducible as a whole." Peer recommendation: doc/implementation_notes/state_tree_reproducibility_overnight_2026_10_22.md (main
  repo). Plan sent (~8 working days + overnight acceptance): 1) K by restart stability (2nd independent fflex, PC1-bin seeds per K,
  light stage A box 64 subsample greedy, pass = NMI after >= seeds and most distinct group Jaccard >= 0.75; K=2 fails -> leaf), 2 d;
  2) tree_compare harness (per-level NMI/ARI, leaf/sample maps matched by FSC vs own half-map FSC, landscape coords, site rankings)
  + 3 full runs on two_50k/s30a/f692, 1 d + runs; 3) site ranking port (2-class soft EM from both latent tails, parent poses) 2 d,
  landscape sampling port (axis = the stable group's latent; needs user OK) 3 d; greedy done (stree16).
- f692st FOCUSED MERGE-UP (stree_mup_f692st 54060315, 59 min, merge-up job 18.9 GB MaxRSS at 643k): 64 leaves -> 28 units, 0 junk,
  merges only WITHIN families (r8 8->1, r7 7 (+r7_s4 alone), r4 5, r6 4+2, r2 4+2, r1 4 (s3/s6/s7/s8), r3 2, r5+r51+r52 8 incl.
  r51_s5 0.98/0.79 = the other arm position SWALLOWED into a 72k unit); r1_s5 SURVIVES alone (1.07/0.51); 18 singletons.
  'Neither position' leaves: r6_s1/s6/s7 one unit, r6_s5/s8 one unit, r6_s2, r6_s4, r4_s6 alone, r4_s1/s7 in r4's unit, r2_s4 in
  r2's unit; none junk-flagged. Report: stree_port/stree_mup.54060315.out.
- SPEED MEASURED (rep3 r1 node, sacct of the parts listed in node_r/<step>/stderrout/outfile.<jobid>): part queue wait ~1 s per
  batch, 0.5 min per node in 21 batches -> single-allocation steps NOT worth it (peer's 6-8 min estimate wrong). The node's
  overhead is IN-PROCESS: flex 6 min + fflex 5.5 min (flex_pca's in-process delivered-state recs), cheap 2.4 + recB 2 min
  (per-state assembly). Trimming node projects would save ~2-3 min per r1 node (needs prune_particles grouped by stkind + a
  hetero_tree row-map refactor). Job names: tree parts are 'simple_afm_cryoem' like klem workers -> use the outfile ids.
  flex fix (needs the user's OK, flex_pca is theirs): the flex master's threads come from SLURM_CPUS_PER_TASK
  (set_master_num_threads), but simple_flex_pca_state_service builds its in-process rec service with params%nthr (8) -> use the
  master's count there + give hetero_tree's flex/fflex step jobs more cpus.
- FINAL MAPS PORTED (uncommitted, build 22 = snap stree17): hetero_tree final phase = plain() (one plain reconstruct3D per node with
  leaves, states = leaf index in LEAF_TABLE order, node_*/plain/) -> final_maps() (state_tree over the plain leaf maps -> new module
  src/main/volume/simple_tree_wiener.f90 = tw_units.py --stage maps at the map box, no padding) -> mergeup(); POLL_SECS 5.
  VALIDATED on f692st (stree_tw_f692st.sbatch MODE=tw, script plain maps linked): raw/own/emp2 corr 1.0000 (rms 0.0003-0.0043),
  tree2 min 0.9992; tree R^2 0.659 identical -> emp2 delivered as the script. 64 units 16 min (768 map writes), 7.7 GB RSS.
- SITE RANKING PORT (build 23 = snap stree18, gated): prg=site_rank (commander in simple_commanders_refine3D, UI heterogeneity,
  count 10) -> state_mergeup%site_em (soft 2-class cross-half EM, deviation x site mask in real space, 10 its, em_lim 15 via
  new(em_lim=), estep pfloor 1e-3 + beta clip [0,2]; N* = argmax hold-out (kb+2..khi) cross-half corr of the masked top-N
  deviation over SITE_GRID) from both tails min(15000, n/6) of infile (one coordinate per row). Validation stree_site_r1.sbatch.
- MASK TEST DECISIVE (stree_fflex_mask.sbatch 54067114): port fflex code + the SCRIPT's rsmask -> vs script fflex NMI 0.95,
  closed group Jaccard 0.98; two seeds 0.93 -> focused flex is ~deterministic GIVEN the mask and hypersensitive to it (port masks
  Dice 0.96-0.98 with the script's -> 0.31-0.61). The script-vs-port gap = each run's own mask (from its own cheap maps).
- site_rank vs prototype r1 z1 (54067416): lo pole Spearman 0.83, hi 0.67; same poles (old-closed enrichment 5.0/4.4/3.2x vs
  5.1/4.6/3.2x; top-2500 maps 1.12/0.42 vs 1.16/0.44, hi 0.46/1.05 vs 0.44/1.16); N* lo 25000 (proto 15000), hi 60000 (60000).
  Self-consistency test (z2, z3 tails) = stree_site_r1z.sbatch 54068198.
- PLAIN final maps on the 3 port r1 runs (stree_final_r1.sbatch): most closed plain leaves 0.92-0.94 / 0.46-0.56 (script r1_s5
  1.066 / 0.507); tree R^2 0.11-0.21 -> emp2.
- SYNTHETIC TRUTH TREES (stree14, old 8-cap/15k/unfocused merge-up; 54058308/09): s30a 20 leaves -> 14, NMI 0.41/0.43, purity
  0.34, swing eta2 0.31 (10-19 A, 2 states), closed 0.93 (61k) + 0.59 (79k); scripts rflex_s30a 34 states NMI 0.42 purity 0.41
  eta2 0.38 (7 swing states, closed 0.96/0.94/0.92 at 15-21k), s30c 0.42/0.38/0.33, s30split5 0.45/0.36/0.22. sim55_200k: port 14
  -> 10, NMI 0.21 purity 0.23, closed best 0.19; klem consolidated 0.07-0.08. s30a took 4 h 19 min.
- site_rank nearest-neighbour rewrite (build 25) gave NaN: empty shells' eps 0 -> 0/0; prototype falls back to eps 1 -> build 26.
  BUILD 26 (snap stree21, gated) MATCHES THE PROTOTYPE: NN Fourier site EM (crop-grid insertion/lookup, SIMPLE rotation row
  convention loc = h R(1,:) + k R(2,:), half-space with h=0 Hermitian, median-density eps x weight sum, beta clip [0,2], pfloor
  1e-3); scale 1.81/1.54, class 0.17 -> 0.26, argmax changes 0.136 -> 0.012, N* 15k; z1lo vs z3lo 0.981 (proto 0.989); each
  endpoint vs prototype 0.965-0.978; closed top-2500 1.14/0.45 (proto 1.16/0.44). Lesson: KB-gridded reconstructor maps (builds
  23/24) changed the EM dynamics (classes grew to 50-79 %); regulariser scaling cannot fix that.
  CROSS-RUN EXTREMA REPRODUCE (stree_site_xrun.sbatch 54069096): site_rank at shared input poses + shared site_02.mrc from each
  run's own fflex latent (script, port1-3): closed pole groups at Spearman 0.999-1.000 (top-1k Jaccard 0.94-0.97), open
  0.998-1.000 (0.95-0.99); closed top-2500 maps 1.052-1.055 / 0.465-0.473 in all four runs; N* 40k closed, 25k open (latent
  signs flip between runs; two tails handle it). Open item: per-run data-derived site masks (factorial v2).
- HARNESS prg=landscape_compare (2026-10-22, build 28, snap stree23; src/main/volume/simple_landscape_compare.f90): runs file
  "name mapdir labels" + pcg_mskfile node mask + lp 15 [+ infile2 "name site_post" poles]. REPRODUCES landscape_cmp.py: r1 5 runs
  (stree_lcmp2 54069733) |diff| 0.000 hdr/(1)/(2), 0.001 (3); seeded p2 exact; p3 exact except (2) where the Python prints a spurious
  null cosine (deviations sum to 0 -> rank K-1; harness keeps non-null components only). Extras: (1b) matched-pair cross-run even x
  odd DEVIATION FSC vs own; (4) extreme maps at both ends of pooled PCs above the noise share sum p(1-rel)|d|^2/sum p|d|^2, matched
  vs noise + membership Jaccard + FSC: r1 8-leaf script vs ports FAIL at PC1-low 0.02-0.25 (Jaccard 0.04-0.09), ports 0.82-0.97;
  seeded p2 1.01 / Jaccard 0.98-0.99; p3 closed 1.01, opens 0.96-1.01. (5) pole rankings grouped at 0.99 = site_xrun Python exactly.
  Own-pixel maps: (4) tests agreement, not extremeness ([[selected-subset-maps-imprint-inflates-extremes-2026-10-22]]).
- POLE INTEGRATION (build 29, snap stree24; peer design 2026-10-22): hetero_tree node = mean, flex, mask, fflex, SITE (site_rank
  infile=fflex/flex_pca_coordinates.txt site_dirs=1,3 [new key], pcg_mskfile = hetero_tree's pcg_mskfile else node mask), SEEDS
  (reconstruct3D nstates=P from site_seeds.txt), stageA (refine3D nstates=P, vol_j = seed maps), recB, stageB, split. merge_poles
  (state_mergeup): endpoints single-linkage at Spearman >= 0.95 (avg ranks, simple_stat avg_ranks), pole posterior = mean of
  members, seed = argmax (peer's poleseed*.txt used ONE endpoint per pole; mean accepted). Peer: every reproducible-pole result used
  the DISCOVERY factor mask (r1 site_02 = 1f1, ~21 %); broad focus mask as site mask gave a less pure closed pole (0.80/0.69).
  Test stree_poles_r1.sbatch 54070208: 3 runs resumed at fflex from port1-3 latents (node table "r 0 4 root -"), maxdepth 0.
  Also fixed: cheap_children read nlines() lines incl. the '#' header -> dropped the last particle.
- INTEGRATED r1 (54070208, stree24): site -> seeds -> stage A/B ran in 3 runs (fflex latents of port1-3); P = 3/3/4 (p3's closed
  endpoints z1lo/z3lo at Spearman 0.929 < 0.95). p1 vs p2 stage-B maps within noise (fit 0.93-1.00, match 0.97-1.01, FSC = own),
  Jaccard 0.81-0.88, NMI 0.69; p3's extra state (2.1k seeds -> 17.8k) unmatched (0.45). Cross-half pole identity (54070882): z1 closed
  across runs 1.00, merged within-run pairs 0.97/0.98, p3 pair 0.95 -> the merge line is a morning decision (peer). state_split crashed:
  see [[gfortran-optional-scalar-speculative-deref-trap]] (fixed build 32, snap stree27; resumed 54071622).
- DISCOVERY PORT (build 31, snap stree26): discover_sites in simple_hetero_tree (shares half_variance with falloff_mask) + prg=site_discover
  dir= infile= mskdiam=. Reproduces discover_sites2.py on disc_port1/2/3 (54071387): same factor splits (cos 0.17/0.21/0.13), every
  candidate's core voxels (+-2), shares, molecule %, centre, repro, seed identical, mask Dice 1.000; 5 s / 0.5 GB. Peer rules built in:
  rank by variance share (1f1 first in all runs), local masks only (<= 30 % molecule), Dice >= 0.9 = one site, siblings never merged.
  Open (peer, by morning): which state set feeds discovery before stage A (fflex vs cheap vs recB, factorial job 54071099).
- WHOLE-OUTPUT ACCEPTANCE (54072183, stree28): units = leaves = stage-B states (maxdepth 0, merge-up merged nothing); p1 vs p2 (P = 3)
  cross-half 0.96-1.00 at stage B / plain leaves / emp2 final maps, Jaccard 0.81-0.88; p3's extra state 0.42-0.47 everywhere. emp2
  half-maps are Wiener posteriors: reliability ~1.00 and FSC 3.7-5.6 A are NOT a noise yardstick. Section (4) now matches each
  extreme to its best cross-half counterpart with its rank (open A/B near-ties showed as -0.3 before). Peer: P = 4 seeding also
  reproduces (z1-closed -> intermediate 0.63/0.45 "both arms weak"); recommends P = 3 at the node, recursion to resolve the rest.
- pole_merge KEY (build 34, snap stree29, gated): site_rank/hetero_tree pole_merge=<Spearman> (default 0.95, unchanged pending the
  user); site_spearman.txt; hetero_tree logs poles + endpoint Spearmans. 0.99 on p1 -> P = 4 (z1-closed 3.6k, z3-closed 11.3k).
  Depth-1 test (peer request): copies poles_r1d1_p1/p2 resumed at r1's split with maxdepth 1 (node table "r 0 10 root -"),
  site_02 for every node, job 54072508; per-step r1 node cost: site 1.6, seeds 1.5, stage A 5.1, recB 1.4, stage B 7.5, split 2.1 min.
- site_disc OPTION (build 36, snap stree31, gated; default no = unchanged): node chain fflex -> stA0 (stage A from fflex labels)
  -> recB0 -> discover_sites on recB0 (in process, top site by variance share = sites/site_01.mrc) -> site -> seeds -> stage A
  (K = P) -> recB -> stage B -> split. Factorial (54071165): recB re-derives 1f1 in every run (Dice 0.92-1.00, pole Spearman
  0.997-1.000); fflex states 2/4 runs wrong component, cheap children never. Step numbers changed again (stA0 5, recB0 6):
  never resume node tables from builds < 35 with stree30+.
- DEPTH-1 (54072508, stree29, site_02 for every node): children correspond (Jaccard 0.81-0.88) but their splits do not
  reproduce (closed child 0.36-0.78, open 0.18-0.74); child poles degenerate (4-518 seeds beside 10-33k). p1 r3 stage A died:
  refine3D needs > PROB_MULTISTATE_MIN_POP = 5 per state. TECHNICAL guard (build 38, snap stree33): merge_poles drops endpoints whose
  pole draws < 6 seeds (refine3D's exported constant; the user rejected size floors, so never EM_MIN/minpop here); a node with < 2
  poles seeds stage A from its fflex states (no lost particles). CHILD DISCOVERY (stA0 -> recB0 -> site_discover -> site_rank,
  54074522/23): closed child's top site reproduces (1f1, 37 % var, Dice 0.90), z1 pole rankings 0.88-0.95 / maps 0.97-0.98, z3
  0.84-0.88 / 0.88, but argmax seeds do not (3.9k vs 10.2k); open children: no reproducible site (Dice 0, Spearman ~0).
  TRAP in quick scripts: parse print_project_field states by token (awk $i ~ /^state=/), never substring 'state=' (other keys
  match); flex_pca does write hard labels into fflex/f.simple (project_states is right).
- 2026-10-22 USER DECISIONS (09:55, push confirmed by the user directly): pushed 8421e1a31 to fork/hetero-tree (gated stree33);
  merge rule = pole identity on independent halves; site_disc=yes default; no edits to src/main/flex (falloff-prior session) or
  the port spec. HALVES RULE (build 40, snap stree35, default pole_merge=0): halves = row pairs ((row-1)/2 mod 2, as
  half_poles.py), site EM per half from the same tails, half version = top-N* deviation map over the hold-out band; merge if
  min cross >= min self (noise-normalised by Spearman-Brown eo reliability); accept a group if every member's best RAW partner
  (both directions) is inside it; seeds = full-node posteriors over accepted poles. Root r1 (54079707): P = 3 in p1/p2/p3 (p3
  by rejecting z3lo: 0.862 vs 0.868 - thin margins); closed child: ALL endpoints rejected at half data (raw 0.0-0.15). <2
  accepted poles -> node refined as ONE state -> leaf (split skipped). site_rank infile2 = extra starts (name + per-row 0/1
  file), site_starts.txt = basins. Depth-1 test with new defaults: job 54079790 (poles_r1d1h_p1/p2; 54079751 failed: copy
  source lacked fflex). QUARTER SE (build 41, snap stree36; site_halves.txt merge/accept lines with margin, SE, ratio; root
  54079805, halves_se/): open A vs open B only -0.74 / -2.12 / -1.66 SE (p1/p2/p3); closed twins +0.71 / +0.12 / -0.05 SE
  (null-like); closed vs open -3 to -21 SE. A 2-SE keep-apart margin would give root P 2/3/2: worse than margin 0 (3/3/3).
  TOY CALIBRATION (54079998 + score 54080065, halves_toy2/; factorial toys_for_port truth): halves merge/accept OK (simsites C
  closed AUC 0.945, twins merge at ~0 SE), but ARGMAX SEEDING FAILS: site EM = minority vs rest, so a binary switch gives ONE
  pole (simsites A 70/30: P=1 -> leaf, switch missed) and 91 % of particles have best pole posterior < 0.5 (C: closed seed 137k
  at precision 0.12). Twins missing the merge can BOTH be rejected (IgM arm 4 lost an extreme). Proposed to the user (not
  implemented): rest state (join a pole only at posterior >= 0.5; RNase P roots get an 8-10 % rest, K 3->4) + twin rescue.
  USER APPROVED BOTH (relayed 2026-10-22): build 42 = seed_poles helper (rest = last label, guard 6) + twin rescue ("best
  partner" = best OTHER endpoint; with self counted no twin case is rescued). CORRECTION (factorial census, 20 starts): simsites
  A has BOTH switch sides as attractors; z1/z3-only P=1 was start COVERAGE. SITE_ITS = 10 (= census iters) is NOT converged:
  argmax change at it 10 RNase P root median 5.6 % (max 10.7 %), toys 0.4-1.4 %; poles are 10-iteration states. Proposed (not
  approved) SITE_ITS 30 / to-convergence test (user's call; do not start until answered).
  BUILD 42 (snap stree37, gated) TOY RERUN (54080711, halves_toy3): C closed seed precision 0.12 -> 0.63; A K 1 -> 2 (deleted
  seed 0.95 pure, recall 0.30); IgM arm-4 lost extreme recovered by twin rescue; swing S1 continuum: rest = the MIDDLE (15.4
  deg between ends 9.7 / 20.7). A 6-seed pole survives the guard (C disc). Roots p1-p3 + p1/p2 depth 1: job 54080712 (rest_r1_*).
  TRAP FIXED (build 43): a queued step writes its .job record only when slurm STARTS it; hetero_tree polled every 5 s, saw
  "not started" and re-ran advance -> duplicate submissions clobbered seeds.txt (old depth-1 54079790 crashed at r2 in both
  runs). Fix: in-session list of submitted dirs (pending) in step_status / plain / mergeup. Old depth-1 children: r2 lost a
  strong pole to TWIN DOUBLE-REJECTION in both runs (3 twins all rejected) -> leaf; twin rescue targets it.
  ROOTS p1-p3 NEW RULES (54082284, build 43 snap stree38, rest_r1_cmp): K=4 in all runs (closed / open A / open B / REST),
  every stage-B state within noise across runs (cross-half 0.94-1.00). Rest: seeds 8-10 %, stage B 21 %, Jaccard 0.71-0.79,
  readout 0.57/0.56 (both arms half) in all runs; NOT the census middle corner (top-N Jaccard 0.03-0.05). Open A/B -1.0..-3.7 SE.
  30-ITERATION census (factorial session, simsites A): 19/20 starts on the two true states (Spearman 1.000) vs 15/20 at 10 its;
  same states, less scatter. Supports SITE_ITS 30 / to convergence; still awaiting the user's go (C closed-class growth untested).
  DEPTH-1 NEW RULES (54082284, rest_r1_d1cmp): no crash (race fix ok), p1 7 / p2 6 leaves; root K=4 reproduces but child splits
  do NOT (closed: leaf vs split; open A sub-maps cross-half 0.55; open B 0.19-0.45; rest child leaf in both). Degenerate 14-43
  seed poles pass the 6-seed guard at children. Final leaf maps match 0.76-1.00 of ceiling (they stay near parent states).
  CONVERGENCE (user approved, relayed): build 44 = stop at argmax change < 0.5 % (cap 50) -> FALSE STOPS at it 2 on simsites
  A/C (assignments stall it 2-4, then take off); build 45 (snap stree40) adds a 10-iteration floor. Converged r1 roots (conv_r1,
  med 20 its): P=2 closed vs OPEN in all 3 runs, open A/B merge (within 1 SE), rest vanishes (9-79 seeds) -> the 4-state root
  relied on 10-it early stopping; swing S2 middle pole and IgM arm-4 third pole also slide onto the dominant split.
  Build-45 toys (halves_toy5): C true unchanged (no over-growth); C disc pole becomes a 194k SINK; simsites A P=0 by a
  NON-CANDIDATE VETO (z3hi halves split noise/real: its half-2 copy out-matches the real pole's own halves -> all rejected).
  3-SE read-off: conservative (S3 rare + C disc lost), roots P=2 converged. Proposals 1-4 sent; awaiting decisions before
  making 3-SE default / running depth 1.
  BUILD 46 (snap stree41; decisions relayed 2026-10-22: 2+3 user-approved, 1+4 coordinator's calls): merge_poles_halves merges
  unless HALVES_KSIG=3 SE apart and accepts at >= 3 SE (quarter-pair SEs inside; module helpers halves_merge_margin /
  halves_accept_margin / halves_se); NON-CANDIDATE FIX (partners = candidates only); REST TEST: rest ranked by 1 - best pole
  posterior, site_em(score=...) gives its top-N* half/quarter versions (NOT raw tiny-set maps: nn_map of an empty/tiny set
  -> deviation ~ -mean, spurious cross-half corr), tested as endpoint 'rest' with l_rest=.false. (rejected -> best pole).
  RUNS ON BUILD 46 (2026-10-22 ~14:50): (a) r1 two-run depth-1 54084886 (s3d1_p1/p2); (b) toy trees swing S2 54084887 (maxd 2,
  minpop 2000), IgM 54084888 (maxd 2); (c) simsites C-disc 54084889 (site_disc=no, site_02 mask, maxd 1) -> toytree_<tag>/;
  FULL 692k TREE (user priority: one tree, all 692,392 rows of f692/subset_all.txt, no pre-filter): job 54084961, full692/,
  maxdepth 20, minpop 3000 (technical minimum of the quarter SE; USER CONFIRMED), step jobs on norm (nparts 16, nodepar 3).
  PROFILE (prof_tree.py, step-folder mtimes): 92k node from stA0 ~47 min = stA0 11.9 + stageA 10.0 + stageB 11.5 (70 %);
  small node (19-33k) 30-44 min: flex 4.9 + cheap 2.5 + fflex 6.0 + refines 13; master write_states / project copies ~4-5 min
  per node (932k-row master); site_rank 2-5 min at 92k (est. ~37 min single-process at 692k). DEV SUBSET: 100k random rows
  (seed 20261022) f692/subset_dev100k.txt + stree_port/dev100k_rows01.txt. Optimisation options 1-6 sent; awaiting choice.
  BUILD 47 (stree42, results-identical speedups): set_parts = min(nparts, max(4, ceil(n/20k))) per node, site step nthr
  max(nthr, 32), write_states = file copy + ptcl3D/ptcl2D write_segment_inside (node_project too). Verify job 54085368 (v47a/b
  vs s3d1_p1), dev100k baseline 54085369 (nodepar 10, nparts cap 32, norm). Build 48 switches (flex_discover, parent_focus)
  were REVERTED: the user moved the result-changing speedups (fixed-pose interior nodes, subsampled flex/site EMs, one flex per
  node) to a NEW chat; mine = build 47 + verify + full692 + dev100k baseline + landscape_sample merge + junk rule; ask the
  coordinator before touching advance()/seed_poles; write .claude/patches/hetero_tree_b47.patch after verify (EXCLUDE the
  landscape_sample patch: temp GIT_INDEX_FILE + git apply -R --cached lsamp patch, then git diff --cached HEAD). BUILD 49
  (stree44) = b47 + landscape_sample (sampling session's patch ~/lsamp_jobs/landscape_sample_8421e1a31.patch, new source ->
  cmake . rerun). BUILD 50 (stree45) = b49 + landscape() final phase (OK'd; after plain(), before mergeup; per finished node
  P>=2 -> landscape_sample in node/landscape, P<2 -> empty pairs + exit 0; summary log once). Gate job 54085866 (lsgate_p1/p2
  skeletons of poles_r1 trees, bitwise vs lsamp/repro_r1/fortran). B46 r1 ROOTS (s3d1): K=2 closed/open in both runs, rest
  merged into closed, stage B 38k/54k, cross-run 1.00 / Jaccard 0.96-0.97. Resume with new nodepar:
  safe only with NO pending step jobs (new master's pending list is empty). full692: switch to b47 + nodepar 10 at the root
  S_DONE if verify passes.
- 2026-10-22 ~15:00: verify47 PASSED (b47 vs b46 seeds 0.999 / stage-B 0.994 = run-to-run level; not bitwise); patch written to
  /Users/atifao/SIMPLE/.claude/patches/hetero_tree_b47.patch (base 8421e1a31, no landscape, no junk). BUILD 51 (stree46) = b50 +
  JUNK SCREEN (spec doc/implementation_notes/per_node_junk_detection_2026_10_22.md in worktree laughing-carson; prototype
  routes/toy55/junk_score.py xd3): state_mergeup%junk_screen(seeds) [12 A band jhk/jsh/jw + f2j per shell taken in new() before
  the CTF; other-half refs cons_o + FSC_g * (state_og - cons_o), no mask; beta/w per half at the labels; NUM/DEN per state, g* =
  argmax NUM - DEN/2, a = NUM/DEN, v = 1/DEN; xd3 EM 500 its; junk -> seeds 0], one call in exec_site_rank before site_seeds.txt,
  key junk_screen (default yes; parameters next to json, parse next to json), hetero_tree advance(): forwards the key, logs
  site_junk_report.txt, K<2 fallback resets seeds only if ks < 1 (K=1 keeps junk zeros; K=0 not screened). Outputs
  site_junk.txt (row flag amp var best post) + site_junk_report.txt. Gate 54087132 OK. Validation job 54087151 (junkval/:
  simsites C disc/true vs factorial/junk_test/simsites_truth_jt.txt, r1 p1-p3 expect 0). Junk patch on b47:
  scratchpad junk_screen_on_b47.patch. BUG FOUND (S2 toy tree, b46): the REST endpoint (quarter SE 0.3-1.7) "merges" with every
  pole and single linkage chains all poles into one group -> K=1 (truth 3 clusters); proposed fix: rest never merges/twins,
  accept-or-reject only (gives K=3 on S2 by the margins). Also: identical endpoints get SE ~0 -> "apart" (twin rescue covers).
- 2026-10-22 ~16:00 HALVES RULE FINAL (user pre-approved via coordinator; replay-validated on 12 nodes): "c+last+twincomplete" =
  COMPLETE linkage (groups merge closest pair first only if EVERY cross pair within 3 SE), merge SE floored at the self-matches' SE
  (halves_floor), strict candidates (raw self >= 3 SE above 0; halves_cand), REST SOLO + LAST (never merged/twinned, no partner in
  the poles' acceptance, tested afterwards vs the accepted poles only), twin rescue only if the twin is within 3 SE of every member
  of the partner's group. Why: single linkage let a noisy endpoint bridge (S2: rest; dev100k: z1lo) -> K=1; rest-as-partner vetoed
  the r1 closed pole; twins chained the s3d1 closed child. merge_poles_halves(..., l_rest, solo, grp, ...) solo REQUIRED (optional
  scalar trap). site_halves.txt: merge lines + floor column, quarter replicates appended; replay = scratchpad rule_replay.py.
  Results: S2 K=3 (truth), C-disc P=2, IgM P=3, r1 x3 closed/open, dev100k z1hi 29k | rest 71k, 692k root P=4 (z1lo 71k |
  z1hi+z3lo 407k | z3hi 181k | rest 15k), s3d1 children leaves. BUILD 53 (stree48) = junk + rule (worktree); 47b tree
  ~/SIMPLE_stree47b (build_stree47b.sh, gate stree_gate47b.sbatch) snapshot s47c = b47 + rule, no junk/landscape. Patches:
  .claude/patches/hetero_tree_b47b.patch (base 8421e1a31) + junk_screen_on_b47b.patch (on top). Junk at the 692k root: 17,346
  (2.5 %, mode 0.37). Full-depth toy trees on stree48 (maxd 20, minpop 3000, junk+landscape): fulltoy_s2 54087725, igm 54087726,
  ctrue (site_disc=no mask_C) 54087727, cdisc (site_disc=yes) 54087728; scorer toy_score.py. full692: restart at root S_DONE on
  s47c (NODEPAR 10, NPARTS 32).
  UPDATE ~16:15: CPU FIX (fast-tree found it): since b47 the site step ran 32 threads on nthr CPUs -> hetero_tree has a 2nd qsys_env
  qsite (qsys_nthr = max(nthr, SITE_NTHR)); submit_in uses it when cl nthr > master nthr. BUILD 54 = stree49 (worktree), 47b-3 =
  s47d; patch hetero_tree_b47b2.patch (+ junk_screen_on_b47b.patch on top). Toy trees on b53 were STOPPED (dirs kept as
  fulltoy_*_b53) and RELAUNCHED on stree49: s2 54087891, igm 54087892, ctrue 54087893, cdisc 54087894. Rule check on b53
  (rulecheck2/): Fortran == replay on all 11 nodes. full692 moved to 47b EARLY (safe: a resumed master waits for RUNNING step jobs
  via hetero_step.exit.job; only queued-not-started ones get resubmitted): master now 54087898 on s47d via stree_full692r.sbatch
  (TRAP: stree_full692.sbatch resume used projfile=root.simple but the tree's root project is rnasep.simple -> instant failure).
  The root's site step 54086775 (build 46, 8 threads) finishes ~17:15; same 4 poles as the new rule (replay).
- 2026-10-22 ~18:20 USER (direct): "you can merge it to what we have. so we only have one branch to move forward with. But no
  pushing" -> committed the gated b54 state as 15cf8c90e on claude/sharp-proskuriakova-b68280 (no trailer) and merged
  origin/master (had moved to 7b06e4441, one commit past 8c423bb52) -> merge 51c8cdac4, no conflicts, NOT pushed. Batch2 tree
  ~/SIMPLE_m8c (build_m8c.sh, gate stree_gatem8c.sbatch): m8c1 = 8c423bb52 + ours (gated), m8c2 = branch state 51c8cdac4.
  This branch is the one to move forward with; fast-tree rebases onto 51c8cdac4.
  m8c1 r1-root validation (vm8a/vm8b vs v47a/b): stage-B 0.992-0.994 identical (noise), site step 4.7 -> 2.2-3.0 min (CPU fix),
  prob_tab table fill only ~5-7 % faster (master pftc gain barely visible). Rest-rule studies (analysis only, default unchanged,
  waiting for the user): combined-z (fixes S2, admits junk-tilted simsites rests), (d) junk-screen-before-rest (same decisions),
  mixture test (calls S2 AND simsites rests "state" -> not a separator). Full-depth toys (stree49): S2 2 leaves (middle split),
  C discovered 2 leaves (domain A yes, closed no), C true 5 leaves (closed recall 0.80 / precision 0.52, domain A yes).
  DEV TWINS (reproducibility + fast-tree reference): devtwin_a 54096776 SIMPLE_SEED=1, devtwin_b 54096777 SIMPLE_SEED=2, m8c2,
  stree_devtwin.sbatch (nthr 4, nparts<=8, nodepar 10, maxdepth 20, minpop 3000, csbdevel). SIMPLE_SEED env = deterministic
  seed_rnd (else /dev/urandom).
- TODO full-f692 port run: check whether the 'arm in neither position' particles (f692st r6_s1/s2/s4/s5/s6/s7/s8, r4_s1/s6/s7,
  r2_s4, ~36k, junk per the user) land together and whether the junk rule flags them.
- Synthetic truth trees launched 2026-10-21 ~20:30 (stree_tree3.sbatch, s30a job 54056062, s55 job 54056063, NP 8, NODEPAR 1):
  comparators scored with the same scorecard: s30a rflex3 (rflex_s30a), rflex7 (s30c), rflex6split (s30split5); s55 klem chain
  synth_55_200k consolidated_12_rule_12 / rule_auto. No rflex13 script tree exists on s30a or sim55_200k.
