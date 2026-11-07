---
name: toy10-latent-night-2026-10-15
description: "Overnight 2026-10-15/16 latent-variable battery on the RNase P toy10 (and toy10c): what runs where (batch2 routes/toy55: model_class_s3/s6/c3, lat/, latc/), the scripts (model_class, ptcl_latent_eval, latent_windows, deform_proto, run_toy10_battery), the readouts, and the first model-class finding"
metadata:
  type: project
---

Plan agreed by the user ("go ahead, run it all overnight", 2026-10-15 ~03:00 EDT): decide which latent-variable
model class can represent the toy landscape and whether it can be inferred per particle at RNase P SNR. No discrete
refine3D_states from random seeds, no paper-workflow replay, no region/multibody basis (the deformation field is global).

Streams (batch2 `/mnt/beegfs/elmlund/afan/routes/toy55/`, scripts mirrored in `scripts/cls_expansion_proto/toy55/`):
- A `run_model_class.sbatch` / `run_model_class2.sbatch` -> `model_class_s3`, `model_class_s6` (10 maps),
  `model_class_c3` (12 maps with the variants): `model_class.py` = leave-one-out linear span (K=1,2,3,5; K=all is
  exact by construction since the consensus is the weighted mean) vs multi-scale symmetric demons deformation vs
  deformation + occupancy, unexplained fraction per band 117-30/30-20/20-12/12-8 A; flow-field PCA effective dimension.
  First pass (sigma_elastic 8 and 14 px) over-regularised: displacements 1-2 A, deform left 0.69-0.74 unexplained,
  linear K=5 0.64-0.67, deform+occ 0.62-0.69; volume-difference PCA effective dimension 3.5-4.7 per band.
- B `run_toy10_battery.sbatch` (hooked by the sim) -> `lat/` (toy10) and `latc/` (toy10c, env TAG=c): 0 oracle
  E-steps vs the TRUTH maps at lp 20/14/10/8 (40k subset; the per-particle ceiling); 1 flex_pca lp 20 neigs 16 gmm
  placement + `ptcl_latent_eval.py` (Spearman vs the conformer axis per PC and cross-fitted, 15-NN accuracy, pairwise
  d', closed AUC, delivered states, junk AUC); 2 `latent_windows.py` rounds at 20/14/10 A: assign (8 quantile
  windows of PC1) -> reconstruct3D nstates=8 -> maps readout (FSC, corr with the dominant truth conformer per band)
  -> rounds 2+ per-window refine3D refine=shc maxits 2 (time-guarded) and merge -> soft cross-half E-step with a
  Gaussian smoothness prior (sigma 1.5 windows) -> posterior-mean coordinate; 3 `deform_proto.py` on the final
  window maps (demons to the middle window, flow PCA, flow-corrected consensus FSC, grid E-step on 20k vs windows);
  4 refine3D_states nstates=4 from the coarse windows (pose_policy local, 20->8 A) + flex_pca per coarse group
  (mixture of local linear models) scored on the same particles as the global latent; 5 solve2D + pose-mean
  cls_expansion -> sub-class-pooled readout of the particle latent.
- Conventions: E-step projections use `view_tangent_toy` helpers with TOY_SMPD=1.464, SIMPLE slice convention
  euler2m([-e1,e2,-e3]) and an automatic calibration (conv, shift sign/scale, contrast sign) per run; `selection`
  stores integer window labels from `infile`; `--axis=$AXIS` (values start with '-').
- Jobs 2026-10-15: sims 53878367 (toy10) and 53878368 (toy10c), model-class rerun 53878369; waiters in the session.

Results so far (2026-10-15 05:50 EDT): toy10 = 199,996 particles (160k conformer, 40k junk, 22.5k closed), refined
poses median 3.98 deg (p90 7.35), shifts 1.44 px. flex_pca lp20 neigs16 on the particles: Spearman of every PC with
the conformer axis |rho| <= 0.20, cross-fitted supervised combination 0.49, 15-NN accuracy 0.167 (chance 0.149),
pairwise d' median 0.88, closed AUC 0.67, junk AUC 0.81 (PC1 is junk/contrast), delivered states 1 giant + 3 tiny:
the per-particle LINEAR latent is near chance at RNase P SNR with 4 deg poses. PC1 windows = random partition
(purity 0.15, window maps 5.3-7.1 A all consensus-like). A supervised-seed arm (lat_sup, SEEDCOL=2 = cross-fitted
truth axis, shared flex run) tests the window machinery independent of the seed.
TRAP fixed 05:46: the E-step must fit a per-shell amplitude scale of the reference projections to the images
(external maps / SIMPLE recvols are in arbitrary units); without it the likelihood ranks references by their own
power ("signal fraction" 1.006 in every shell, max posterior 1.0, everything assigned to the two smoothest maps).
The first oracle numbers in toy10_battery.53879061.out are from the unscaled code; the valid oracle is the
toy10_oracle job (lat/oracle_lp*.txt rewritten).

ORACLE (corrected, toy10_oracle.53879242.out, 40k-particle subset, truth conformer maps as references, cross-half not
needed): per-particle accuracy over 10 conformers 0.53 / 0.57 / 0.58 / 0.59 at lp 20 / 14 / 10 / 8 A (chance 0.15),
Spearman(posterior-mean axis, truth) 0.81 / 0.83 / 0.84 / 0.84, closed recall 0.70-0.77, mean max posterior 0.8;
signal fraction per shell 0.4 (lowest) -> 0.23 (20 A) -> 0.05 (12 A). => per-particle state inference IS possible at
RNase P SNR with 4 deg poses when the references are right; most of the information sits at <= 20 A. The linear
flex latent (0.17 accuracy, 0.49 supervised Spearman) is the bottleneck, not the data.
Window E-step from the PC1 seed (random windows): purity 0.156, Spearman 0.0 at lp 20 (consensus-like references
cannot bootstrap). Arms added: lat_sup2 (supervised seed, rank-based windows, junk coordinates drawn at random) and
lat/pooled (toy10_pooled job: K=8 k-means on the pose-mean sub-class-POOLED latent -> hard EM, waits for exppm).

EM CLIMBS from a weak seed (lat_sup2, 05:58): supervised-seed windows Spearman 0.485 / purity 0.208 -> after ONE
cross-half E-step against the 8 window reconstructions at lp 20: posterior-mean coordinate Spearman 0.722 (flat
prior) / 0.639 (smoothness prior sigma 1.5), purity 0.263, Hungarian 0.228; oracle ceiling 0.81-0.84. The flat-prior
posterior beats the smoothed one, so the carried coordinate is now the flat-prior one (latent_windows.py patched
06:05; earlier-started processes still carry the smooth one). Rounds 2-3 (14, 10 A, per-window pose refinement)
decide whether it converges to the ceiling.

toy10c (composition, latc, 06:30): oracle accuracy over 12 states 0.49 / 0.52 / 0.53 / 0.53 at lp 20/14/10/8 (chance
0.12), Spearman 0.78-0.81, closed 0.71-0.77; the deleted-domain state 61 recall 0.56-0.58, the added-blob state 62
(2.6 % mass) 0.41-0.43 = the hardest state. flex lp20: 15-NN 0.144 (chance 0.122), supervised combination 0.46;
PC1 puts the deleted-domain state at its extreme (composition IS linear in volume space) but junk shares that end.
Arms: latc (PC1 seed), latc_sup2 (supervised seed, job 53879827), latc/pooled (toy10c_pooled 53879325).

TRAP (06:45): per-window pose refinement inside the latent loop (refine3D refine=shc maxits 2 at lp 14 against a
25k-particle window map) made the poses WORSE on the toy: median 3.98 -> 6.54 deg (p90 7.35 -> 16.6), shifts 1.44
-> 2.09 px; at lp 10 it recovered only to 5.16 deg. The PC1-seeded arm stayed at chance through 3 rounds (purity
0.16). Arms without the pose step: lat_sup3 (toy10, 53880351) and latc_sup2 (toy10c, 53880362, REFINE_WINDOWS=0);
lat_sup2 keeps it (comparison). Refinement against window references needs a tight local bound / more particles.

lat_sup2 rounds (06:52): posterior-mean Spearman 0.72 -> 0.74 -> 0.75 (lp 20/14/10), window purity 0.26 -> 0.29 ->
0.30, Hungarian 0.23 -> 0.26 -> 0.28, max posterior 0.49 -> 0.81; poses 3.98 -> 4.69 -> 4.25 deg (refinement hurts
less with coherent windows). Slow convergence below the oracle (0.84 / 0.53-0.59): the 8 window references are
mixtures along ONE coordinate while the conformers span a 2-3-D deformation space; next lever = a 2-D latent
(flow mode 1 x mode 2 grid of references) and more / softer windows.

DEFORMATION PROTOTYPE on lat_sup2 windows (07:00): demons of the 8 window maps to the middle one: displacements p90
0-5 A, flow PCA effective dimension 1.3 (mode 1 = 87 %), mode-1 score monotone with the window index (Spearman 1.00).
Grid E-step (ONE canonical map warped along mode 1 at 9 latent values, 20k subset, lp 12): Spearman 0.746, purity
0.275, closed recall 0.42 vs the 8 reconstructed window references 0.752 / 0.304 / 0.22: a single map + one flow
mode reproduces the discriminative content of 8 reconstructions and extrapolates to the closed state. Flow-corrected
consensus FSC 4.4 -> 4.2 A (windows are mixtures, little to correct). 2-D grid (mode 1 x mode 2, 35 refs, 30k):
mode 2 of windows made along ONE coordinate is noise (Spearman(z2, truth) 0.22, d' 1.42 vs 1.40): a 2nd latent
dimension must come from a 2-D state space in the EM itself, not from 1-D windows. Per-conformer z1 means order the
10 conformers into the truth groups (27/28 > 48/49 > 1/51/18/46 > 10/37).
TRAP: refine3D_states in the kernel build rejects vol1..N ("use classify3D_refs"); it derives the state volumes from
the labels -> call it without vols (battery patched; standalone run_r3s.sbatch for finished arms: lat, lat_sup2).

lat_sup3 (supervised seed, NO per-window pose step, 07:10): Spearman 0.72 -> 0.77 -> 0.76 (free) / 0.76 (smooth),
purity 0.30 -> 0.31, Hungarian 0.27 -> 0.28, 4/8 window maps match their dominant conformer at 20-8 A (sup2 with the
pose step: 1/8, Spearman 0.75). Both plateau at ~0.75-0.77: the 1-D window limit, not the pose step.

latc_sup2 (compositional toy, supervised seed, no pose step, 07:12): Spearman 0.68 -> 0.73 -> 0.72, purity 0.17 ->
0.25, Hungarian 0.23; final windows: the deleted-domain state 61 dominates window 7 (share 0.31, next to its parent
49 in window 6), conformer 27 window 8 (0.35), junk collects in window 1 (65 %): the 1-D coordinate partly separates
the compositional state because an 8 % mass loss is a strong low-resolution signal.

Deformation grid vs reconstructed windows, replicated (07:20): lat_sup3 grid 0.767 / purity 0.30 / closed 0.37 vs
windows 0.758 / 0.31 / 0.28; latc_sup2 grid 0.741 / 0.24 / 0.36 vs windows 0.719 / 0.25 / 0.21. Flow PCA on every
arm: mode 1 = 85-88 %, effective dimension 1.3-1.4, monotone with the coordinate. Standalone refine3D_states jobs:
lat 53880893, lat_sup2 53880892, lat_sup3 53882136, latc_sup2 53882147.

latc (compositional toy, PC1 seed, with pose step, 07:22): Spearman 0.21-0.22 (PC1 carries a little composition
signal), purity 0.16; the pose step again degraded poses (3.97 -> 7.08 -> 5.68 deg, p90 20 deg).

2-D grid on lat_sup3 (07:35, 30k subset, lp 12): 1-D grid 0.771 / purity 0.30 / closed 0.48; 2-D grid (7x5) z1 0.767,
z2 uncorrelated with the main axis (0.06) but it separates {1, 48, 49, 46, 51} (z2 +20..+35) from {10, 18, 27, 28,
37} (~0): pairwise d' median 1.57 in 2-D vs 1.33 in 1-D. A modest second axis exists even in 1-D-window flows when
the windows are cleaner (sup3), none when they are noisier (sup2).

POOLED 2D readout (lat, 07:13): pose-mean cls_expansion = 1600 sub-classes (median 54 members), majority-conformer
purity 0.27 (936 with >= 20 members); junk enriched x49.6 (84 % of junk in junk-majority sub-classes). The flex lp20
latent POOLED per sub-class: Spearman with the truth axis -0.42/+0.45/-0.38 on PCs 2-4 (PC1 = junk), 5-NN accuracy
vs the majority conformer 0.342 (per particle 0.167), pairwise d' 1.74 (per particle 0.88): pooling doubles the
linear latent's discriminability but it stays far from the per-particle oracle (0.53-0.59).
Pooled-seed EM (lat/pooled, step 7): K=8 k-means on the pooled latent -> seed purity 0.19 (Hungarian 0.185, closed
recall 0.49 in one state); after one E-step at lp 20: purity 0.21, Hungarian 0.20, closed 0.37 (states unordered, no
Spearman). Weak seed, weak climb: the k-means states follow view/junk structure, unlike the supervised axis.

FINER LINEAR LATENT (lat/flex12, lp 12, neigs 16, 07:31): PC5 alone Spearman +0.63 with the truth axis (PC2/3 -0.41/
-0.45; PC1 still junk), cross-fitted combination 0.713 (lp 20: 0.49), 15-NN 0.254 (lp 20: 0.167; chance 0.149),
pairwise d' 1.54 (0.88), closed AUC 0.75 (0.67), ladder 18 10 37 1 46 51 27 28 49 48 ~ truth order. The lp 20 basis
was simply too coarse: the discriminative information is at 20-12 A (where the oracle gains too). Truth-free arms
launched 07:55: k-means K=8 on the lp 12 latent -> hard EM (lat/kseed, job 53883132) and the same on lp 20
(lat/kseed20, 53883133); per-group local flex evals are re-run inside the kseed job (the battery's printed nothing).

Pooled-seed hard EM, 3 rounds (lat/pooled, 07:35): purity 0.19 -> 0.21 -> 0.22 -> 0.23, Hungarian 0.185 -> 0.20 ->
0.215 -> 0.23, closed recall up to 0.57 (r2): a truth-free seed from the 2D sub-class-pooled latent climbs slowly
but monotonically (the PC1 seed never leaves 0.15).

refine3D_states K=4 from the PC1 (random) coarse windows, local poses 20->8 A, 6 its (lat/r3s, 07:40): Hungarian
0.085 -> 0.113, weighted purity 0.158 -> 0.191, closed recall 0.50 in one state at purity 0.16, poses 3.98 -> 4.39
deg: discrete EM from a random seed moves little and mostly on junk/view structure.

refine3D_states K=4 from the SUPERVISED coarse windows (lat_sup2/r3s, 07:45): Hungarian 0.166 -> 0.179, weighted
purity 0.254 -> 0.254 (flat), closed recall 0.34 -> 0.54 at purity 0.19; poses IMPROVED 4.25 -> 3.72 deg (p90 7.1),
shifts 1.31 px: the full local multi-state refinement (6 its, 20->8 A) fixes poses when the states are coherent,
unlike the 2-iteration per-window shc step, but it does not sharpen the partition.

refine3D_states K=4 from the lat_sup3 coarse windows (07:51): purity 0.268 -> 0.285, Hungarian 0.172 -> 0.188, closed
recall 0.47 -> 0.60 at purity 0.20, poses 3.98 -> 3.52 deg (p90 6.6): best poses of the night; the partition gains little.

refine3D_states K=4 on the compositional toy from the latc_sup2 coarse windows (07:51): purity 0.209 -> 0.185,
Hungarian 0.133 -> 0.139, closed 0.45 -> 0.58, poses 3.97 -> 3.58 deg: poses improve, the partition does not.

TRUTH-FREE k-means seeds + hard EM (07:59): lp 12 latent K=8 (lat/kseed): purity 0.220 -> 0.236 -> 0.263 -> 0.280,
Hungarian 0.213 -> 0.276, closed recall 0.35 -> 0.53 in one state (purity 0.33); lp 20 latent (lat/kseed20): purity
0.180 -> 0.208 -> 0.243 -> 0.272, Hungarian 0.15 -> 0.26. Both still rising at round 3 and already at the level of the
supervised-seed window arms (0.30-0.31): the cross-half E-step does the work once the seed is not a single junk axis.
Continuation of lat/kseed for 4 more rounds at 10/8/8/8 A submitted (toy10_kseed_cont, KSTART=3). The per-group
"local flex" runs never produced an embedding (flex_g1..4 lack flex_pca_embedding.bin) - not pursued.

toy10c tail (08:14): lp 12 latent 15-NN 0.210 (chance 0.122), cross-fitted 0.70, d' 1.59, closed AUC 0.73; PC1 at
lp 12 puts the deleted-domain state 61 at one extreme (composition + motion mixed). Pooled 2D: sub-class purity 0.23,
pooled 5-NN 0.32, d' 1.99. Pooled-seed hard EM: purity 0.16 -> 0.19 -> 0.20 -> 0.21, Hungarian 0.21, closed 0.50-0.53.
toy10c k-means EM on latc/flex12 submitted (toy10c_kseed).

toy10c k-means EM on the lp 12 latent (latc/kseed, 08:36, truth-free): purity 0.205 -> 0.206 -> 0.226 -> 0.239
(chance 0.12), Hungarian 0.24; final states: the DELETED-DOMAIN state 61 owns a state at share 0.37, two states are
pure junk (95-98 %), 46 / 28 / 49 / 37 dominate others at 0.18-0.27: the truth-free EM isolates junk and the
compositional state by itself; the conformational states remain mixtures.

k-means EM continuation on toy10 (08:47, rounds 4-7 at 10/8/8/8 A): purity 0.286 / 0.282 / 0.286 / 0.286, Hungarian
0.283, closed recall 0.56 at purity 0.41 (the closed state crystallises in one of the 8 states). PLATEAU ~0.29 with 8
hard states for a 10-conformer continuum + junk. Night complete 08:47; report in
scripts/cls_expansion_proto/toy55/NIGHT_REPORT_2026-10-15.md.

**How to apply:** read `lat/` logs (toy10_battery.*.out) before any RNase P run; the morning deliverable is one table
per step with the truth readouts and the model-class verdict. Related: [[rnasep-toy55-design-2026-10-15]],
[[no-region-basis-for-rnasep]], [[prototype-in-python-before-fortran]].
