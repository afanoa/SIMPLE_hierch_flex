---
name: igm-like-synthetic-armsigm1-2026-10-19
description: "IgM-like synthetic armsigm1 (big asymmetric core 67 % power, 5 Fab arms 6.6 % each, RNase P absolute noise) - poses 5.1 deg, per-image arm d' 2.2 at refined poses, 2D sub-classes best 0.14; arms3 at 1.4 deg poses only 0.54"
metadata:
  node_type: memory
  type: project
  originSessionId: a718a461-105d-44de-b65c-c91952cccfae
  modified: 2026-10-19T22:03:07.499Z
---

Built 2026-10-19 because arms3 is unrepresentative of IgM: its 36 A core is too small to align on (refined poses 7.5 deg
median, 15 % > 20 deg, 4-fold pseudo-symmetry flips) and simulate_particles sets noise from EACH IMAGE'S OWN variance
(add_gauran: sdev = sqrt(var/snr)), so a heavier particle gets no SNR gain unless snr is scaled.
- Generator toy55/make_igm_maps.py (filled textured asymmetric core 54x46x40 + RNase P density; 5 textured Fab arms
  36x22x20 on the rim at 72 deg, theta +-25 / phi +-30 independent; writes snr_scale.txt = projection power vs the
  RNase P consensus, x5.80). run_igm_sim.sbatch MODE=armsigm1 (one job: maps -> sim at snr 0.0012 x 5.8 -> refine3D
  prob lp 8 mskdiam 220 from 10-deg perturbed -> old chain inline -> per-arm readout). igm_arm_masks.py (sector masks).
- Results: refined poses 5.1 deg median, p90 9.1, no flip tail. arm_llr ceiling per arm (> 30 deg far): true poses
  lp 20 d' 2.4, lp 10 1.6; REFINED poses lp 20 d' 2.2 (arms3: 2.5 true -> 0.96 refined). The data carry each arm.
- 2D route (abinitio2D 100 + cls_expansion 5 classes): 90 sub-classes, best per-arm std ratio 0.14, many 0.25-0.47,
  median 0.58 (arms set: best 0.46).
- arms3 old chain at truth+3 deg poses (1.4 deg median): best leaf 0.54-0.59 (refined poses 0.66-0.69), board corr
  0.64; phi barely moves, most regions lock onto arm 1 (2-3 regions per family for 4 arms) -> poses help modestly; the
  structure (too few per-arm regions) is the limit. Oracle per-arm board at true poses: median 0.50, corr 0.73.
- Running at write time: old chain on armsigm1 (rnasep_klem/synth_armsigm1), new dtree_v10 on halves (igm_dtree.sbatch
  -> simarmsigm1_A/_B, lat_armsigm1_A/_B, dtree_armsigm1_dtree, readout dtree_arm_readout.py).
- NEW TREE (dtree_v10 on independent halves, igm_dtree.sbatch): root K* 7 with junk isolated (19.9k / 20k junk in one
  child, both halves), then ~30 nodes per half, leaves 0.8-5k; every leaf narrows several angles a little: best per arm
  angle 0.53-0.67. A partition cannot represent 5 independent arms (10-D product space). Junk branch still splits
  (K* 5 -> 2 -> 2 on 2 junk volumes): the shuffle null may miss view sorting in junk.
- FACTORISED (igm_armboards.sbatch: one focused board per arm, sector mask, random start, refined poses, 7x3 + 1 free):
  every board tracks its OWN arm, |corr(z, angle)| on the dominant angle 0.56-0.75 (arm1 phi 0.75, arm2 theta 0.66,
  arm3 phi 0.56, arm4 theta 0.71, arm5 phi 0.68), other arms <= 0.17; per-cell std ratio best 0.43-0.60. Second angle
  weaker (one lattice axis dominates). The single free component does NOT catch junk in focused boards (purity
  0.12-0.5): clean junk first (tree root). Masks were truth-geometry sectors: a real run needs >= NARM data-derived arm
  regions. Design for IgM: tree root (junk + coarse groups) -> per-arm focused boards -> per-arm angle distributions.
- OLD CHAIN on armsigm1 (data-derived regions, refined poses): 5 motion regions per family this time; region boards
  track arms 1/3/4/5 at |corr| 0.44-0.56 (arm 2 never; 3 boards at 0.08-0.22), leaf states weak (std ratio mostly
  0.63-0.93): the continuous per-arm board coordinate is the deliverable, clustering cells into states loses it.
- TRUTH-FREE factorised (igm_factorised.sbatch -> factorised_armsigm1/): junk = dtree root child on the root board's
  free components (1.00 vs 0.00): 98.7 % of junk at 99.9 % purity. Cleaned whole-particle board (7x7, 20->10 A, random
  start, 32 cells) -> derive_regions2 k 8 -> 8 regions, most centred 16-27 A IN THE RIGID CORE; only region 6 = one arm
  (0.94 overlap), others mixed or core. 5x5 focused boards, CV R^2 of the arm angle from (z1, z2): arm 2 0.38/0.31,
  arm 1 phi 0.34, arm 3 0.25, arm 5 0.09, arm 4 none (no region). Geometry masks gave all 5 at |corr| 0.56-0.75 ->
  region derivation is the bottleneck (core variability = pose/blur sorting dominates the motion profiles).
  Next lever: keep the core out of the region step (derive_regions2 --within / consensus-core exclusion), finer-band
  board for the regions (the old chain's family feedback boards at 14->8 A found 4/5 arms).
- REGION STEP v2 (igm_factorised2.sbatch -> factorised2_armsigm1/): periphery_mask.py (consensus envelope minus the
  high-density core; core 0.40 / env 0.05 picked WITH a truth check on this set) + 14->8 A board from the coarse
  board's coordinates (14 cells) + derive_regions2 --within: 8 regions centred 44-62 A, single-arm overlaps arm 2 0.97,
  arm 3 0.84-0.93 (x3), arm 5 0.99 (x2), arm 1 0.58 (+ one 1/2 mixed); arm 4 STILL no region (also missed with the
  geometry-free v1; the geometry-mask board finds it at 0.71). CV R^2 best: arm 2 theta 0.49, arm 3 phi 0.36, arm 1
  phi 0.24, arm 5 phi 0.23, arm 4 0. Motion-pattern k-means splits one arm's two angles into separate regions.
  Next: residual-periphery discovery (board on periphery minus found regions -> regions there) and merging regions
  whose board coordinates correlate (one balanced board per arm).
- REGION STEP v3 (igm_factorised3.sbatch -> factorised3_armsigm1/): residual periphery (periphery minus found regions)
  = 99 % of arm 4's zone and nothing else; residual board -> 2 regions on arm 4 (0.88 / 0.97), R^2 0.48/0.23. Merge
  (region_merge.py, CCA of board coordinates): different arms 0.12-0.39, same arm 0.51-0.68, nuisance regression
  changes nothing; gap rule's 2-sd guard REJECTED the 0.39 -> 0.51 jump at 12 candidates (sd 0.075) -> no merge (1.5 sd
  would merge). Fresh balanced 5x5 boards on all 12 regions (new seeds): ALL FIVE ARMS truth-free, both angles:
  CV R^2 arm1 0.19/0.19, arm2 0.41/0.22, arm3 0.39/0.29, arm4 0.40/0.34, arm5 0.34/0.33 (|corr| 0.43-0.64; geometry
  masks 0.56-0.75). SEED SENSITIVITY: same mask, other seed: arm 5 0.08/0.17 -> 0.34/0.33, arm 3 0.18/0.15 -> 0.39/0.29
  -> multi-start per region (keep the most confident board, as dtree does) is the next lever.
- v4 MERGED + MULTI-START (igm_factorised4.sbatch -> factorised3_armsigm1/merged4; region_merge gap:1.5 -> threshold
  0.45, 12 -> 8 regions = arm 3 x3, arm 4 x2, arm 5 x2 merged; igm_mstart_board.sbatch array 4 starts; start_select.py):
  mean max responsibility picks the truth-best start in 6/6 informative regions (good starts maxr 0.25-0.44, failed
  0.05-0.15). Dominant-angle |corr| per arm 0.60-0.66 for ALL FIVE (CV R^2 0.45-0.53, arm 5 phi 0.41) = the geometry-mask
  boards (0.56-0.75); second angle weak (R^2 0.09-0.22; v3 single starts had 0.19-0.34 on it). 2 of 8 regions (mixed,
  no arm) track nothing (maxr 0.14-0.17): need a truth-free drop rule (cross-half reproducibility of the coordinate).
- HARDER SET armsigm2 (2026-10-20, igm2_pipeline.sbatch: make_igm_maps --pairs --splay 18 --theta 40 --phi 75 = 10 Fabs
  in 5 hinge pairs, +-58 in / +-75 out, core 50 % power, Fab 5 %; route unchanged): refined poses median 8.1 deg but
  p90 142 deg (pose FLIPS, pseudo-5-fold; armsigm1 had p90 9); junk 98.9 % at 0.96 purity; periphery keeps 82-91 % of
  arm zones; regions + residual cover all 5 hinge pairs (overlap 0.92-0.98 per pair); but EVERY board fails: per Fab CV
  R^2 <= 0.10 (|corr| <= 0.30). Suspects: pose flips, 2 Fabs x 2 angles per pair region (4-D on a 5x5 lattice), wide
  swings. Next: arm_llr ceiling at true vs refined poses; core-masked refinement (pcg_mskfile = the periphery step's
  data-derived core) as a route stage.
- armsigm2 CEILING + CORE REALIGN (igm2_realign.sbatch, 2026-10-20): whole-image conformer id (200) true lp20 acc 1.00
  d' 3.9, true lp10 0.46 / 1.57, refined lp20 0.59 / 1.42, realigned lp20 0.50 / 1.24 -> the current poses cost most of
  the evidence. Fab-SPECIFIC share is tiny: LLR vs arm-far (>30) ~= vs arm-near (<5, others free), 1-2 % (armsigm1 ~14 %),
  so arm_llr's per-arm d' is not arm-specific with 10 independent Fabs. Core-masked refine3D (pcg_mskfile = periphery
  step core 344 nm^3, prob, lp 8, 5 its) made poses WORSE: median 8.1 -> 9.8, > 60 deg 16 % -> 22 %, > 120 8 % -> 11 %;
  factorised route on realigned poses still R^2 <= 0.08. Likely cause: the synthetic core is a near-D2 ellipsoid with
  weak texture (armsigm1's arms broke the symmetry; wide-swinging Fabs cannot) -> a synthetic artefact, real IgM Fc+J is
  asymmetric. Next: route at truth+3 deg poses (separates pose loss from the per-Fab signal limit); asymmetric core.
- armsigm2 at NEAR-TRUE poses (igm2_truepose.sbatch, truth+3 deg = 1.41 deg median, eo by parity in the import oritab
  survives import; route unchanged, OSFX=_tp3): every one of the 10 Fabs tracked on one angle, mostly out-of-plane phi,
  best CV R^2 0.13-0.27 (in-plane theta <= 0.16); a pair region's board can carry both Fabs (region on pair 5: arm 9 phi
  0.26, arm 10 phi 0.22). Refined poses <= 0.10, core-realigned <= 0.07 -> poses cost 2-3x in R^2, but even near-true
  poses stay far below armsigm1 (0.41-0.53): the per-Fab-specific evidence limit (1-2 % of an image's score) is real.
  Merge did not fire (no gap, 12 regions). Implication: real IgM needs good poses (asymmetric Fc+J core) AND pooled
  evidence (2D sub-classes) for in-plane Fab angles; per-particle boards give each Fab's out-of-plane tilt.
- armsigm2 POOLED 2D ROUTE (igm2_route2d.sbatch: synth_2droute MSK 220 + rd2d stage1 print + subclass_angle_readout.py):
  99 sub-classes (median 98 ptcls, 5 expanded classes only); EVERY Fab pinned by some sub-class in BOTH angles: in-plane
  theta best 0.20-0.49 (8/10 <= 0.42), out-of-plane phi 0.16-0.37; 15 labels pin a theta < 0.5, 27 a phi < 0.5. A no-
  information label of ~40-100 particles would sit near 1 (min over 99 ~ 0.7), so this is signal. Same readout on armsigm1:
  theta 0.14-0.52, phi 0.39-0.57. 2D alignment is immune to the 3D pose flips. Per-particle boards at near-true poses
  reached only R^2 0.13-0.27 on one angle. Open: coverage (only 5 views expanded) and the 2D -> 3D link (per-Fab labels
  across views; earlier global stitching failed on the arms sets).
- armsigm2 FULL EXPANSION + BAG PRIOR (igm2_bags.sbatch, run_kernel_cls_all.sh: 97 classes -> 1941 sub-classes in 5 min;
  igm2_bagboard.sbatch: Python klem.py boards, --bags + --bag_parents, 5x5 + 1 free, per merged region): NO GAIN. Near-
  true poses, same start: beta 2 vs beta 0 best per-Fab R^2 phi 0.09-0.21 vs 0.14-0.25 (beta 0 as good or better; theta
  mixed, arm 3 0.15 vs 0.08); refined poses with bags <= 0.11 (= plain Fortran <= 0.10). Pooled log-likelihoods only
  sharpen what the cells already resolve. TOY ARTEFACT: with 200 discrete configurations some full-expansion sub-classes
  (30-43 ptcls) hold ONE configuration (std ratio 0.00 for every Fab), so pooling is optimistic here; bulk: 152/1288
  labels pin a theta < 0.5, 175 a phi < 0.5. The 2D -> 3D cross-view link stays the unsolved piece for IgM.
- armsigm2 SUB-CLASS E-STEP (igm2_avgs.sbatch / igm2_avgboard.sbatch: klem.py --avgs on the 1940 full-expansion averages,
  ptcl2D dump as bags -> class-frame renders, calibration ncc 0.84 rot sign clear; near-true poses): pool COLLAPSES (all R^2
  <= 0.02); unmix phi 0.01-0.27 / theta <= 0.08; unmixc worse; per-particle control (same Python board) phi 0.14-0.25 /
  theta <= 0.13. No mode beats per-particle; in-plane stays <= 0.1. Second failed 2D->3D transfer (bag prior, sub-class
  E-step) on this set, as on the light arms set: likelihood pooling under the board's cells cannot carry "one Fab pinned
  per view" into 3D. Next idea (not run): explicit per-Fab rigid-body angle fit to each sub-class average given its view.
- STATUS 2026-10-20 (user): IgM PAUSED as the paper's limits case. Refined-pose sub-class pool also collapsed (<= 0.02).
  Works: junk removal, finding the moving parts (periphery + regions + residual), per-arm coordinates on 5 narrow arms,
  out-of-plane Fab tilt on 10 wide Fabs at good poses. 2D only: in-plane Fab angles (pooled sub-classes). Unsolved:
  cross-view 2D -> 3D link (stitching, bag prior, sub-class E-step all failed), pose flips on the symmetric synthetic
  core, 1-2 % per-image Fab-specific evidence. Reopen with real EMPIAR-11077 (asymmetric core) or a parametric per-Fab
  rigid-body fit to sub-class averages (user's call: departs from the data-derived design).
TRAPS: background watcher loops must use `ssh -n` (without it the loop's ssh silently fails and nothing is submitted);
cls_expansion projects (stage1 build) need ~/SIMPLE_stage1/build to print (main build: unsupported particle record
width) - rd2d.sbatch. Real IgM is harder: [[igm-pdh-papers-heterogeneity-baseline]] (Fabs +-50-60 in, +-70-80 out).
Related: [[igm-arms-ceiling-pose-limited-2026-10-19]], [[klem-synthetic-benchmarks-arms-tethered]], [[klem-data-derived-tree-v3-rule-2026-10-19]].
