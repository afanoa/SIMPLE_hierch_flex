---
name: klem-organization-and-merge-up-2026-10-20
description: "f692e region boards are redundant (CCA 0.81-0.88: ~2 motions per family, not 6) -> closed conformers scattered over 4 boards; user goals: merge UP the tree, states as the particles support, sample extremes + intermediates, organization must be meaningful (composition top, conformation deeper) and general; merge-up / landscape-sampling / v2 tools"
metadata:
  node_type: memory
  type: project
  originSessionId: 25d47be4-7d1d-48d8-abdf-97b1c3c946ab
  modified: 2026-10-20T19:03:37.784Z
---

User 2026-10-20 afternoon (after f692r reproduced f692e): improve the old build with the later work; "instead of pruning
reproducible components we should merge up the tree"; 218 zoom states hold several closed conformers but too few particles
per state -> low resolution; "support as many states as particles" + "sample the extreme and intermediate states to say
anything about biology"; the organization must be MEANINGFUL (expects composition in the top layers, conformation deeper)
and must work on other datasets (no RNase P-specific masks/references; closed-arm readout only as a check).
Findings (f692r = f692e rerun):
- Per-half discovery (f692hA3/hB3) is uninformative: 2 of 63 f692e states found in both halves; half projects' region
  states eo-skewed ([[halves-split-by-eo-break-refinement-trap]]); f692e/f692r states are eo-balanced (0 of ~150 skewed).
- REDUNDANT REGION BOARDS: region_merge_check.sbatch (CCA of board coordinates) on f692r fam1: region_1, motion_1/2/3 at
  0.81-0.88, motion_5 0.69-0.73, motion_4 alone 0.17-0.36 (fam2 0.77-0.94 among four) -> ~2 motions per family. Closed
  leaves s1 of four boards share 35-56 % of particles (board_overlap.sbatch); consolidation hands each particle to one board
  -> closed conformers (0.93 in region_1 s01, 0.92 in motion_1 s01) in different branches. region_merge's gap rule merged
  nothing (it assumes most pairs unrelated); fixed 0.5 separates both RNase P (0.69+ vs <=0.36) and IgM (0.51-0.68 vs
  0.12-0.39). region_merge.py now reads OUTPARENT/region_merge_rule.txt if present (opt-in; IgM unaffected).
- MERGE UP (merge_up.sbatch -> mu_prep.py / mu_tests.sbatch (splitbench8 --pair identity --xem 3 --softem 0 --focus) /
  mu_collect.py / mu_rec.sbatch / pooled_analysis): f692r zoom level 218 -> 185 (3 blend-pruned, 15 drained by the
  cross-half EM at <= 20 A, ~15 pairwise), region level -> 157 (both closed maps still separate). mu_family.sbatch (all units
  of a family, periphery focus, --xem 0 --arch 0 = no re-assignment/blend pruning): fam1 61 -> 30.
- LANDSCAPE SAMPLING (landscape_sample.py/.sbatch): leaf maps -> deviation from the family parent in the periphery, PCA,
  fold the smallest group into its nearest until every group >= NMIN; NMIN from the leaves' periphery resolution vs ln N at
  TARGET A. User: only meaningful if the tree organization is (pure leaves, one conformation per branch).
- v2 route (routes/toy55_f692v2 = guarded copy, csbdevel): f692_regions_v2.sh = periphery (periphery_mask.py) -> candidate
  regions (derive_regions2 --within) + boards -> residual periphery regions + boards -> region_merge (RMTHR 0.5) -> 4 random
  starts per merged region (keep max mean maxr) -> f692_region_one3; f692v2_chain.sbatch (SRC run's family stage linked).
  Test: f692s55v2 (synthetic 55, SRC f692s55d, CTFABS 0 PSCALE 0), job 54027289. Synthetic halves: sim55_700k_A = odd
  rows, _B = even rows (meta_rows.txt written). f692s55d zoom consolidation only 9 states (1/55 strict) - too thin to calibrate.
Next agreed shape: composition at the top, ONE board per distinct motion (merged regions), zooms for purity inside, merge up
+ sampling on top; measure organization by leaf purity and branch coherence (synthetic truth; real: closed proxies).
Related: [[klem-tree-vs-product-organisation-2026-10-18]], [[klem-leaf-states-are-redundant-across-regions]], [[f692e-landscape-pruning-test-2026-10-20]].
User 2026-10-20 15:10: drop composition (core borrowing) - "if we merge correctly it should work"; deliver plain all-particle
maps of the merged states (states_dir). The masked-periphery FSC budget failed (median 2k-particle leaf "4.3 A" -> no
pooling): use fixed particle budgets (e.g. 10k / 20k) for landscape sampling. Family-level merge-up f692r: 218 -> 84
states (median 3.8k, max 28k particles); the two closed maps stay distinct (15: 0.921, 6.1k; 16: 0.931, 5.5k).
User 15:12: "although weighted is fine" - soft/weighted state reconstructions are OK (not composition): per-particle posterior weights over the final merged states (half A weights from half-B maps), via SIMPLE state_weight_set + reconstruct3D; report native weighted half-map FSC.
GOAL (user 2026-10-20 18:20, work autonomously): "a tree organization which actually represents the heterogeneity so any merge will
be legitimate". Fast-loop findings on the synthetic benchmarks (simswing30 / simswingbreath / simprodi120 / simprodc120, truth in
gt_particles.txt): (1) on a mixed pool region boards sort COMPOSITION, not motions (per-bin true swing flat at 15 A); (2) oracle
composition (swing family only) -> one factor tracks the swing but compressed (6.5-23 of 0-30 A) + spurious noise factors;
(3) the free component on synthetic boards swallowed 57-64 % of real particles and 1 % of junk (real f692r boards: 7-11 %) -> NFREEB 0;
(4) composition is NOT separated by the board families: f692sw30 parent (full chain) families 42-48 % apo each; cells only 60 %
pure (chance 41 %) -> next: purify composition by seeded K-state refine3D_states (fast_comp.sbatch), then factors per class,
then focused (solvent/pcg_mskfile masked) refine3D_states along each factor to sharpen ends. User: refine3D already has a
solvent mask and pcg_mskfile, so focused refinement is easy.
HIGH-SNR fast tests (SNR 0.005 = 4x the calibrated 0.0012, 100k; sims swing30hs/swingbreathhs/prodihs/prodchs; ~6 min each):
- realistic SNR at 100k is BELOW the data's support (families ~30k each; f692e worked at ~300k per family): every set failed there.
- product sets: the dominant class separates -> one factor tracks the swing (Spearman +0.71/+0.77; bins 7.0 -> 24.1 A of 0-30, breathing
  flat); COUPLED set -> one factor carries both (k 7.8 -> 23.7 A with j 2.8 -> 8.0 %) = the right answer; independent: breathing factor
  MISSING (|rho| <= 0.09 on both axes of every region board and of a whole-particle family board) although breathing maps differ more
  than swing maps (8 A corr 0.90 vs 0.92) and 3D classification separates the breath family -> a board limitation, open.
- composition: board cells -> spectral seeds (eigengap K) -> refine3D_states re-classification: purity 0.41-0.72 -> 0.48-0.84; with
  KSEED=8 over-clustering: 0.50 -> 0.86 (swing30hs), 0.39 -> 0.86 (swingbreathhs): apo / holo / swing family / breath family / junk
  each get their own states (apo in 2-3 = its conformers). refine3D_states takes pose_policy local|global only (not fixed).
- "composition" mixing on the swing sets was partly the ARM: with each particle's map closed-arm occupancy as truth, factors order
  closure across apo (open) and swing-family (closed) particles (0.18 -> 0.62); within the closed family the 0-30 A swing stays flat.
- the representative board = max mean maxr can pick a less specific board (whole-particle) -> lower correlation; revisit.
Next: factors inside the K=8 refined states (ff3_*), reproducibility filter, focused refine3D_states sharpening.
