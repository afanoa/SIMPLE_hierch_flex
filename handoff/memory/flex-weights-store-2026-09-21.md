---
name: flex-weights-store-2026-09-21
description: "flex per-state weight files (flex_weights_state_NNN.bin, out-segment imgkind flex_weights, sigma2-store pattern) implemented on master 2026-09-21; batch2 test tree ~/SIMPLE_flexw; the Epanechnikov kernel still has live callers"
metadata: 
  node_type: memory
  type: project
  originSessionId: f68e7939-b211-41bf-a8cf-50e0616a32a4
  modified: 2026-09-21T16:15:43.324Z
---

2026-09-21: the flex state-weight store landed on `master` (base 2b6f6c833), uncommitted at the time
of writing. ONE FILE PER STATE (`flex_weights_state_NNN.bin`): the user argued for it because the
workflows select states, and `vol_flex` is already one out-segment entry per state. My first cut was
one multi-column file plus a projinfo key; both were replaced. Files: `src/fileio/simple_flex_weights_file.f90`
(bytes), `src/main/flex/simple_flex_weights_state.f90` (domain, producer `flex_weights_deliver`,
consumers `flex_weights_consumable/load_state/load_all`), `sp_project%add_flex_weights2os_out` /
`get_flex_weights` (imgkind `flex_weights`, state s, in `simple_sp_project_out.f90`; also in the
state-associated removal list), and the producer call in `simple_flex_pca_model.f90` (before
`write_discrete_state_project`, NOT gated by rec_states). The user objected to my touching the sigma2
registration routines ("I didn't ask you to touch sigma file stuff"), so sigma2 code is untouched.
Unit tests are written, run on batch2, then DELETED at the user's request ("we know it works").
Plan + part 2 (consensus volume, NOT done):
`doc/implementation_notes/flex_weights_store_and_canonical_consensus_plan_2026_09_21.md`.

Format follows the canonical sigma2 store AFTER e647a3705 (2026-09-21: sigma2 dropped its per-row
checksums, writes rows contiguously): per file header(512) | scalars real64 (mass, neff, pop, bw,
target) | weight real32 (nptcls) | hard-label flag int32 (nptcls); full physical layout, layout digest
= `sigma2_state_project_layout_digest`; the files of one delivery share generation/digest/nstates and
are validated as a set before any is published. `weight_kind` is inferred (partition = every nonzero
row sums to one across files; the default GMM AUTO path is a partition, kernel rows are not).

**Kernel verdict:** the Epanechnikov kernel is NOT dead code. It is delivered on the reliability-path
fallback of the DEFAULT placement (when diffusion k-center fails), on `state_axis>0`, on external
targets, in `cv_select_bandwidths` (`nbins>1`, default 1) and under `SIMPLE_COV_GMM=0`; its `dist_out`
feeds the bandwidth CV. Removing it means removing those paths first; the user asked to remove it only
"if it truly has no callers", so it stayed.

**How to test on batch2:** tree `~/SIMPLE_flexw` (worktree of `~/SIMPLE` at origin/master + the patch,
Release/ARCHOPT/tests build in `~/SIMPLE_flexw/build`, built on the login node; `cmake .` before make
when sources are added). End-to-end: `~/flexw_resume.sh` = states-only resume (`rec_states=no`,
~5 min incl. the sigma2 bootstrap) from the 10180 `S_P_pcg` embedding cache into
`/mnt/beegfs/elmlund/afan/10180/flex_verify_20260915/W_store`, run via `srun -p csbdevel -c 8`;
verify with `print_project_field oritype=out` and a python read of the five files against
`flex_pca_state_weights.txt` / `flex_pca_coordinates.txt` (labels).
**Consensus experiment (part 2, 2026-09-21):** `src/main/flex/simple_flex_pca_consensus.f90` builds seven
candidate weightings (plain, dominant, viewbal, medoid, density, chi2rob, dom_rob) under
`SIMPLE_COV_CONSENSUS=1` and reconstructs them as `flex_pca_consensus_NNN.mrc` (+even/odd, NOT registered;
`reconstruct_flex_weighted_states` gained `register=`). Run: `~/flexw_consensus.sbatch` -> `~/flexw_consensus.sh`
(10180 S_P_pcg resume, rec_states=yes, box_rec=320, nparts=10) into `.../10180/flex_verify_20260915/W_cons`;
score with `~/flexw_score2.py` ALWAYS through `srun -p csbdevel -c 4 --mem=32G` (box-320 FFTs; the user
objected to it running on the login node: "it's just going to get killed"); numpy 1.14 on batch2, so no
`np.fromfile(offset=)`. RESULT 2026-09-21: eo-FSC favours N monotonically (plain 6.7 A beats every
subset); chi2-residual robust weight is blind to the junk the latent isolates; AND the run's frame was
wrong (3_selection poses vs 4_abinitio3D vol1: FSC 0.01), so the embedding must be refitted on the
4_abinitio3D project before the weighting question is answered (plan doc section 2.5). Deviation-space scorer `~/flexw_score3.py <rundir> <target_state>` (mean-subtracted FSC per RECOVAR SI S.J + our
merge gate; local res p90 needs win=64): on the first run only the dominant-state family moves the map toward the
target (dev FSC 0.96-1.0), medoid/density move elsewhere, chi2rob moves nowhere (dev res 68 A). REFIT on the
4_abinitio3D project (`.../10180/flex_verify_20260915/W_refit`, 25 min wall for fit+states+candidates) VERDICT:
the random matched-count control reaches 8.49 A vs the dominant state 8.63 A, so selecting a conformation buys NO
global resolution, only specificity (dominant deviation reproducible to 21 A, random's is noise); medoid/density
agree with the dominant direction at 0.55-0.67; chi2rob does nothing on either frame. That subset conclusion was OVERRULED by the user's design: "I want a canonical
consensus that's a weighted combination of ALL the particles instead of 1.0 weights. The mass doesn't make sense
to me" and "the idea is that you can refine this to get better poses than an equally weighted volume." So:
full-support weights in (0,1] (posterior-metric Gaussian, population Gaussians at 1x/2x median, component-covariance
Gaussian, smooth kNN density, medoid-centred), z* = dominant population's latent mean; the decisive metric is a
PAIRED refine3D from the plain vs the weighted map, not a map statistic. Then the user objected to Gaussians ("a canonical point could be in the middle of nowhere and the gaussian
would drop off before it hit particles"): v3 candidates = mean-shift mode + kNN-adaptive kernel, Cauchy/Huber IRLS
location M-estimators (centre = fixed point of its own weights), quantile (rank) weights, posterior-metric at the
mode, plus the fixed Gaussian as the collapse reference; each centre logs a void diagnostic (kNN radius vs median).
Resume dirs `.../W_refit2` (soft Gaussians) and `.../W_refit3` (v3) for 10180 and 10028; scorer
`~/flexw_score4.py <dir> <target>` (per-shell deviation vs resolution kept); pose test `~/flexw_pose.sbatch`
(env ARM, VOL, MAXITS; the dataset's own refine3D recipe: refine=greedy lp=8 trs=10 mskdiam=500 filt_mode=nonuniform
sigma_est=global; ribo variant `~/flexw_pose_ribo.sbatch`, mskdiam=320). RESULTS so far: on 10180 (7.5 A, 47k) the
soft weightings move the map only 5-11% and their deviation is weak (0.5-0.6 repro at 40 A); on 10028 (3.4 A, 105k)
every soft weighting keeps 3.37 A AND its deviation reproduces to 9 A (better than the hard state's 13 A) = the
design works at high SNR. The soft weights on 10028 emphasise states 4-5-6-7 (the density mode), NOT the largest
GMM state (2); read direction against the state nearest the mode (`~/flexw_state_emphasis.py <dir>`). RESULT (run 7, 10180): the
responsibility floor ladder w = eps + (1-eps) r(state densest at the mean-shift mode) keeps the hard state's full
deviation quality (repro 0.80-0.82 at 40 A) at every eps while the resolution cost vanishes at eps=0.5 (7.45 A =
plain, map moved 8%); it beats the mode kernel at matched resolution (16%/0.82 vs 11%/0.57 at 7.88 A). The
canonical state must be chosen by KERNEL EMPHASIS (densest at the mode), not nearest latent mean (that picked a
diffuse junk state on 10180). Pose arms in `.../W_pose/{plain,dominant,mshift,quantile,eps25,eps50}`; compare with
`~/flexw_pose_compare.sh <W_pose>` (refine3D log summary lines).
Second specimen (10028 ribo, `.../10028_simple/flex_verify_20260915/W_refit`) submitted to confirm. The 3D wiring is NOT a loop: the user's
design (2026-09-22) is flex TWICE, at the start and toward the end, weights frozen in between, refine3D
reconstructing with the frozen weights ("we don't need a loop, too expensive; marginalising has been done, we are
trying something different"). Plan doc sections 3 (stages A, B), 4 and 5. STAGE A IMPLEMENTED 2026-09-22 (uncommitted): `simple_flex_pca_consensus.f90` =
`build_canonical_weights` (mean-shift mode, state by kernel emphasis, w = eps+(1-eps) r, CONSENSUS_EPS=0.5 internal);
the state stage writes `flex_weights_consensus.bin` (imgkind flex_weights_cons) and `flex_pca_consensus.mrc` (+halves,
imgkind vol_flex_cons); the rec3D `register=` switch and `flex_weights_deliver_column` support it. Built in the SECOND
batch2 tree `~/SIMPLE_flexw2` (the first served running jobs); verification run `.../10180/.../W_prod`. Stage B (weights
into refine3D reconstruction) not started.
Related: [[flex-pca-state-weights-are-gmm-not-kernel]], [[build-on-batch2-login-node]].

**Family rule (2026-09-22, ~00:30):** on the 10028 continuum the single canonical tile (state 5) reproduces its
departure from plain only to 19 A while the 10%-mode kernel reproduces to 8.8 A, so `build_canonical_weights` now takes
the deconvolution's macro-cluster labels (`macro=deconv_labels`, present on cache resumes via the legacy
`flex_pca_deconv_labels.txt` next to infile) and the canonical population is the canonical state's whole macro-cluster
(state -> macro by responsibility-weighted vote; r_i = min(1, sum of the family's responsibilities)); on discrete data
(every state its own macro) it reduces to the single state. Logged as `canonical population = states 4 5 6 7`; the
consensus flag column marks family members. Verification run `.../10028_simple/flex_verify_20260915/W_fam` (job
53287169, build bfc97d021afd in `~/SIMPLE_flexw2`), scored on the node by `~/flexw_fam_score.py <famdir> <refdir> [ctrl]`
against the ribosome ladder in `W_refit6` (ribo6s; `W_refit7` was cancelled at 11 min and has NO candidate maps) plain/eps25/eps50/ms10 (targets: family states, their mass-weighted union, controls 2 and 8).
Pose test on 10180 (`W_pose`, 6 arms x 3 iterations, `~/flexw_pose_compare.sh`): resolution is shell-quantised
(8.36/8.50/8.63 A) and every arm sits at 8.50 A with mean particle score 0.373 at iteration 3 = NO measurable pose
gain from any weighted reference at 10180's SNR; final table in `~/flexw_pose_final.txt` (watcher).
**Worktree caveat:** the user switched `/Users/atifao/SIMPLE` to `flex2D` mid-session (they commit there every ~10 min);
the uncommitted consensus files ride along in the worktree. Commit them on master only, and note master now has the
`Refactoring of sp_project` commit (c057c4d74, adds `add_entry2os_out`) that the batch2 trees (08bd48202) predate. The
Mac `build/modules` set is stale (Sep 18), so local -fsyntax-only of anything using sp_project fails with a vtype
mismatch; the batch2 build is the check.
**Stage B (2026-09-22 01:00, uncommitted):** the user rejected the inject-the-map pose test ("not propagating the
weighted reconstruction is probably not the best way to test this"): from iteration 2 refine3D rebuilt an equal-weight
reference, so all arms converged. Now `flex_weights=yes|no` on refine3D/reconstruct3D (params + parse + phases
`l_flex_weights` + both UI files) makes `simple_matcher_3Drec` (calc_3Drec AND calc_projdir3Drec) load the committed
`flex_weights_cons` column once via `flex_weights_load_consensus` (store module; consumers gained optional
`which_imgkind`) and scale each particle's Fourier numerator+CTF^2 plane in `prep_imgs4rec(..., scales=)` before
insertion; the search is unweighted. Verify the key reached `build/lib/simple/simple_args.f90` (args-generator trap).
Propagation test launcher `~/flexw_prop_ribo.sbatch` (env ARM VOL FLEXW MAXITS; project = W_fam/p.simple whose out
segment points at W_fam's column; ends with an equal-weight reconstruct3D in `eqrec/` = pose-quality comparison);
arms in `.../10028_simple/flex_verify_20260915/W_prop/{plain,weighted,wplain}`. The family column on 10028 came out as
states {4,5} only (K=6 macro-clusters split the continuum): dev res 16.1 A vs 19.3 (state 5) vs 8.8 (mode kernel), so
the macro-cluster family is the wrong instrument on a continuum; a regime rule (island -> responsibility column,
continuum -> mode kernel) is the reading, not yet implemented (plan section 2.17).

**Propagation test trap (2026-09-22 01:50):** flex_pca's `write_discrete_state_project` leaves ptcl3D multi-state
(labels 1..nstates), so refine3D on that project dies in every worker with "previous best state outside boundary"
(ERROR STOP inside an OMP region: the message only shows with GFORTRAN_UNBUFFERED_ALL=1 when running the worker
command under srun by hand). Fix in the launcher: `selection projfile=p.simple oritype=ptcl3D states=1,...,8`
merges them to state 1 and keeps the *_cons out entries. Stage B verified: workers log
`RECONSTRUCTION weighted by the flex canonical consensus: neff 94936.1 floor 0.50 min w 0.500 max w 1.000`.

**Propagation test result (2026-09-22 03:00, W_prop plain vs weighted):** NULL. Both arms converge in 3 iterations
(greedy lp=8 moves poses 0.25 deg: cdrgn10028 poses are already converged), equal-weight maps from the final poses
are identical (eo-res 3.769 A both; FSC(plain,weighted) 0.985-0.997 at 4-8 A); the weighted reference resolves one
shell lower (3.99 vs 3.95, neff 95k vs 105k). The machinery works; the claim ("better poses than the equal-weight
reference") needs a test from WORSE starting poses (abinitio3D output or perturbed poses) with the mode-kernel column.
wplain (plain start + weights inside every reconstruction) reproduced plain's poses EXACTLY (eqrec FSC 1.000 at 4-8 A, identical per-particle stats): the propagated weighting changes nothing at converged poses; only the starting reference made the 0.985-0.997 plain-vs-weighted difference.

**Worse-poses ribosome test (2026-09-22 11:55, W_prop2, mode-kernel column):** from +-5 deg / +-2 px perturbed poses all three
arms (plain, weighted, wplain) converge in 3 iterations to identical poses (eqrec FSC 1.000 plain vs wplain; 0.990-0.999
plain vs weighted = one shell from the kernel-map start), eo-res 3.74-3.77 A. The weighted reconstruction inside
refine3D does not change the ribosome poses at any tested start. The membership is now ALWAYS the mode kernel (the
island/continuum switch was removed at the user's request: "a point on the whole embedding"). Fast harness: the toy2c
quick loop (`~/flexw_toy_chain.sh`, `flexw_toy_loop.sbatch` + `flexw_toy_arm.sbatch`, scorer `toy_pose_err.py` vs
toy2c/all.txt ground truth, ~30 min end to end).

**VERDICT 2026-09-22 (peaked toy toy2p, ground truth):** canonical weighting cannot improve population poses. Even the
EXACT volume of the peak state (11% of particles + neighbours) refines poses worse than the blurred consensus (5.66 vs
4.75 deg mean rotation error; shifts better 0.40 vs 0.81 px); every weighted reference sits on the plain floor; the
mode finder itself is correct (mode = true peak state 83). Uniform toy: exact central state 8.87 vs blur 6.96. The
blur is the best SINGLE reference; only per-particle references (mean + U z_i, the flex E-step projections) can beat it.
Harness: toy2p (`toy2p_sim_state.sbatch` array + `toy2p_finish.sbatch`, counts must be even), `flexw_toy_loop/arm/ladder`
with TOYDIR=toy2p, refine=greedy lp=10, scorer toy_pose_err.py; ~30 min end to end.

**FIRST POSE GAIN (2026-09-22 15:10, branch flex-pose-marginal, worktree /Users/atifao/SIMPLE_posemarg, batch2
~/SIMPLE_posemarg):** the mixture-marginal pose polish (`simple_flex_pca_em_polish.f90`, SIMPLE_COV_POSE_POLISH=1 after the
paired fit; MCFA K=16 marginal -2logL + data terms, 156 local trials) takes the peaked toy from 4.84 (start) / 4.75 (plain
refine3D floor) to 3.85 deg mean rotation error. Caveats: moves saturate at the 6-deg radius; from TRUE poses it drifts
2.9 deg (score bias); rigid-tangent deflation kills the gain (4.89). Launcher `~/flexw_toy_polish.sbatch`
(env TOYDIR NAME START KMIX DEG PASSES; SIMPLE_COV_POSE_POLISH_DEBUG/MAXP for diagnostics). Plan section 10.

**Polish round 4 (16:10):** centred search (each pass moves once to the best trial, step/3 per pass; the first
version compounded steps inside a pass) at lp 12, radius 3, 3 passes: 2.78 / 2.48 / 4.42 deg (view-dir 2.16) vs plain
floor 4.75 / 4.17 (3.89) = -41% / -45%; truth drift 2.18 at lp 12. Band is the lever (full band 3.75, lp15 2.98,
lp12 2.78; lp20 too coarse; rank 24 worse); smaller steps win; a second fit on hill-climb poses HURT (3.96) -> march
only from low-bias poses. Frequency-march rule: score up to the band where the basis is faithful (basis cross-half
FSC), refit, raise.

**Polish rounds 6-7:** lp 8 radius 3 and lp 10 radius 4.5 both 2.55 deg mean (plain floor 4.75, -46%); truth drift 1.97
at lp 10; the band ladder pays down to 8 A on the toy (Nyquist 4.9); radius must reach the perturbation (2 too small).
**Polish round 8:** lp 8, radius 4.5, 3 centred passes: 2.37 / 2.05 deg (view-dir 1.87) = half the plain floor; truth drift 1.61 at lp 8. Best recipe so far on toy2p.
**Polish lp 6, radius 4.5:** 2.11 / 1.75 deg (view-dir 1.69), -56% vs plain floor, 60% of particles within 2 deg. The band ladder is monotone to 6 A (finer = better); the early 'full band' arms were at flex's derived default band, not Nyquist.
**Toy band ladder closed:** lp 5 (Nyquist) radius 4.5: 1.95 / 1.53 deg (66% within 2 deg) vs floor 4.75; monotone 16->5 A. Ribosome polish arms (box_crop 128, lp 8) running 2026-09-23 as W_polish/ribo_conv + ribo_pert; the perturbed-pose fit's het-FSC is 14 A vs 8 A from converged poses.

**RIBOSOME GAIN (2026-09-23 05:40):** marginal polish on the converged cdrgn poses (box_crop 128, lp 8, radius 3): poses
move 0.42 deg mean, the equal-weight map improves 3.769 -> 3.521 A (FSC at 4 A 0.25 -> 0.34). The perturbed arm was void:
10028's W_prop2/start has 4 px gaussian shift errors (pre-fix oriops) and the polish has no shift search -> 10 A map;
angle-only start W_prop2/start_ang made, rerun + a second polish round (ribo_conv2) submitted. Launcher
`~/flexw_ribo_polish.sbatch` (NAME START LP DEG PASSES), ~7 h per arm.

**Shift search (2026-09-23):** the first form (all rotation passes, then one shift search) cannot fix rotations when 5 px shift errors are present (toy2p: 4.83 -> 4.91 deg, shifts 5.2 -> 2.5 px). The passes must interleave: shift grid first on a per-thread working plane, re-phase it in place when a delta is accepted, then rotations. Rotation search on a mis-shifted plane moves every particle by the full radius (it absorbs the shift error). refine3D greedy baseline from the same start: shifts 0.9 px but mean rotation 6.84 deg (tail flies off).
**Batch2 prototype trees:** SIMPLE_posemarg/posemarg2 are an older master than the local worktree (201265f41); a partial file overlay breaks the build (missing simple_flex_weights_state). Sync the whole worktree (rsync -rc, exclude .git/build), never a file subset. Worker entry points called from simple_flex_pca_model must be in the em parent's `public ::` list.
**Distributed polish verified (2026-09-23 11:50):** PCA_STAGE_POSE over 4 parts reproduces shared memory on toy2p (4.14 vs 4.13 deg). Three traps on the way: the worker stage validator in simple_flex_pca_strategy has a fixed case list (add new stages there); the master never initialises the MCFA mixture in distributed runs (E-step-side), so mixture state cannot be saved from the master; a per-thread fplane copy must be a full derived-type assignment (ctfsq_plane is read by the former). Interleaved shift search still stalls at ~half the shift error (5.2 -> 2.6 px); re-phasing convention verified against gen_fplane4rec, so suspect latent absorption of translations or the rotation/shift trade-off (mean-only diagnostic queued).
**Shift score must be mean-only (2026-09-23 12:30):** the marginal likelihood is flat in shift because the basis spans translations (toy2p shifts-only: marginal 2.55 px vs mean-only 0.56 px after 5 passes); the mean alone is too weak for rotations (4.50 vs ~3 deg). Hybrid: shifts by the mean-projection score, rotations by the marginal. Also: the grid search always accepts a move (min over 78 noisy trials < centre), so a biased score turns into full-radius drift.
**Rotation floor is systematic (2026-09-23 14:30):** from truth the marginal polish moves ~all particles by ~1.1 deg at lp 8/box 72 (toy2p), grid and quadratic-fit steps alike, shifts exact or not: the objective's minimum is displaced (latent partly explains rigid rotation of displaced states), not noise acceptance. Angle-only polish 4.84 -> 2.38 (lp 8) / 1.95 (5 A). Angles+shifts loss (4.1) was pass ordering: run the full shift ladder before each pass's rotation trials (SHPASSES). Env switches on the prototype: SIMPLE_COV_POSE_POLISH_{SHIFT,SHPASSES,SHIFTSCORE=mean|defl|marg,STEP=grid|quad,SIGTOL,ROTSCORE=marg|fixz,MEANONLY}.
**Polish floor = per-particle angular noise; posterior-mean step fixes most of it (2026-09-23 21:00):** truth-start drift tracks the search radius below ~1 deg at lp 8 (flat landscape), axes random (not a convention), state-map scoring drifts the same (not misspecification). STEP=soft (trial rotation vectors weighted by exp(-(s-smin)/2), move taken unconditionally) cuts the drift 1.11 -> 0.40 deg while keeping the angle-only gain (2.38). Split-half gate and quadratic fit did not help; fixz/state-map/deflation are not fixes. Ribosome round 2 with argmin degraded 3.52 -> 3.83 A; soft-step rerun pending. Launcher trap: flexw_ribo_polish.sbatch unsets ALL SIMPLE_COV_* (use flexw_ribo_polish3.sbatch).
**Trap (2026-09-24):** the prototype's shift search defaulted ON (3 px) on the posemarg3 tree; the ribosome soft round applied 0.42 px shifts unasked and came back 4.05 A (confounded). Default now off; always set SIMPLE_COV_POSE_POLISH_SHIFT explicitly and print the "shift search radius" line check.
**Ribosome verdict (2026-09-24 03:30):** round 2 of the polish degrades regardless of step rule (argmin 3.83, soft 3.80 vs round 1 3.52 A) because a box_crop-128 / lp-8 model cannot inform sub-0.2 deg moves; the fix is band (band-split score: marginal <= 8 A + consensus-map data term beyond; or a box-256 fit). Soft step is still the right step rule (toy floor 1.11 -> 0.40 deg).
**CORRECTION (2026-09-27):** the ribosome "3.769 -> 3.521 A" polish gain was a selection artifact: the ribosome launcher's post-run `selection states=1..8` keeps 51k of 105k particles (flex assigns 16 states) and the consensus poses on those same 51k give 3.471 A. The marginal polish (box_crop 128, lp 8) LOSES 0.05 A in round 1 and 0.3 A in round 2 on the ribosome, with every step rule and the band-split score. Only the toy gains (real pose errors). Any ribosome eqrec comparison must use identical particle sets; check state counts in the exported ptcl3D table first. oriops pose import: `oritab=... oritype=ptcl3D outfile=<existing project>` with NO projfile (with projfile the table is ignored).
**gen_fplane4rec needs memoize_ft_maps for ITS grid (2026-09-28):** the flex prep memoizes the cropped grid per batch; any full-box plane prep (band-split polish, pose_cv) must call memoize_ft_maps([boxpd,boxpd,1], smpd) first (outside OpenMP) or it reads wrong/unallocated maps. pose_cv (prg=pose_cv, held-out per-particle score) exists on the prototype branch; ~/pose_cv_compare.py compares two tables.
**Held-out verdict (2026-09-28):** on identical 51k particles, pose_cv (band 20-3.5 A) says: marginal-only polish rounds lose (-1.9, -2.4 vs consensus, 25-28% improved); the BAND-SPLIT polish (flex marginal <= 8 A + consensus-map data term 8-3.5 A, soft step, no shifts) gains +1.6 +- 0.02 over the consensus with 68% of particles improved, at identical FSC-0.143 (3.471). Use pose_cv + the clean-51k protocol (flexw_ribo_clean51k.sh) for every ribosome comparison; never FSC alone.
**Band-split final (2026-09-28):** from the consensus directly +1.76 +- 0.02 (74% improved); a third round adds nothing (+1.64, moves 0.11 deg): one round, consensus start. Gain sits at 9-4.4 A; FSC-0.143 unchanged (3.471) while shell FSCs rise ~0.005.
