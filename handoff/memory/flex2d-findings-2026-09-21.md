---
name: flex2d-findings-2026-09-21
description: what the flex2D tests established: low-rank family cannot align (40-57 deg), cluster templates retain poses to 0.9 deg, random-projection init is 1000x too small, class search must exclude the mean
metadata:
  type: project
---
- The pure low-rank family mu + U z is a weak in-plane aligner on the 10028 ribosome: retention spread
  from the reference poses 57 deg (rank 16), 52 deg (48), 40 deg (96) at 20 A; shifts fine (1-2 px).
  Rank-vs-views (n-width) limit, confirmed, not a convention bug (FFT sweep == direct evaluation).
- Cluster templates (free class images) as alignment references retain the reference poses to 0.9 deg /
  0.76 px (t12) and give class resolution 11.2 A vs 10.5 A for the reference at equal subset (8000, box 182).
- The random-projection basis init (random z) is ~1000x too small in power; PPCA EM inflates it too slowly.
  Fixed by prior-free (ALS) probe iterations (`n_probe_iters`, prior_lambda=0) before the PPCA E-step.
- Template validity must come from populations snapshotted at solve time (`pop_solved`), not the running
  accumulators, which the per-iteration reset zeroes first.
- With class search, the noise-free global mean must NOT be a candidate: a young single-particle template
  always loses on its own noise energy, so nobody gets a label. Choose among templates only (as cluster2D).
- Random-subset templates (~270 randomly rotated ptcls each) are symmetric blobs and never break symmetry;
  cluster2D's `cls_init=ptcl` (one particle per template) is the start to copy.
- The abinitio2D schedule (6 stages: halve to lpstop, maxits over 4 marching stages, nits_per_stage after,
  trs=0 in stage 1) is now the flex2D default; MSA = van Heel's multivariate statistical analysis.
Env switches: FLEX2D_DEBUG, FLEX2D_FIX_LABELS, FLEX2D_CLS_SEARCH (default on). See [[flex2d-branch-and-test-setup]].
- Embedding as a (class, rotation) SEARCH INDEX works: at reference poses/labels (240 templates) the
  index (K+1 correlations) puts the true class in the top-8 for 78% at its index-best rotation (random
  8-of-240 = 3%), and the index-guided search (shortlist 8 + current class, FLEX2D_SHORTLIST) retains
  the reference solution exactly (g12: 1.0 deg, 0.6 px, 10.49 A = the reference restoration).
- Embedding as a TEMPLATE PRIOR is negative: rank-32 predictions of class averages have FRC 0.15 to the
  reference averages at 20 A, so the shrinkage weight stays 0.99 and shrunk == raw (FLEX2D_TMPLTEST).
- Embedding as MEMBERSHIPS is negative: k-means labels give 15 A class averages vs 10.5 A for cluster2D's
  labels at the same poses/count (m16/m32).
- ALS is scale-invariant: rescale U by the latent's sample second moment after each prior-free solve
  (model%solve(l_rescale)), else PPCA shrinks z to zero afterwards.

**Autonomous session 2026-09-21 evening (ab initio dynamics), all on the 8000-particle subset:**
- Rotation EM is sound: reference labels + randomised angles converge back to the reference (1.0 deg, rr1).
- Random-start failures were dynamics, fixed one by one: hard-EM freeze (0 changes by it 30) -> extremal
  randomisation + power-sampled rotation/class (decaying over 20 its); class collapse -> stochastic
  neighbourhoods + re-seeding; mean must not be a candidate; single-particle seeds (cls_init=ptcl);
  Gaussian stage-1 templates; leave-one-out templates (FLEX2D_LOO) gave the best NMI.
- Best random start (c60: LOO, Wiener 0, shortlist 8, latent off): NMI 0.516 / 62.8 deg / 9.1 A, all 60
  populated, vs abinitio2D b60 0.60 / 59.1 deg / 9.1 A. Speed 6 s vs 22 s per iteration.
- OPEN DEFECT: shift recovery. Reference shifts are rms 21 px (picking off-centre); abinitio2D recovers
  19-21 px from zero, flex2D only 6-9 px even with cluster2D's 5x5 coarse grid (s7: fixed labels/rotations,
  zeroed shifts -> 5.8 px). Increments stay ~1 px/iteration; cause not found (suspects: centering/common
  offset handling, search-window vs template-smear width, minimize with new_fixed).
- Negative: Wiener constant in the template ratio (hurt NMI), template FRC filter (no effect), white
  sigma (no effect), latent term in the evidence (destabilises, and SIGBUS runs were a rebuild under a
  running job, not code). Never rebuild while sbatch jobs use the binary.
- Bookkeeping trap: with nstages=2, ALL maxits iterations go to stage 1 (trs=0); use maxits small +
  nits_per_stage large for shift tests.

**2026-09-22 (shift defect resolved, full-dataset runs):**
- Shift search is correct (known +-3/10 px perturbations undone); large-offset failure is physics: ring
  averages over random offsets carry a J0(2 pi k d/box) factor (slow decay, sign flips) -> smeared templates
  have HIGH even/odd FRC (0.9 at 16 A) yet no coherent shift gradient; self-included templates pin, LOO
  ones random-walk. abinitio2D itself does NOT recover shifts from zero in isolation (b30z: 10.6 px
  residual, rotations scrambled to 62 deg) -> the isolated test is the wrong yardstick.
- The stored shifts ARE physical: raw low-passed centroids correlate -0.81 with them (rms 11 px/axis).
- PARTICLE PRE-CENTRING (FLEX2D_PRECENTER, default on; calc_shiftcen_serial needs
  construct_thread_safe_tmp_imgs) took the random-start subset to parity: pc60b NMI 0.545/62.6deg/9.75A
  vs abinitio2D b60 0.60/59.1deg/9.1A; shifts rms 16.5-17.9 of 21 px. LOO shift refs in the marching stages.
- Full-dataset f300 (83.6k ptcls, 300 templates, K=32): 29 min for 40 iterations (abinitio2D ~90 s/it)
  but NMI 0.30: the index shortlist at iteration 2 used the unscaled init basis -> non-finite scores ->
  maxloc 0 -> OOB write -> no class assignments. Fixed (index only with prior_lambda>0, finite guard).
- NEVER rebuild+install while a job runs on build/bin/simple_exec (SIGBUS); use build without install +
  snapshot binaries in ~/SIMPLE_flex2D/snap/simple_exec.<commit>.

**2026-09-22 (b), full-dataset gap closing (a300 = abinitio2D 61.7 deg / 5.4 px / NMI 0.761 / purity 0.585 / 5.02 A, 74 min):**
- f300c (shortlist 8->32->all, 2nd shift round, ml_reg delivery): 71.0 deg / 10.97 px / 0.620 / 0.372 / 6.11 A, 41 min.
  Its shift residual equals the ZERO-shift level (b30z 10.57 px): flex2D's shifts carry no information. Log
  signature: mean |shift incr| plateaus at ~2.5 px for 30 iterations (a300 decays 3.4 -> 0.3); pt3/pt10 prove
  the search+frame conventions (10 px perturbation -> 2.5 px, increments -> 0.01) so it is the joint dynamics.
- Structural differences found vs abinitio2D: (1) cluster2D CENTRES every class average each iteration
  (prep2Dref: img shift + add_shift2class) and flex2D never did; (2) abinitio2D runs calc_group_sigmas after
  EVERY iteration while flex2D kept the raw-power bootstrap sigma2 for the whole run (weights the signal band
  down by its SNR); (3) abinitio2D assigns classes through the prob table (population-balanced frontier) from
  stage 3 on; (4) flex2D's stage-1 templates are Gaussian-crushed (0.03 at the band edge) + FRC-filtered.
- Implemented (commits 0de8556fd, 5b82301ea, 05f973a52, 645fb6bcf): template centring in the polar frame
  (FLEX2D_CENTER), per-iteration global sigma2 from the residuals (FLEX2D_SIGMA_UPDATE, gen_sigma_contrib
  normalisation: half-plane sum/(2 pftsz), by parity, persists across stages), frontier global assignment from
  a 32-candidate table (FLEX2D_ASSIGN, stages >= 3). IMAGE-BASED polar centring FAILS: the band starts at
  ring 2 which is past the first zero of the ribosome disc transform, so the "brightest point" runs to the
  search boundary (27 px rms); the members' MEAN reference-frame shift (cluster2D center_type=params) is the
  right gauge (1.2 px rms per iteration). sigma2 update alone did not stop the random walk (f300f).
- cavg_centres.py: f300c cavgs 7.1 orig px rms off-centre vs a300 4.8 vs reference 6.1 (estimator-limited).
- 2026-09-22 (c) MECHANISM: the reference shifts have within-class rms 20.9 px (global 21.1, class-gauge
  scatter 2.5); pre-centring takes the global residual to 13.1 px; every flex2D shift-search variant then
  makes it WORSE in stage 2 (13.1 -> 15.7 by it 14) with increments of 3-4 px and cos ~0 to the truth,
  while class churn stays 40-90 %/iteration (abinitio2D churns just as much: a300's in-plane dist is
  ~90 deg/iteration for 25 iterations, yet its shift increments decay 3.4 -> 0.3). Under churn a shift
  search relative to a sequence of different templates converges only if all templates share ONE gauge:
  cluster2D's seg centring anchors every class average to its segmentation centre; flex2D templates
  carried the arbitrary CoM offsets of their founders (+ the members'-mean gauge drifts with churn,
  f300d/e/g/h/i/j all ~15.3-15.9 px). Image-anchored polar centring (afa327c8d: polar->Cartesian
  gridding + calc_cavg_offset, sign calibrated) is the cluster2D-faithful fix (f300o).
- The LOO shift reference was built from the RUNNING accumulators (zeroed every iteration): a zero
  reference for the first batch, a partial mean afterwards (f1ae9fd1c fixes it with solve-time
  snapshots; LOO now opt-in). power_pick weights exp((ev-emax)/wsq) are ~uniform over the top nb.
- 2026-09-22 (d) RETRACTED: the in-run 'shift residual vs true' diagnostic captured xy_true AFTER the
  fresh start zeroed the shifts, so it reported the rms SPREAD of the run's own shifts, not the error
  (its '13.1 after pre-centring' = the pre-centring rms 13.13 itself; the zero-start runs' '3.4-3.8 px'
  = shifts that never moved). On the real scorer the zero-start runs (f300r/s/t/u) score 17.9-19.1 px
  (zero-shift level 20.9), the pre-centred ones 10-11 px: PRE-CENTRING IS A REAL GAIN, keep it. Fixed in
  the commander (capture before the zeroing). Trust only compare_poses on test_ptcl2D.txt.
- Real full-set standings (compare_poses, a300 = 61.7 / 5.4 / 0.761 / 0.585 / 5.02 A): f300f (sigma2
  update) 71.8 / 10.29 / 0.649 / 0.405 / 5.99 is the best flex2D; f300m/n (LOO off or corrected,
  new search order, assignment, no centring) 77-81 deg / 9.96 / 0.61; every centring variant ~15 px.
- 2026-09-22 (e) confirmation: pure defaults with cluster2D's shift-first/refine-before-choice order
  (f300x) 77.5 / 10.42 / 0.604 / 0.354 / 6.35; + assignment (f300y) 78.0 / 10.38 / 0.594. Both below f300f
  (71.8 / 10.29 / 0.649 / 0.405 / 5.99): per-candidate shift freedom lets wrong templates win. Search order
  made opt-in (FLEX2D_SHIFT_FIRST), defaults = f300f's configuration + LOO fix (commit after e46b69602);
  f300z is its confirmation. The corrected true-shift diagnostic: pre-centring 21.1 -> 13.2 px, then the
  shift search gains only ~0.05-0.1 px/iteration (12.8 -> 12.3 over 25 its) -- THAT is the open gap.
- f300z (defaults at 4b8cd83c2): 74.8 / 10.04 / 0.594 / 0.348 / 6.43. f300f (snapshot 5b82301ea) stays the
  best measured run; whether its buggy LOO reference or run-to-run variance explains the 0.05 NMI is open.
- 2026-09-22 (f) index A/B on the full ribosome set: WITHOUT the embedding shortlist (g300f, full class
  search) 78.5 / 9.72 / 0.631 / 0.422 / 5.64 A in 55 min vs f300z (shortlist 8) 74.8 / 10.04 / 0.594 /
  0.348 / 6.43 A in 42 min; joint rotation-shift searcher (g300j) 77.7 / 10.02 / 0.590: no shift effect.
  The index costs labels/resolution, not only candidates. 10076 (Davis 50S intermediates): the raw
  stack is protein-DARK, SIMPLE wants protein-bright (10028 core-ring +0.25) -> use
  L17Combine_weight_local_neg.mrc; abinitio2D crashes in simple_class_frcs (set_frc 64/128 vs 80) at box
  256 / 1.31 A even with autoscale=no, box 320 runs. Dev set: 10076_flex2D/sub20k (20k, seed 20260922).
- 2026-09-22 (g) 10076 20k dev subset (10076_flex2D/sub20k, box 256 from the _neg stack, 60 classes, mask
  260): abinitio2D (box 320) median 15.0 A / best 5.2 / 40 of 60 classes near-duplicates (deviation
  correlation > 0.9 at 20 A), ~20 featureless sink classes; flex2D soft EM (5f2409457, now default)
  10.8 / 5.9 / 8 duplicates in 8 min; robust w-estimator weights neutral (10.8/6.6); class basis rank 3
  neutral (11.0/6.1, modes = residual pose error); 90 iterations neutral; a 6-8 A polish (labels fixed,
  b3fed8763 makes shifts searchable without class search) is NEGATIVE (14.0/7.5): above 8 A the
  alignment against own averages degrades. Scripts: cavg_diversity.py, mrc_montage_ranked.py.
- 2026-09-22 (h) POSE TEST (10076 20k): same abinitio2D class members restored with abinitio2D poses 14.2 A
  median / 5.18 best vs with flex2D poses 20.5 / 7.76: the gap is POSE accuracy (rotations agree to 9 deg
  within joint cells, shifts differ 7 px rms). At fixed good labels, 10 flex2D in-plane iterations WITH
  the template FRC filter 13.1/5.74, WITHOUT it 9.33/5.24: the per-ring FRC weighting of the templates
  blunts the in-plane objective (cluster2D's ml_reg refs keep full amplitude). In the ab initio loop
  FLEX2D_TMPL_FILTER=0 gives uncapped eo-FRC median 7.45 vs 8.28 (soft EM) but 16 vs 8 near-duplicate
  classes; without the Gaussian too the run collapses to 30 classes. abinitio2D uncapped: 12.3 median,
  3.9 best, 40 duplicates. Restoration recipe for pose tests: import_particles stk+oritab (deftab lines +
  e3/x/y/class/state/eo, ONE source only), update_project qsys=local inside the job, make_cavgs nparts=1.
- 2026-09-22 (i) CODE STUDY RESULT: cluster2D's marching stages use gauref=yes (Gaussian at the stage lp) and
  ml_reg=no, i.e. NO FRC weighting of the alignment references until the final stages; ml_reg (tau=1) then
  damps by ~FRC/(FRC+1-FRC) above ring max(6,hp). flex2D had Gaussian in stage 1 only + FRC weight
  2f/(1+f) in every stage. Port (f7bb5f50b, FLEX2D_MARCH_FILTER=cluster2d default): 10076 20k uncapped
  eo-FRC median 7.62 / best 5.24 / 6 near-duplicates (soft EM 8.28 / 8 dup; unfiltered 7.45 / 16 dup);
  ribosome true residual 11.56 px at it 22 vs 12.3 (first config moving the shifts). Units/frames between
  the programs verified identical (pftc shifts in original px, prepimg4align -x*crop_factor).
- 2026-09-22 (j) RIBOSOME with the cluster2D filter schedule (h300m, f7bb5f50b defaults): 55.0 deg / 9.28 px /
  NMI 0.737 / purity 0.536 / 5.88 A (q75 12.4), 47 min -- beats abinitio2D on in-plane (61.7), near on labels
  (0.761/0.585), shifts still 9.3 vs 5.4, tail of weak classes. 10076 full set 200 cls: abinitio2D 9.02 A /
  best 3.88 / q75 16.8 in 121 min; flex2D soft EM (pre-schedule) 9.86 / 5.24 / 12.7 in 84 min; f200m =
  full set with the new defaults.
- 2026-09-22 (k) ZERO-START with the new defaults (h300m0): 21.1 -> 19.75 px over 10 shift iterations,
  abinitio2D from the same start reaches 5.4: offset sub-clustering is NOT the mechanism flex2D lacks. The
  de-rotate/de-shift conventions of the M-step data (conjg(shmat), inverse rotation index) were re-read
  and are consistent with the reference sweep. Per-iteration shift refinement in flex2D contributes ~0.1
  px/it vs ~1 px/it in abinitio2D; cause NOT found after testing every structural difference. Pragmatic
  route queued: flex2D partition -> abinitio2D cls_init=prev for the shift polish (a200p on 10076 full).
  10076 full uncapped: abinitio2D 7.29 A median / 3.33 best / 145 of 200 near-duplicates; flex2D (soft EM)
  8.06 / 5.24 / 76 duplicates.
- 2026-09-22 (l) ISOLATION on the 8k subset (reference labels, corrected residual): +-20 px start, rotations
  right: flex2D 16.4 -> 13.0 and stalls (steps 2.1 -> 0.4), scorer 12.5 px / 23 deg; abinitio2D from the
  SAME start (cls_init=prev, 17 its) 9.8 px / 60 deg. +-10 px + random rotations: rotations freeze wrong
  (62 deg), shifts 8.2 -> 7.4. +-10 px + class search on: recovers (3.7 px). So the smeared-template fixed
  point is GENERIC (both programs stall from a smeared start); abinitio2D's 5.4 px comes from its
  zero-shift/random-group dynamics never building smeared templates, which flex2D's stage 1 does not
  reproduce (h300m0: 21 -> 17.9). cluster2D's bpgau2d = exp(-3.56 (k/klp)^2): the ported Gaussian is
  exact. Pre-centring (13 px centroid noise) is a smeared start for either program.
- 2026-09-22 (m) STAGE-1 SORTING TEST (8k, zero-shift random start, 60 cls): within-class rms of the true
  offsets after stage 1 = 20.8 px for BOTH abinitio2D and flex2D (global 21.3): neither sorts by offset.
  abinitio2D's full-run shift recovery (21 -> 5.4) happens in stages 2+ from smeared averages, and is NOT
  reproduced by its own cls_init=prev schedule from a +-20 px start (9.8 px). Mechanism unknown; every
  code-level difference tested. 10076 full set with the new defaults (f200m): capped median 8.14 / best
  5.24 / q75 9.75, 86 min (a200n 9.02 / 3.88 / 16.8 in 121 min).

(n) 2026-09-23, THE shift gap resolved: cluster2D does not recover 20 px offsets by shift search. On the 8k
ribosome subset from zero shifts, abinitio2D with center=no stalls at 20 px exactly like flex2D; with the
ptcl3D field zeroed it recovers to 9 px as usual. The recovery is the per-iteration class-average centring
(center_type=mass: calc_shiftcen_serial at cenlp on the low-passed masked cavg, applied to every member
rotated into its frame via add_shift2class/map3dshift22d; off in stage 1, which_iter > 2, uncapped). Stage-1
classification at zero shift rotates particles so their off-centre blobs coincide, so one class offset
corrects each member (first stage-2 iteration: 11 px jump, cos 0.82 to truth; the searcher adds ~3 px/iter).
The fixed-rotation shift landscape is flat in BOTH programs (grid argmax ~1-2 px, cos 0.2; smeared refs have
no power beyond 40 A and ring 2 = mask self-overlap). flex2D port: FLEX2D_CENTER=mass (default from
2aaf8f504) -> 8k fresh 21.2 -> 10.5 px (abinitio2D 9.1). Not the cause: extremal/sampling, soft EM, robust,
sigma (also flat sigma), Wiener eps, embedding (neigs=1 stalls too), searcher class/limits/frames,
polar-vs-Cartesian refs (cluster2D's own cavgs loaded into flex2D's slots still stall). Diagnostics kept:
FLEX2D_SHDBG (grid scan + per-ring dump), FLEX2D_CARTREFS/CARTLP, FLEX2D_SEED3D, CLUSTER2D_SHDBG
(snhc_smpl pose print). Trajectory harness: 10028/flex2D_tests/traj (var*.sh, per-iteration true residual).

(o) 2026-09-23 later: BOX CAP. Every flex2D "best class" number was the Nyquist of the cropped box its final
make_cavgs ran at (5.24 A = 160 px/2.62 A on 10076, 5.30 A = 182 px/2.65 A on 10028); abinitio2D restores
its LAST iteration at the full box (320 px). Fixed from fd175a51e (make_cavgs at box/smpd, autoscale=no);
all earlier flex2D resolution comparisons understate flex2D. ML-reg-from-stage-2 schedule
(FLEX2D_MARCH_FILTER=frc): ribosome median res 5.48 vs 6.07 but poses worse (9.65 px, 67 deg, NMI 0.66);
10076 worse everywhere (8.06 A, 107 duplicates) -> Gaussian marching stays default. FLEX2D_SH3 (cluster2D's
final alternating rot/shift search) neutral on 8k (11.1 vs 11.0); FLEX2D_SHIFT_FIRST worse (12.2).
DISTRIBUTED flex2D (nparts>1, from d2859a1d4..0a4f24466): master writes flex2D_model_master.bin +
flex2D_sigma2_master.bin + project, schedules `flex2D part= which_iter=` workers via qsys (workers replay
the stage schedule, run one E-step over their range, write flex2D_acc_partN.bin, flex2D_stats_partN.bin
(incl. re-seed pool), algndoc_N.simple, JOB_FINISHED sentinel); master merge_algndocs + read_acc_add.
Gotchas fixed: qsys resolves the worker binary as SIMPLE_PATH/bin/<basename> (snapshots need a bin/
symlink and SIMPLE_PATH exported), algndoc must be written with isegment=PTCL2D_SEG, workers must call
qsys_job_finished, model%nacc is nacc(NHALF) in the reduce, the project compenv gets overwritten by the
run's cline so qsys_name must be put on the cline for the nested make_cavgs. Nested sbatch works on batch2.
8k/4 parts: identical trace to shared memory, 401 s for 40 iterations (overhead-bound at this size).
(p) distributed snapshots need BOTH binaries under snap/<H>/bin: simple_exec (flex2D workers) and
simple_private_exec (make_cavgs parts); export SIMPLE_PATH=snap/<H>, SIMPLE_QSYS=slurm,
SIMPLE_QSYS_PARTITION=csbdevel in the master job and run update_project qsys_name=slurm first.
Full-set timing with 8 parts x 24 threads: ribosome 300 cls 40 iterations in ~17 min (19 s/round),
10076 200 cls in ~20 min (25 s/round) vs 47/85 min on one node. Full-box finals (fd175a51e):
ribosome median 5.99 A best 3.86 q 4.51/10.05 (poses unchanged 8.98 px); 10076 median 7.76 q 6.35/8.56
("best 2.62" = full-box Nyquist artefact of the eo_res best-class statistic; use median/quartiles).
(q) 2026-09-23 final distributed full-set numbers (ca9c91efb, 8 parts x 24 thr, mass centring, Gaussian
marching, full-box finals): ribosome h300x 40 iterations in ~17 min: 9.01 px / 61.7 deg / NMI 0.713 /
median 5.19 A q 4.55-7.20 (abinitio2D a300: 5.40 px / 61.7 / 0.761 / 5.02, 74 min single node);
10076 f200x ~20 min: median 7.76 A q 6.35-8.22, 72 near-duplicates (abinitio2D a200n 7.29, 5.11-13.97,
145 dup, 121 min). Distributed make_cavgs parts fail on 10028 ("sigma2 values must be finite and positive"
in resample_sigma2, inactive particles' zero sigma rows); single-process restoration works -> restore
without nparts for projects with state=0 particles until fixed. cavg_diversity.py at 320 px takes ~15 min.
(r) 2026-09-23 shift remedies on 8k (end residual, default 11.0, abinitio2D 9.1): band-limited shift
search (FLEX2D_HPSH 100/150, +LPSH 30) WORSE 12.8-13.7; bilinear polar->Cartesian gridding for the
centring (FLEX2D_T2I=bilinear) worse 13.1 (angle convention unverified); FLEX2D_SH3 on top worse.
What helps is a coarser/longer stage 1 (cleaner blob rotations -> bigger first centring gain):
lpstart=25/30/35/40 -> 10.6/10.1/9.9/10.1; FLEX2D_NITS1=16 (stage-1 iterations 8->16) 10.2;
lpstart=30+NITS1=16 9.5, lpstart=35+NITS1=12 9.6. Full-set distributed runs h300y (30/16), h300z
(35/16), f200y (10076, 30/16) launched 2026-09-23 with snapshot 1609ed340.
(s) full-set stage-1 results (distributed, 20 min each): ribosome lpstart=30+NITS1=16 (h300y) 8.16 px /
57.0 deg / NMI 0.733, lpstart=35+16 (h300z) 8.15 / 59.7 / 0.731 (default h300x 9.01/61.7/0.713;
abinitio2D 5.4/61.7/0.761). Class resolution: plain median over classes WORSE (6.65/7.42 vs 5.19) because
~120 classes end with <100 particles (default 31; abinitio2D 32), but the PARTICLE-weighted median is
better (4.55/4.47 vs 4.68; abinitio2D 4.47) -> the poses/classes of the median particle now match
abinitio2D; the remaining issue is class-population balance (no split/rebalance in flex2D). On 10076
(f200y) the coarse/long stage 1 is neutral-to-worse (7.76 median, q 6.99/9.32 vs 6.35/8.22).
Use particle-weighted resolution (from cls2D pop/res) as the class-quality metric from now on.
(t) REBALANCING (FLEX2D_REBALANCE=0.35 of mean pop, FLEX2D_REBALANCE_STAGES=n; 7caeb9ba5..adb14a2e2):
classes below the floor take a random half (hash of iptcl, recipient, iter) of the largest class + its
template. Ribosome full set, coarse stage 1 (lpstart=30, NITS1=16) + rebalance at stages 1-4 (h300s):
class-median res 5.02 A q 4.47/6.35 worst 14.2, 1 class <100 ptcls, 8.13 px, 54.0 deg, NMI 0.735,
purity 0.548 -- matches abinitio2D's class median (5.02) with better rotations; particle-weighted median
4.87 (abinitio2D 4.47, h300y 4.55: rebalancing spreads members out of the sharpest classes).
10076 full set, default stage 1 + rebalance stages 1-4 (f200s): median 7.11 A q 5.05/7.76 (abinitio2D
7.29, 5.11/13.97) but 120 near-duplicates (72 without; abinitio2D 145) -- split copies stay similar.
Stage-1-only rebalance (h300r/f200r): 80 tiny classes / 10076 7.62 q 5.90-8.06, 85 dup.
(u) EXTERNAL TRUTH for 10028 shifts = the project's ptcl3D field (RELION-refined). Within-class rms of
(2D shift - 3D shift): abinitio2D a300 5.51 px, flex2D h300x 10.33, h300s 9.75. So flex2D's shift deficit is
real, not an artefact of scoring against abinitio2D's answer. Test plan agreed 2026-09-23: (1) poses vs
ptcl3D, (2) 10076 class purity vs deposited state labels (cryoDRGN EMPIAR repo), (3) abinitio3D from each
program's class averages vs the deposited map. Mean+basis templates ("pred") were measured earlier:
FRC 0.15 to the reference averages at 20 A (rank 32) -> the linear subspace cannot represent the views.
(v) two-stage rebalancing (REBALANCE_STAGES=2): ribosome h300t 8.15 px / 59.1 deg / NMI 0.737 / class
median 5.08 q 4.43/6.52 / 23 tiny classes / particle-weighted 4.78 (4-stage h300s: 5.02, 1 tiny, 4.87);
10076 f200t median 7.35 q 5.44/7.76 (4-stage f200s 7.11 / 5.05-7.76 / 120 dup). JUNK: on 10076, particles
abinitio2D parks in >15 A classes are separated by flex2D's per-particle evidence 'corr' (AUC 0.78, worst
10% = 63% disc-class) better than by the latent norm (AUC 0.70; junk has SMALLER |z|, near the mean);
combining adds nothing. cryoDRGN EMPIAR repo has no 10076 state labels (poses/ctf/ind only).
(w) 2026-09-23 evening: INDEX SPEED: 8k/300 classes, FLEX2D_SHORTLIST=0 vs 8: 600 vs 583 s -> the class
sweep is NOT the cost; the index buys ~3%. GENERATIVE-MODEL ALIGNMENT (SOLVAR-style): latent term in
the evidence (LATENT_EV=1) 11.9 px vs 9.4 base; + shift search against template+U*zbar (FLEX2D_GENREF)
15.9 -> the basis absorbs misalignment, removes the pose gradient. AB INITIO 3D FROM CAVGS (10076):
abinitio2D averages 6.08 A; flex2D averages 52 A with the FRC-derived lp ladder (starts at 20 A instead
of 70) but 6.45 A (interim) with lpstart=70 forced -> flex2D's stored class FRCs are too optimistic for
abinitio3D_cavgs; must fix the frcs contract. Ribosome pair: both 43.9 A (protocol issue, uninformative).
POLISH abinitio2D from flex2D's h300s partition: 6.18 px / 53.7 deg / NMI 0.752 / purity 0.582 / median
4.87 A -- best on every count. STRUCTURAL SPLIT (3332df100): ribosome 7.95 px / NMI 0.743 / median 5.30;
10076 7.49 A q 5.44/7.91, duplicates 96 (from 120). SYNTHETIC SUITE: scripts/flex2D/synth
(make_truth.py, run_synth.sh, score_synth.py), cluster dir /mnt/beegfs/elmlund/afan/synth2D; maps in
10076_flex2D/emdb (emd_8440 50S 320@1.31, emd_2660 80S 360@1.34, 8441-8449 fetching).
(x) SYNTHETIC SUITE first results (50S emd_8440, 20k, box 160/2.62 A, 60 cls, SNR 0.05, +-10 px; scores vs
truth: rot spread / shift rms / view spread / keep): baseline abinitio2D 13.9 deg/3.43 px/17.9 deg/0.66 vs
flex2D 7.9/3.57/10.1/1.00; cone70 8.8/2.89/12.4/0.71 vs 4.0/3.03/5.2/0.93; junk20 14.8/3.62/19.6/0.60 vs
8.4/3.65/11.2/1.00 (noise junk AUC ~1 both). Shifts EQUAL on synthetic data; flex2D wins rotations/view
purity/retention. Both ~7-8 min on 4 parts. DEFLATE_POSE (04b3dec9a): lev 11.9->11.4, gen 15.9->14.5 (8k);
coupled alignment from the REFERENCE start (prev_gen_dfl) drifts to 6.9 px -> the coupled alignment
harms even a good start; keep the embedding decoupled from the pose search.
(y) SUITE round 2 (abinitio2D vs flex2D; rot/shift/view/keep, + state majority & NMI): rare view 1%:
12.6/3.31/15.6/0.73 vs 8.9/3.51/11.2/0.98; 200 cls: 5.0/3.10/6.3/1.0 vs 4.5/3.52/5.5/1.0; SHIFT +-30 px:
abinitio2D FAILS 76.7 deg/29.3 px/57 deg/0.11 vs flex2D 8.1/3.64/10.2/1.0 (pre-centring + centring);
SNR 0.02: 13.6/3.50/16.8/0.67 vs 8.3/3.58/10.4/1.0; MIXED STATES (emd_8440 vs 8445 50/50): abinitio2D
majority 0.547 NMI 0.004 (= chance) vs flex2D 0.710 NMI 0.126; 90/10: 0.900/0.005 vs 0.922/0.053;
5 states 40/30/15/10/5: 0.410/0.007 vs 0.541/0.120. abinitio2D does not separate the 10076 states in
2D at all; flex2D partially. DRIFT from the reference start happens entirely in the first shift iteration
(5.44 -> 7.81 px, rotations/labels unchanged) -> the coarse-grid shift refinement adds ~2.4 px random
error on aligned particles; FLEX2D_COARSE=0 (4c44d76c5) under test. Ablations prev_nocen/nosoft/norob/
nosig/notrs launched from the reference start.
(z) DRIFT CAUSE = CENTRING: from the reference start (5.43 px vs 3D truth) flex2D ends at 8.09; with
FLEX2D_CENTER=0 5.51 (no drift); no shift search 5.44; nocoarse/nosoft/norob/nosig 7.9-8.1. The mass
centring on the gridded polar template is BIASED by several px per class (a converged template measures
~5 px off) and label churn turns per-class offsets into per-particle error. Tests from the reference:
CENTER=seg, T2I=bilinear, CENMAX=3/1 (61e605d87). FROZEN-POSE EMBEDDING (global mean+16 basis at
abinitio2D poses, 50/50 synthetic states): supervised AUC 0.55, 2-means NMI 0.00 -> the global embedding
carries no state signal (views dominate); flex2D's classification separated states via templates
(NMI 0.126). Per-class basis (FLEX2D_CLS_RANK=2) at frozen poses under test for within-class state signal.
(aa) 2026-09-24 centring estimators from the reference start (5.43 px vs 3D truth): mass 8.09, seg 13.5,
bilinear 10.2, cap 3/1 px 8.1, PHASE (real-symmetry) diverges 25.7, ENERGY (squared image centroid)
diverges 26.3, RAW (members' raw-centroid mean, 4df8fd872) 9.58 -- first pass applies 8.8 px offsets on
converged classes. Mechanism: any per-class gauge offset x label churn = per-particle scramble; only an
estimator consistent with the run's own alignment centre is safe (flex2D's mass centring IS
self-consistent on its own templates: offsets -> 0 late), so the remaining 8 vs 5.4 px on the ribosome
comes from the early high-churn iterations while offsets are large. abinitio2D polish from the flex2D
partition (6.18 px) is the practical last-2-px recipe. Keep mass as default; phase/energy/raw are dead.
(ab) 2026-09-24 centring POLICY (b9d493b34, FLEX2D_CENTER_STAGES / FLEX2D_CENTER_TOL) does not help: from the
reference tol1/st2 7.93/7.94 px (default 8.09, off 5.51); from scratch 11.6 (full centring 9.4). The first
centring pass does all the damage: it applies 15 A rms on classes that cluster2D's exact estimator (lp 30 hard,
mask, clip<0, masscen; replicated in scripts cavg_centres_exact.py) puts at 7 A rms / 2.7 A median on the
Cartesian averages -> the gridded polar template biases the estimate ~2x. Untested fix: centre on restored
Cartesian averages (cluster2D's path). abinitio2D's own final classes are NOT at its estimator's fixed point
either (6.9 A rms), and a2d_nocenter drifts 4 -> 15 A and stays (the 20 px stall).
(ac) 2026-09-24 STATE SPLITTING on the 50/50 synthetic mix (t4_mix5050, EMD-8440 vs 8445): the pure classes come
from the structural split in rebalancing (28 rebalanced -> 22/60 classes >80% pure; abinitio2D 0/60, purity
0.549 = chance). ncls=120 + FLEX2D_REBALANCE=0.8 + REBALANCE_STAGES=4: 74/120 >80% pure, mean purity 0.825,
majority 0.807. FLEX2D_SPLIT=n (41e5d2684): n most bimodal classes (d' along dominant latent) split into the
n smallest; scripts/state_montage.py (cluster synth2D/scripts) draws purity-annotated montages.
(ad) 2026-09-24 split-all (FLEX2D_SPLIT, d'-ranked 34/120 and population-ranked 45/120 >80% pure) is WORSE than
plain rebalancing (74/120): forced splits of pure large classes = duplicates + churn. Supervised ceiling
(scripts ceiling.sh + truth_labels.py: truth partition view x state = 240 cls, cls_init=prev, REBALANCE=0):
poses free -> purity decays to 0.879 (170/240 >80%, majority 0.864); 25% of labels flip in the FIRST E-step at
lp 30 A with 43% rotation changes, i.e. the coarse stage scrambles states. FLEX2D_FIX_POSES also freezes
labels (mode 1 = embed only), so a poses-fixed/labels-free ceiling needs a pose-locked class search.
Bug fixed 2cc2033ce: cls_init=prev class count was per worker partition (accumulator mismatch).
(ae) 2026-09-24 RESOLUTION IS THE STATE LEVER: on t4_mix5050 (120 cls, REBALANCE=0.8 x4) lpstop=8/20/25/30 ->
majority 0.807/0.882/0.880/0.934 (104/120 classes >80% pure at 30 A; the 30 A ceiling is 0.922). Finer
lpstart (15/10) on the ceiling run made it WORSE (0.67/0.65): the subunit difference is low-res, high-res
noise dominates the fit. Full-band pose polish with labels fixed (scripts/polish.sh) keeps purity (0.881 from
the stop-20 run) and buys almost nothing in poses. Recipe: labels at 20-30 A, averages at full box.
(af) 2026-09-24 generalisation of the 30 A label band: 90/10 (8440 vs 8445) minority recall abinitio2D 0.00,
flex2D 120 cls 8 A 0.20, 30 A 0.50, 200 cls 30 A 0.65 (36 pure minority classes). 5-state 40/30/15/10/5
(8440..8444): only 8440 separates (300 cls 30 A recall 0.96); 8441-8444 (consecutive intermediates) never,
even with 300 cls and 6 rebalance rounds -> E-step discriminability, ceiling runs at 30/40 A pending.
scripts/state_recall.py gives per-state pure-class count + recall (use it instead of majority for
unbalanced mixes; NOTE it was once run on the login node by mistake - keep it inside jobs).
(ag) 2026-09-24 5-state CEILING (corrected truth_labels.py: class = nstates*(view-1)+state; the 2-state formula
collided labels for >2 states): from the perfect 600-class partition, held at 30 A or 40 A, 8441-8444 collapse
in the first E-step (recall 0.00-0.06) while 8440 holds (0.94-0.99). So flex2D's unsupervised result on t4_mix5
(8440 only) EQUALS the ceiling: consecutive 10076 intermediates are below the per-particle discrimination
threshold at SNR 0.05 for this E-step; more classes/rounds cannot help. Separable pairs: 8440 vs 8445 (50/50
0.934, 90/10 minority recall 0.65 at 200 cls).
(ah) 2026-09-24 CLASS ENRICHMENT -> 3D STATES (scripts/ab3d_state.sh, abinitio3D nstates=2 particle-level,
class-restrained sampling) on t4_mix5050 at iteration ~32: from the flex2D 30 A partition per-particle 3D state
purity 0.914 (9370/8913 correct, balanced pops); from flex2D 8 A 0.658; from abinitio2D 0.501 (one state
collapsed to 1.4k ptcls). Real 10076 two-phase (30 A labels + FIX_LABELS full-band polish): 200 cls median
7.23 A (vs 7.11 at 8 A, 9.64 unpolished), NMI vs published labels 0.085 (8 A: 0.040, abinitio2D 0.031);
400 cls 8.06 A, NMI 0.104. Published labels: 10076_flex2D/labels/empiar10076/analysis/major_classes/
published_labels_major.pkl (131899; -1 junk 26575; 0=A 1=B 2=C 3=D 4=E 5=F), scorer labels/real_recall.py.
The deposited-map projection-bank scorer (mapscore.sh) is NOT discriminative (view dominates) - drop it.
Random-split control: 96/120 vs 104/120 with the global-latent split -> the global latent adds little;
class-conditional latents (FLEX2D_CLS_RANK, FLEX2D_SPLIT_CLS, 90fcf18c3) under test with zc_auc.py.
(ai) 2026-09-24 CLASS-CONDITIONAL BASIS (FLEX2D_CLS_RANK) at 30 A, 120 cls, reb 0.8 x4 on t4_mix5050: rank 1/2/4
-> majority 0.949/0.948/0.929 (107/110/106 of 120 pure) vs 0.934 without; and the split direction NOW matters:
random split with the class basis 0.856 (87) vs latent split 0.948 (110) [without class basis: 0.904 vs 0.934].
So the embedding earns ~9 points once the class basis absorbs within-class variability. zc keys are dropped by
the fixed ptcl2D record (N_PTCL_ORIPARAMS, ori2prec) -> class latents travel in the worker stats block (4c8761fd2)
and are dumped to flex2D_zcls.txt (scorer scripts/zc_auc.py).
(aj) 2026-09-24 abinitio2D CONTROLS at 120 cls on t4_mix5050: default 0.552, lpstop=30 0.557 (1-2 pure classes,
recall 0.01-0.03) -> abinitio2D does NOT separate states even at the coarse band; the state separation is a
flex2D property (rebalancing structural split + class-conditional model), not a resolution-schedule artefact.
3D jobs cancelled by the user at iteration ~32 with the interim evidence (0.914/0.658/0.501) as final.
(ak) 2026-09-24 class-latent STATE SIGNAL measured (zc_auc.py on flex2D_zcls.txt, rank 2, 30 A, 120 cls): in the
21 classes still mixed at the end, AUC of zc1 0.71-0.74, best-of-2 0.82, 12-14 of 21 classes > 0.8 -> the
class basis carries a per-particle state score in the residual mixtures. Class split vs global split with the
class basis: 0.940/107 vs 0.924/103 (same run pair), i.e. equal within noise (run-to-run +-0.02).
Next lever: second-pass split of the remaining mixed classes by class-latent bimodality (d' on zc, min pop).
(al) 2026-09-24 real 10076 with the class basis (f200x: 200 cls, 30 A, CLS_RANK=2): NMI 0.088 vs 0.085 without,
recall C/D/E 0.05/0.06/0.14 -> no gain on real data. Class-latent modality on synthetic 50/50: BIC-GMM picks 3
modes for pure AND mixed classes (heavy-tailed per-particle noise), i.e. the states are separable in aggregate
(AUC 0.82) but not as visible modes (d' ~1.3 < 2); n-modality is not directly usable for choosing n.
(am) 2026-09-24 class-latent n-MODALITY: Ashman D of a 2-comp fit is ~3.5 in EVERY class (pure and mixed), not
even/odd (AUC 0.54), not pose/defocus (|r| <= 0.08), weakly evidence (0.2): the per-particle class-latent
noise is platykurtic and dominates; state AUC 0.68 over all classes / 0.82 in mixed ones. So the latent is
state-informative in aggregate but NOT n-modal per particle at SNR 0.05 -> mode counting cannot pick n;
denoise first (neighbour pooling / higher SNR) or use the aggregate score. Scorers: zc_modes.py, zc_eo.py,
zc_what.py in synth2D/scripts.
(an) 2026-09-24 cross-group 5-state (t4_mix5x: 8440/8441/8445/8446/8449 at 40/30/15/10/5), 300 cls 30 A: B and C1
separate (recall 0.84-0.91), D1/D2/E never (0 pure classes) with or without the rank-4 class basis. BUT in the 72
mixed classes the class latents carry the D1/D2/E signal: best-of-4 AUC 0.889, 54/72 classes > 0.8 -> the
representation knows the state where the E-step/split cannot resolve it; an unsupervised readout (2-means on
the latents) has not yet beaten noise (BIC modes = 5 for pure classes). SNR ladder test launched (t4_mix5050
at SNR 0.2) to see where per-particle modes appear.
(ao) 2026-09-24 WHY the class latent is not n-modal: its dominant direction is the WITHIN-CLASS VIEW SPREAD
(R^2 on the true tangent-plane view offset 0.90 at SNR 0.05, 0.87 at SNR 0.2, 0.76 on mix5x rank 4; ~10 deg
view spread per class), state R^2 0.04. The flat "bimodality" is the uniform view distribution. SNR x4 does
not change it (107/120 pure, AUC 0.90 in mixed classes, modes still 4-5 for pure classes). Deflating the top-1/2
PCA directions of the class latents lifts the residual's state AUC 0.71 -> 0.76 (mix5x) but still no modes.
Scorers: zc_view.py, zc_deflate.py. Real-data ceiling (labels/real_ceiling.sh) + oracle junk removal
(labels/real_nojunk.sh) running/queued.
(ap) 2026-09-24 rank-4 class basis + deflating the top-2 class-latent directions: unsupervised state AUC of the
residual's dominant direction 0.839 in the 22 mixed classes (13/22 > 0.8) vs 0.70 raw; top-3 deflation
over-shoots (0.695). Modes still never appear (Ashman D ~3.2 in pure and mixed). FLEX2D_SPLIT_DEFLATE=n
implements the deflated split direction (needs SPLIT_CLS=1 and CLS_RANK > n).
(aq) 2026-09-24 FLEX2D_SPLIT_DEFLATE=2 with rank 4 + class split at 30 A: 50/50 110/120, 0.954 (= best, noise);
cross-group D1/D2/E still 0 pure classes (recall B 0.84, C1 0.82). The split direction is not the bottleneck
for the hard states; the E-step erodes them (ceiling). Real 10076 ceiling first attempt was VOID (the f200v
ptcl2D.txt dump had been doubled with class=0 lines -> 6 classes imported); scripts now dump the view run
fresh (view_ptcl2D.txt). Oracle junk-free run + ceiling rerun in flight.
(ar) 2026-09-24 real 10076 ORACLE JUNK REMOVAL (labels/real_nojunk.sh; import_particles IGNORES a state column ->
use prg=selection oritype=ptcl2D infile=<0/1 per particle>): 105,324 labelled particles, 200 cls, 30 A ->
NMI 0.094 vs 0.085 with junk, recall C/D/E 0.07/0.15/0.15 -> junk is a minor loss. Two real ceiling attempts were
VOID (writer used the last parsed class for every particle -> 6 classes); third attempt running with the fix.
(as) 2026-09-24 REAL 10076 CEILING (valid, rceil200: published label x 200 views = 1143 populated classes, poses
from f200v, relabel at 30 A, rebalance off): purity 1.0 -> 0.674, NMI 0.193; recall A 0.34, B 0.13, C 0.36,
D 0.48, E 0.54, F 0.28; churn 1-7%/iter (gradual erosion, not a first-step collapse). Unsupervised best is
NMI 0.104 (400 cls) -> the search reaches ~half of what the E-step can hold, and the E-step itself holds only
about half the labels. Both levers matter on real data. Chain launched: ceiling at 40 A, evidence-pruned 300 cls
(labels/real_prune.sh), ceiling at 20 A.
(at) 2026-09-27 real 10076 band ceilings: 40 A NMI 0.151 (worse), 30 A 0.193, 20 A 0.192 -> 20-30 A is the band.
Practical evidence pruning (labels/real_prune.sh: drop bottom 20% by flex2D corr; 37% of the dropped are
unlabelled vs 20% base) + 300 cls at 30 A: NMI 0.111, purity 0.547, recall C/D/E 0.10/0.20/0.20 = best
unsupervised real result (400 cls unpruned 0.104; oracle junk-free 200 cls 0.094). Gap to the ceiling ~2x.
Next chained: pruned 400 cls at 30 A, pruned 300 cls at 20 A.
(au) 2026-09-27 real 3D test queued: labels/real_ab3d.sh (abinitio3D nstates=4, 10x8) from rprune300 (flex2D
pruned 300 cls 30 A) and from a200n (abinitio2D), scored per-particle vs published labels; chained after the
prune ladder (rprune400, rprune300_lp20). Interim scoring recipe: copy the project, print ptcl3D, confusion.
(av) 2026-09-27 real ladder: pruned 400 cls 30 A -> NMI 0.125, purity 0.584, recall C/D/E 0.17/0.24/0.28 (new
best; 300: 0.111); pruned 300 at 20 A 0.080 (worse; 30 A stays). Next: pruned(0.3) 600 cls, then labels-frozen
polish of rprune400 (rprune400p) for the montage. ab3d_prune300 (4 states) running, ab3d_a200n queued.
(aw) 2026-09-28 real ladder continues: pruned(0.3) 600 cls 30 A -> NMI 0.148, purity 0.632, recall C/D/E
0.20/0.27/0.39 (ceiling 0.193). 300/400/600 = 0.111/0.125/0.148: class budget is the real-data lever; 800
queued. rprune400p (labels-frozen polish) median 8.22 A. NOTE ab3d_a200n dir predates this campaign (old
abinitio3D_cavgs log inside) - the new job overwrites it.
(ax) 2026-09-28 real 3D interim (ab3d_prune300, 4 states, ~iteration 115 stage-lines): states are B/C/D/E-enriched
(E-rich 15.1k E vs 6.0k C; C-rich 8.3k; D-rich 9.1k; one mixed), majority purity 0.424 vs chance 0.347.
Control from abinitio2D (ab3d_a200n) queued; scorer labels/interim3d.sh <run> gives interim confusions.
(ay) 2026-09-28 real ladder flattens: pruned(0.3) 800 cls -> NMI 0.154, purity 0.641, recall C/D/E 0.21/0.31/0.39
(600: 0.148; ceiling 0.193 = 80% reached). Next: 800 cls, prune 0.4, 6 rebalance rounds (rprune800r6).
(az) 2026-09-28 rprune800r6 (prune 0.4, 800 cls, 6 rebalance rounds, 30 A): NMI 0.164, purity 0.674 = the
ceiling's purity (NMI 85% of ceiling 0.193); recall C/D/E 0.23/0.30/0.41 (recall denominators include the
dropped particles, so recall is understated for pruned runs). Pruning depth + rounds still pay. Next rung:
prune 0.5, 1000 cls, r6.
(ba) 2026-09-28 real 3D from the pruned flex2D partition at iteration ~196: E-state 20.9k E / 5.5k C, C-state
12.0k C, D-state 13.0k D, 4th state = A/F/junk sink; majority purity 0.557 (chance 0.347), rising from 0.424.
abinitio2D-start control (ab3d_a200n) queued behind it.
(bb) 2026-09-28 rprune1000r6 (prune 0.5, 1000 cls, r6, 30 A; 66k particles kept, 6.7k unlabelled): purity 0.720,
NMI 0.183 (ceiling 0.193, 95%), recall over ALL labelled C/D/E 0.21/0.31/0.45. Pruning trades coverage for
purity. abinitio2D control on the same pruned subset (1000 cls) launched to separate pruning from flex2D.
(bc) 2026-09-28 DECISIVE 2D CONTROL on real 10076: abinitio2D on the identical pruned subset (66k ptcls, 1000
cls requested, 526 populated): purity 0.515, NMI 0.065, recall C/D/E 0.00/0.07/0.02 vs flex2D 0.720/0.183/
0.21/0.31/0.45 on the same particles. Pruning alone is not the gain. Fair 3D pair queued: ab3d_prune1000
(from rprune1000r6) vs ab3d_a2dprune1000 (from a2d_prune1000), behind ab3d_prune300 -> ab3d_a200n.
(bd) 2026-09-28 REAL 3D FINAL (ab3d_prune300, 4-state abinitio3D from the pruned flex2D 300-cls partition, 6.6 h,
10x8): E-state 21.3k E / 5.5k C / 4.5k D; C-state 12.4k C / 3.5k B; D-state 13.3k D; state 4 = A/F/junk sink;
majority purity 0.564 over labelled (chance 0.347); state volumes 6.0-15.5 A at FSC 0.143. Controls queued:
ab3d_a200n (from abinitio2D, running), then the fair pruned-subset pair.
(be) 2026-09-28 real 3D CONTROL from abinitio2D (ab3d_a200n) interim at iteration ~139 stage-lines: all four
states are B/C/D/E mixtures (best: 16.7k E with 10.1k C + 6.3k D), majority purity 0.355 = chance 0.347, vs
flex2D-start 0.424 at ~115 and 0.557 at ~196 -> class enrichment carries into real 3D states.
(bf) 2026-09-28 HETEROGENEITY MODE shipped (1f8651c2d): flex2D nstates>1 -> rebalance 0.8 x 4, labels marched to
FLEX2D_LPLABEL (30 A), FLEX2D_PRUNE (0.3) of the particles deactivated by evidence after stage 2, then 3
labels-frozen polish stages (20, 12, lpstop). nstates added to the flex2D UI. Tests: synthetic t4_mix5050
f_hetero2 (nstates=2, 120 cls), real rhetero800 (nstates=4, 800 cls, full project).
(bg) 2026-09-28 hetero mode validated on synthetic: t4_mix5050 nstates=2 ncls=120, no env flags: 108/120 pure,
majority 0.951, with polished averages, 523 s. Real rhetero800 (nstates=4, 800 cls, full project): pruned
39,570 after stage 2, polished median 8.38 A; label score pending.
(bh) 2026-09-28 real hetero mode (rhetero800: nstates=4 ncls=800, one command, no env): NMI 0.147, purity
0.623, recall C/D/E 0.19/0.27/0.39, polished median 8.38 A; vs the manual two-run pipeline 0.154/0.641
(in-run pruning uses stage-2 evidence). abinitio2D on the full set 0.031.
(bi) 2026-09-28 TRAP: hetero-mode ladder marched COARSER on PfCRT (default lpstart 15 A at mskdiam 180 < lp_label 30):
stages 15 -> 22.5 -> 26 -> 28 A. Fixed: lp_stages(1) = max(lpstart, lp_label) in hetero mode. PfCRT run
relaunched (pfcrt_from_scratch/flex2D_hetero; the inverted one kept as flex2D_hetero_inverted). PfCRT sieved
project: 343,106 active of 592,544, box 300 @ 1.035 A, mskdiam 180, post-sieve abinitio2D used ncls=200.
(bj) 2026-09-28 REAL 3D CONTROL FINAL (ab3d_a200n, from abinitio2D's 200-cls partition, 7.2 h): purity 0.418
(chance 0.347) vs 0.564 from the flex2D partition; its states: one D-rich (12.2k D / 5.8k E), one C+E lump
(17.1k C + 23.6k E), two sinks -> C and E separate only from the flex2D start. Fair pruned-subset pair
(ab3d_prune1000 running, ab3d_a2dprune1000 queued).
(bk) 2026-09-28 PfCRT flex2D hetero run (pfcrt_from_scratch/flex2D_hetero: sieved project, 400 cls, nstates=2,
mskdiam 180, 10x8): 45 min, pruned 102,932 -> 240,174 kept, polished median 6.75 A q 5.75/7.40, 3 empty
classes; all label stages at 30 A (lpstart = lp_label). FLEX2D_WRITE_TMPL=1 stacks cavgs_flex2D_iterNNN.mrc.
Comparison montage vs the post-sieve abinitio2D (6_abinitio2D, 200 cls) made as pfcrt_a2d_ranked.png.
(bl) 2026-09-28 EMAN2 CRASH on flex2D output = all-NaN class images: classes with an empty even or odd half get
a 0/0 even/odd FRC and the ML-regularized average becomes NaN (3 of 400 PfCRT classes, 270k NaN voxels).
Fixed in the class averager (simple_classaverager_restore.f90: such classes zeroed, FRC 0). Existing stacks
patched with pfcrt_from_scratch/sanitize.py (NaN -> 0); same pass queued for flex2D_selclust/pass2 outputs.
The cavgs_*_part*.mrc files are raw sums (not viewable) and the per-iteration cavgs_flex2D_iterNNN.mrc are
gridded templates at the crop box.
(bm) 2026-09-28 PfCRT follow-ups (400 cls, nstates=2): flex2D_selclust (10_selclust 162k -> 113.5k kept):
median 7.58 A q 7.06/9.71, 28 min; flex2D_pass2 (second pass on the pruned 240k -> 168k): median 7.06 A
q 5.97/7.58, 40 min; first pass on the full sieved set: 6.75 A q 5.75/7.40. Re-pruning/cluster selection
does not sharpen further. Both had 1-2 NaN classes (sanitized; averager fix 7d4234dfd for new runs).
(bn) 2026-09-28 PfCRT "useful?" test = particle SELECTION at matched count: refine3D_auto (start map
/mnt/beegfs/elmlund/pfcrt/9_refine3D_auto/recvol_state01.mrc, 300 @ 1.035, the 3.93 A run) on 40k particles:
top-40k by flex2D evidence vs top-40k by abinitio2D corr (sieve-active) vs random 40k sieve-active; chained
one at a time in pfcrt_from_scratch/subsets/ref_*. Pair montage (pairs/pfcrt_pairs_crop.png): abinitio2D
classes show more high-res core detail; flex2D partners smoother; eo-resolution is self-consistency, not detail.
abinitio3D_cavgs from both stacks (ab3dc_flex, ab3dc_a2d) + EMD-20806 download in flight for a map-based score.
(bo) 2026-09-28 PfCRT GOAL REFRAMED by the user: the deposited 17k is one state; they want the state the depositors
threw out of the full data. Matched-count refine tests cancelled (also failed: refine3D_auto wants e/o vols at
the autoscaled box 256). Running: subsets/st2_flex (abinitio3D nstates=2 on a random 100k of the flex2D-active
240k, flex2D labels) then st2_a2d (same from the abinitio2D partition, random 100k sieve-active). Score: state
maps vs EMD-20806 (pfcrt_from_scratch/emdb/emd_20806.map) after docking, and state-vs-state difference.
(bp) 2026-09-28 deliverables: guide doc/implementation_notes/flex2D_testing_and_development_guide_2026_09_28.md,
driver scripts/flex2D/autotest/flex2D_autotest.sh (commit e11250be0), session patch of 14 commits since
4df8fd872. PfCRT abinitio3D_cavgs maps docked to EMD-20806 (score_vs_emd.sh): flex 20.7 A at FSC 0.143 (first
pass, names collided -> rerun as ab3dc_score2 with distinct names incl. the old 16_abinitio3D map).
(bq) 2026-09-28 GOAL (user): flex2D exists to drive better class-restrained importance sampling in abinitio3D ->
heterogeneous initial models; the primary metric is the 3D per-particle state purity from the partition.
A GPU flex2D port is wanted (user is CPU-limited at home, has a big GPU) to speed testing; all generated code
must follow the .claude SIMPLE conventions. Guide sections 0, 10, 11 (commit after e11250be0).
(br) 2026-09-29 PfCRT 2-state abinitio3D from the flex2D partition (st2_flex, 100k random of the 240k active,
2.7 h): states 42.7k/57.3k at 18-19 A; docked FSC vs EMD-20806 38.8 A / 77.7 A at 0.143; state1 vs state2
differ beyond 78 A -> no meaningful second state at this resolution. abinitio2D-start control (st2_a2d) done,
scored the same way. Handoff zip in ~/Downloads/flex2D_handoff_2026-09-28.zip.
(bs) 2026-09-29 REAL 10076 FAIR PAIR, flex side: 4-state abinitio3D from the pruned 1000-class flex2D partition
(66k ptcls, 5.9 h): purity 0.710 (chance 0.347; 0.564 from the 300-cls partition): C-state 9.5k C / 1.5k B /
2.3k E; D-state 10.6k D; E-state 20.2k E; sink 0.5k F. abinitio2D control on the SAME particles
(ab3d_a2dprune1000) running. PfCRT st2 pair: both partitions give two 19 A states that match EMD-20806 only
at 39-78 A and differ from each other only beyond 78 A -> no discarded state resolvable by a 100k ab initio;
needs per-state refine3D on the flex2D labels or a better 2D partition.
(bt) 2026-09-29 10076 fair 3D pair DONE: 4-state abinitio3D on the identical 66k pruned particles: from flex2D's
1000-cls partition purity 0.710 (clean C, D, E) vs from abinitio2D's 0.503 (C+E lump, 2 near-empty states).
HEMOGLOBIN (EMPIAR-13532 H2_off, Cs 0.01, box 100 @ 1.3 A, 521,364 ptcls, deposited EMD-76806 C2, NICE workspace
/mnt/beegfs/elmlund/dibyajyoti/NICE_Projects/Cluster/Hemoglobin/.workspace_2; 8_abinitio2D = 100 cls mask 110):
run dir /mnt/beegfs/elmlund/afan/hemoglobin_flex2D/flex2D_n2_400 (job 53440283, nstates=2, 400 cls, mask 110).
TRAP: the 30 A label band is only rings 2-4 on a 130 A box -> label band should scale with box/particle size.
(bu) 2026-09-29 hemoglobin 30 A label run THRASHED and was cancelled at iteration 21: band [2,4] on the 130 A box,
257/309 of 400 classes below the floor, rebalancing moved 0.6M/1.1M ptcls, 90 % class / 83 % rotation changes
per iteration, templates = 30 A blobs. Label band must scale with the box (rings), not be a fixed 30 A.
(bv) 2026-09-29 e2display SEGFAULTS in flex2D run dirs although every .mrc header and pixel is clean; it works in
abinitio2D dirs. The difference is flex2D's binary intermediates (flex2D_*_partNN.bin, model/sigma2 .bin, algndoc
.simple) that EMAN2's browser probes. Fix used: intermediates moved to work/ (cancelled run), view/ with symlinks
to cavgs*.mrc for a running one (refresher job). The user wants fixes, not diagnosis loops: "just fix it".
Permanent fix to do: flex2D writes its part/model files into a subdirectory.
(bw) 2026-09-29 e2display SEGFAULT ROOT CAUSE: flex2D's per-iteration template stacks (write_templates, one
img%write per class) end up with header dmean = rms = 0 from iteration ~20 on (header stats not recomputed for
the whole stack); e2iminfo reads them fine but e2display divides by sigma 0 and segfaults. The earlier
"binary intermediates" theory was WRONG. Fix: sanitize.py rewrites dmin/dmax/dmean/rms from the pixels; a
refresher job (hb_fixhdr) does it for each finished stack of the running hemoglobin job. Permanent fix owed in
write_templates / the stack writer. EMAN2 on the box: /home/atifao/eman2-sphire-sparx/bin (not on PATH).
(bx) 2026-09-29 night: HEMOGLOBIN TOY SUITE. Maps /mnt/beegfs/elmlund/afan/hemoglobin_flex2D/toy/maps (make_hb_maps.py:
EMD-76806 from the NICE 3_scale/outvol.mrc resampled to 100 px @ 1.3 A; hb_T, hb_R6/hb_R12 = one alpha-beta dimer
rotated 6/12 deg about an in-plane axis, hb_dimer = one dimer). Diff power vs T: R12 0.015/0.116/0.407 at 30/12/8 A,
R6 0.004/0.034/0.143, dimer ~0.5 -> quaternary changes live at HIGH resolution, label band must be fine. Tests in
synth2D/hbA_TR12, hbB_TR6, hbC_Tdim (50k, SNR 0.05, junk 0.1, shift 5, a2d 100 cls); scripts run_synth_a2d.sh,
run_flex_hb.sh (MSK/NCLS/SNAPH env), ab3d_state_hb.sh. Build b2da8aa8c: lp_label = mskdiam/9 default, acc parts deleted.
Fixed hemoglobin real rerun: flex2D_n2_400_lpl12_fixed (a64028225 build, job 53450790).
(by) 2026-09-29 HANG BUG: evidence pruning emptied worker range 10 (synthetic junk is contiguous at the stack end) ->
the worker stops at startup ("no active particles in range", empty err/out, FAILED 0 s) and the master waits
forever at the next iteration. Fixed 237f816b4 (pruning restores the 20 best particles of an emptied range). Same
hang would hit any flex2D run whose selection empties a whole worker range. abinitio2D hemoglobin toy baselines:
state majority hbA T/R12 0.523, hbB T/R6 0.516 (chance 0.5), hbC tetramer/dimer 0.956 (chance 0.7); the rot/view
metrics are meaningless for the C2 particle.
(bz) 2026-09-30 HEMOGLOBIN TOY VERDICTS (50k, SNR 0.05, box 100 @ 1.3, mask 110; toy is CLEANER than real: a2d
classes 3.2 A median at 663/cls vs real 7.6 A at 5.2k/cls): T/R12 rotation: a2d 0.523, flex2D 0.549 (12 A labels)
/ 0.556 (8 A); CEILING from the true partition 0.552 (12 A) / 0.538 (8 A) -> a 12 deg dimer rotation is NOT
per-particle separable in 2D; flex2D already sits at the ceiling. T/R6 a2d 0.516, flex2D 0.533. Tetramer vs dimer
(70/30): a2d 0.956, flex2D 1.000 (198/198 classes pure). Compositional yes, quaternary no.
(ca) 2026-09-30 T/R12 at SNR 0.2: a2d 0.528, flex2D 0.556 -> NOT noise-limited: a rigid dimer rotation is absorbed
by the pose search (half the molecule rotated ~ whole molecule seen from a nearby view); needs 3D. B-factor 150 toy
still too clean (a2d 4.1 A at 900/cls vs real 7.6 A at 5.2k/cls); hbH = tetramer/dimer at SNR 0.01, bfac 250.
(cb) 2026-09-30 tetramer/dimer at bfac 150: purity a2d 0.961/0.979 (70/30, 90/10) vs flex2D 1.000/1.000, BUT flex
recall 0.77-0.79 vs a2d 0.88-0.98: the fixed 30 % pruning dropped ~10k GOOD particles. Evidence separates junk
cleanly (median - 3 MAD drops all junk, <=0.6 % good) -> aa087c7e2: hetero pruning default = robust cut
(FLEX2D_PRUNE=auto|fraction|0). hb_3q map (tetramer minus one subunit, 25 % mass, diff power 0.22) -> toys
hbI (50/50 bfac 150) and hbJ (80/20 SNR 0.01 bfac 250); f3_auto reruns on hbF/hbG.
(cc) 2026-09-30 COMPOSITIONAL HEMOGLOBIN TOYS WIN (build aa087c7e2, auto pruning): tetramer/dimer 70/30 purity/minority
recall a2d 0.961/0.88 vs flex 1.000/0.98; 90/10 0.979/0.78 vs 1.000/0.92; ONE SUBUNIT MISSING (25 % mass) 50/50
a2d 0.743 (recall 0.39/0.47) vs flex 0.998 (0.99/0.98); 80/20 SNR 0.01 bfac 250 a2d 0.902/0.53 vs flex 0.988/0.93.
Controls launched: a2d with 200 classes on hbI (a2d_200, run_a2d_hb.sh), 2-state 3D pair hbI ab3d_flex/ab3d_a2d.
Real hemoglobin on the final recipe: /mnt/beegfs/elmlund/afan/hemoglobin_flex2D/flex2D_v3.
(cd) 2026-09-30 CLASS-COUNT CONTROL: abinitio2D with 200 classes on hbI (subunit loss 50/50) majority 0.764 (100 cls
0.743) vs flex2D 0.998 -> flex2D's compositional win is not the class budget.
(ce) 2026-09-30 real hemoglobin flex2D_v3 (aa087c7e2: auto pruning, lp_label mskdiam/9=12.2, 400 cls): 69 min, auto pruning
removed only 89 of 521k (real junk is NOT an evidence outlier), eo median 4.06 A (q 3.61/4.64). Calibration: real a2d on
a random 50k subset, 100 cls: median 8.7 A at 500/cls vs toy SNR 0.01 bfac 250 6.2 A at 585 -> real is noisier;
hbK = T/3q 72/18 + 10 % hb_pair neighbour junk + 5 % noise, SNR 0.005, bfac 300.
(cf) 2026-09-30 shape classifier (hemoglobin_flex2D/cls_shape.py: 400 Fibonacci views/map, 2x downsample, 12 A, 36 rot x
FFT shifts) validated on hbI: flex pure classes T->T 91 %, 3q->partial(3q|3qB|dimer) 75 %; 3q vs 3qB indistinguishable
-> use as tetramer-vs-partial call with ~9 % false-partial baseline. Running on real flex2D_v3 and 8_abinitio2D classes.
(cg) 2026-09-30 GOAL RESULT ON TOY hbI: 2-state abinitio3D from flex2D partition per-particle purity 0.998 (clean T and 3q
states, junk to state 0) vs from abinitio2D 0.699 (both states mixed). Matched-noise toy hbK (SNR 0.005 bfac 300, 10 %
neighbour junk, 5 % noise): majority a2d 0.807 / flex 0.904; 3q recall 0.50 / 0.85; neighbour junk into pure classes
0.00 / 0.31. Real shape calls: 8_abinitio2D 45/100 classes "partial" (baseline false-partial ~9 %). Real 2-state
abinitio3D launched: hemoglobin_flex2D/ab3d_flex (from flex2D_v3) and ab3d_a2d (from 8_abinitio2D).
(ch) 2026-09-30 real shape calls particle-weighted: flex2D_v3 T 35 % / 3q+3qB 58 % / dimer 7 %; 8_abinitio2D 54 / 33 / 12 %.
Montages ambiguous (flex dimer calls look smaller; 3q calls look like asymmetric tetramers) -> not claimable. Guide
section 12 written. Awaiting real 2-state abinitio3D ab3d_flex / ab3d_a2d.
(ci) 2026-09-30 REAL hemoglobin 2-state abinitio3D (ab3d_flex 8.7 A, ab3d_a2d 7.8 A; score3d_hand.sh, 12 A docked corr):
non-mirror hand better; EVERY state best matches hb_T (flex 0.54/0.56, a2d 0.55/0.69); none dimer/3q-like; state-state
corr 0.44-0.60. Projections (real_states_proj.png) show tetramer-sized blobs. -> no evidence real hemoglobin carries
compositional species; its heterogeneity (if any) is T/R-like, which 2D cannot separate (ca).
(cj) 2026-09-30 hbK matched-noise 3D pair with RANDOM state seeding (standard abinitio3D nstates=2): flex partition purity
0.721 (both states T+3q mixed) vs a2d 0.812 (a2d state 1 = 80 % 3q, 3q recall 95 %). The 2D win (3q recall 0.85 vs 0.50)
did NOT carry to 3D here. MECHANISM: abinitio3D independent mode seeds ptcl3D states with gen_labelling 'uniform';
2D classes only stratify sample4update_class, they never seed states -> state-pure classes help only indirectly and
n=1 per arm is seed-noisy (hbI flex 0.998 vs a2d 0.699 went the other way). Direct lever = cavg_ini_ext route:
abinitio3D_cavgs nstates=N assigns a state per class, map2ptcls copies it to particles, abinitio3D cavg_ini_ext=yes
keeps those states (validate_cavg_ini_ext_states). Script synth2D/scripts/ab3d_cavgseed.sh; hbK cseed_flex / cseed_a2d.
(ck) 2026-09-30 CLASS-SEEDED 3D on hbK (synth2D/scripts/ab3d_cavgseed.sh + cseed_prep.py; SEED_ONLY=1 = class level only):
first try without class selection spent both states on pure-noise junk vs protein. With noise classes deselected (user-
style) and ptcl3D synced to the 2D selection: abinitio3D_cavgs seeds, 3q state purity/recall: a2d 0.94/0.50 in 3 of 4
(1 collapsed; nstates=3 leaves a state empty -> cavg_ini_ext refuses), overall 0.72-0.81; flex n2 0.63-0.64/0.66-0.79
(overall 0.81-0.83), flex n3_r2 0.75/0.90 + a pair-junk state (0.863), n3_r1 near-collapsed. a2d recall is capped at
the 50 % of 3q sitting in 3q-pure classes. TRAP fixed in 82c91b4a5 (snap built): flex2D pruning only zeroed ptcl2D
state; map2ptcls skips 2D-deselected particles, so their stale ptcl3D state leaked into cavg_ini_ext seeds; flex2D now
mirrors 2D deselection into ptcl3D. Particle stage (ab3d_cseed_cont.sh) running on cs_flex_n2_r2 vs cs_a2d_n2_r2.
(cl) 2026-09-30 REAL hemoglobin class-level seeds (hemoglobin_flex2D/hb_cseed.sh, cs_{flex,a2d}_n{2,3}): every state still
best matches hb_T at 12 A (0.52-0.68); extra states catch junk (a2d n3 state 2 = 30k-particle ring/donut in all three
projections; flex n3 state 3 = 3.8k blob); only a2d n2 state 1 (19k) is borderline (T 0.591 vs 3q 0.582). No
compositional species in real hemoglobin, consistent with (ci). Driver: flex2D_autotest.sh synth-cseed / synth-cseed-cont.
(cm) 2026-09-30 hbK class-seeded PARTICLE stage (nstates=2): a2d 0.841 (3q state 0.95 pure, recall 0.69) vs flex 0.842
(0.70 pure, recall 0.84) -> TIE; with random seeding a2d 0.812 vs flex 0.721. At real-hemoglobin noise flex2D's 2D
win does NOT carry into 3D by either route. hbI 0.998 vs 0.699 is CONFOUNDED (a2d arm kept 5k junk active, it fell
into the 3q state; flex had pruned it) -> fair replicate hbI_T3q_b150/rs_flex, rs_a2d (RANDOM_SEED=1: noise classes
deselected in both, standard random-state abinitio3D). Guide section 12.1 written.
(cn) 2026-09-30 real hemoglobin class-average look (radial dark spokes, bright off-centre wedges, lopsided backgrounds):
spokes are in 8_abinitio2D too (data property); flex2D's final make_cavgs ml_reg=yes reads the CANONICAL sigma2, which
flex2D never updates (its residual sigma2 lives only in memory + flex2D_sigma2_master.bin; canonical = calc_pspec
image-power bootstrap, ~1.7x residual at ring 2) -- BUT remaking v3 with ml_reg=no (flex2D_v3_remake, nparts=1; make_cavgs
needs nparts on the cline; rank_cavgs is simple_private_exec) looks the same, so the restoration is not the cause.
Hypothesis under test: flex2D classes group by NEIGHBOUR arrangement inside the 110 A mask (particle ~65 A);
flex2D_m80 (mskdiam=80, FLEX2D_LPLABEL=12.2 held, snap 82c91b4a5). 846cf9cb6: flex2D writes cavgs_flex2D_ranked.mrc
via rank_cavgs (abinitio2D ranking). Per-class res: a2d median 7.6 A at 5.2k/cls, flex 8.1 A at 1.37k/cls.
(co) 2026-09-30 DECISIVE: fair hbI replicate (noise classes deselected in BOTH arms, standard random-state abinitio3D,
NPARTS=9 because the contiguous junk block emptied part 10 and hung refine3D): a2d 0.994 vs flex 0.893. The old
0.998 vs 0.699 was entirely the junk confound. With matched selection, flex2D's state-pure 2D classes do NOT give a
better 3D split than abinitio2D's mixed ones on any toy (hbI 0.893 vs 0.994; hbK 0.72-0.84 vs 0.81-0.84).
Real v3 population churn: every rebalance (stages 1-4) found 144-189/400 classes below the 0.8*mean floor and moved
227k-317k particles; class choice treats templates as noise-free (as cluster2D), so big classes win. a2d pops min 3641
median 4875 max 12231 (100 cls); flex v3 min 1 median 1312 max 2399 (400 cls).
(cp) 2026-09-30 DRAIN diagnosed on hbK with the per-iteration population report (e59dc68a0): within a stage classes
COLLAPSE (it 23: min 1, median 30, max 2301 of mean 250; frac < mean/2 0.67); the stage-end rebalances hide it.
FLEX2D_TNOISE=1 (template-noise-corrected class choice, derived bias 4 sum_k w_k cbar2 k tw^2 sum_p G/(G+eps W)^2)
only halves it (median 62, max 3065). Suspect #2: the embedding shortlist (FLEX2D_SHORTLIST=8 in every run so far)
ranks templates by global-basis reconstruction, fitted to the big classes -> small classes fall off shortlists.
Arm d_sl0 (SHORTLIST=0) running. New switches: FLEX2D_TNOISE=<scale>, FLEX2D_LLR=1 (needs FLEX2D_CLS_RANK; ref
mu_c+U_c zc, zc fixed at pre-refinement pose). Residual sigma2 now committed before make_cavgs (e59dc68a0).
Dev runner synth2D/scripts/dev_flex.sh (SHORTLIST, MSK, NCLS, SNAPH env) -> <test>/dev_scores.txt. Crowded toy
hbL_crowd = hbK + 0-2 contact neighbours per particle (scripts make_hbL.sh, crowd.py, make_neigh_oris.py).
(cq) 2026-09-30 hbK A/B (dev_flex.sh, ncls 200, scores majority/NMI/rot/3q recall/pair recall): a2d 0.807/0.096/42.5/0.50/0.00;
base 0.891/0.198/24.5/0.83/0.24; TNOISE=1 0.886/0.195/27.0/0.79/0.23; SHORTLIST=0 0.896/0.206/25.1/0.82/0.28; PRICE=0.5
0.859/0.167/29.0/0.68/0.08; PRICE=1.0 0.845/0.151/33.5/0.62/0.08 -> forcing balanced populations HURTS (final
frac<mean/2 0.10 vs 0.14); the in-stage drain is not what limits purity. Keep all three off by default. Next arms:
CLS_RANK=1 vs +LLR, SPLIT_LDA (snap 1579d2660); hbN_noise (20k pure noise, hbK CTFs) = reference-bias test set.
(cr) 2026-09-30 REAL flex2D_c100 (ncls 100, mask 110, snap c28bda891, 44 min): pops min 1474 median 5099 max 7051 (a2d
3641/4875/12231); churn small (19/10/0 classes below floor). BUT the classes look WORSE than 8_abinitio2D at the same
populations: rotationally smeared, streaky, blurred low-res features, while REPORTING better class resolution (median
5.9 A vs a2d 7.6; eo median 5.2) -> inflated e/o FRC = reference bias / noise alignment in flex2D's pose refinement
(polish to 6 A with self-included templates is the prime suspect). The same flex2D beats a2d on toy poses (24.5 vs
42.5 deg) -> the toy lacks what real data has. hbN_noise bias test (N_a2d, N_flex) running.
(cs) 2026-09-30 LLR A/B on hbK (CLS_RANK=1 both): base 0.891/NMI 0.204/25.5 deg/3q recall 0.81; +LLR 0.859/0.161/28.4/0.62
(shift 2.47 vs 2.68 px) -> embedding-weighted references HURT, like GENREF before. PURE-NOISE bias test (hbN_noise, 20k
noise particles, 80 cls): class res a2d median 11.8 A (q25 10.0, best 7.6); flex2D median 10.0 A (q25 8.7) -> BOTH align
noise into "structure", flex2D more; flex residual/bootstrap sigma at ring 2 = 0.44 on pure noise (templates fit noise).
LOO in the E-step is dead code (l_loo_here never set true). FLEX2D_GOLD=1 (fc7798d3c): each half aligned/classified
against its own half templates x tw -> honest e/o FRC; N_gold running.
(ct) 2026-09-30 GOLD NOISE TEST (hbN_noise, 80 cls): reported class res median / q25 / best: a2d 11.8/10.0/7.6 A; flex2D
merged 10.0/8.7/5.7; FLEX2D_GOLD=1 21.7/18.6/11.8 (10 classes at 130 = none) -> gold removes the high-resolution noise
alignment and makes the class FRC honest; low-res noise blobs remain (inherent). Real gold run hemoglobin_flex2D/
flex2D_c100_gold (ncls 100, mask 110, snap fc7798d3c) launched; hbK d_gold (purity check) running.
(cu) 2026-09-30 GOLD on hbK at 200 cls (~125 ptcls per half): rot 45.1 vs 24.5 deg, majority 0.843 vs 0.891, 3q recall 0.62
vs 0.83, class res median 21.7 A (q75 130 = none) -> at low per-half populations the honest filter kills the signal
and alignment degrades. Testing the a2d-like regime: d50_base vs d50_gold (ncls 50, ~1000/cls); real gold c100 running.
(cv) 2026-09-30 REAL gold c100 (flex2D_c100_gold, 49 min): spokes/smear gone, honest class res 7.2-13 A (median ~10.8),
diverse shapes (two-lobed, rings, U) + 3 junk (lattice), BUT the final ml_reg classes are blurrier/less detailed than
8_abinitio2D; iteration-51 templates looked crisper than the final restoration. Toy hbK ncls 50: gold 0.822/rot 56.2
vs merged 0.867/40.8 -> gold costs alignment accuracy on the toy. Hybrid test: flex2D partition -> abinitio2D
cls_init=prev (a2dprev.sh, d_base_a2dprev) running. User: "looks better but a long way to go"; "there's probably a
more interesting way to use the latents".
(cw) 2026-09-30 ALIGNMENT-FREE prototype (synth2D/scripts alignfree_proto.py; af_feats.py + af_analyse.py + af_grid.sh
cached): CTF phase flip, low-pass 12 A, translation removed (centre = low-passed centroid, acf = autocorrelation),
polar + angular FFT, steerable PCA per m (J=6), cross-power + bispectrum invariants, 40 PCs. hbI (bfac 150, 50/50):
acf state AUC 0.916 view 0.643, k-means 100 majority 0.957 recall T 0.95 / 3q 0.96 (flex2D 0.998, a2d 0.743) --
works with ZERO alignment; centre 0.902/0.872 view-dominated, k-means 0.782. hbK (matched noise): centre AUC 0.866
(info present) but k-means majority 0.731 (T baseline 0.72), 3q recall 0.08; acf 0.681. Grid running (lp 12/18/25,
Wiener sPCA, knn-averaged invariants) in hbK/afg.
(cx) 2026-09-30 alignment-free on hbK (matched noise) quick grid: every unsupervised variant sits at the T baseline
(k-means majority 0.731-0.738, 3q recall <= 0.08) for lp 12/18/25, Wiener sPCA, knn-averaged invariants; local
detrending (minus knn mean) residual PCs state AUC 0.55; view patches x 2-split 0.738-0.745. Supervised state AUC stays
0.85-0.87 -> per-particle invariants of SNR-0.005 images carry the state weakly; nothing can be averaged without
alignment. Real-data check launched: 10076 (published labels, 40k subset, box 320 -> 96 px, lp 20) in
10076_flex2D/alignfree (af_real.py).
(cy) 2026-09-30 alignment-free on REAL 10076 (40k labelled, centre, lp 20): pairwise state AUC 0.58-0.65 among B/C/D/E,
0.81-0.86 vs 30S; k-means 100 NMI 0.012 (chance) -> per-particle invariants fail on real data too; alignment-free
parked. User: "proceed ... the model needs to differentiate view and composition". Plan: (1) FLEX2D_GOLD_LP=<A>
(e61ed9281) shared low band + independent high band; d_goldlp20 / N_goldlp20 running; (2) tied-loading mixture of
PPCA per class: K state templates + one within-class Gaussian subspace (view spread), state by marginal evidence;
prototype on flex2D nstates=1 view classes (d_views100 running).
(cz) 2026-09-30 FLEX2D_GOLD=1 FLEX2D_GOLD_LP=20 (shared band below 20 A, halves independent above; e61ed9281) = FIRST
change better on every count: hbK 200 cls purity 0.897 (merged 0.891, full gold 0.843), rot 26.8 (24.5 / 45.1), 3q
recall 0.82, pair 0.28, honest class res median 7.2 A; pure noise median 18.6 A (merged 10.0, full gold 21.7, a2d
11.8). Real run hemoglobin_flex2D/flex2D_c100_goldlp20 launched. TIED MIXTURE prototype (tied_mppca_proto.py on
flex2D nstates=1 ncls 100 view classes d_views100: purity 0.867, 3q 0.68): tied / untied MPPCA / 2-means K=2 per class
all leave purity at 0.867 (3q 0.68 -> 0.71) for lp 15/20, npc 12/20, q 1-3; BIC prefers 2 in ~52/96 classes even
though splits are not along state -> within-class mixtures do not find the state here. Toys: shared memory now.
(da) 2026-09-30 tied mixture on hbI (50/50 T/3q): flex2D nstates=1 ncls 100 view classes ALREADY 0.963 purity (a2d 0.743
at 100), recall 0.98/0.93; only 4 mixed classes (supervised AUC 0.957 inside them); tied/untied/2-means K=2 all leave
0.963 -> flex2D's own class search already separates composition from view at the class level on both toys; the
within-class mixture adds nothing. The remaining real-data problem is class-average quality (reference bias), not
the view/composition separation of the partition.
(db) 2026-09-30 ROOT CAUSE of weird real hemoglobin flex2D classes: PRE-CENTRING (FLEX2D_PRECENTER, default on for
random starts, added for the ribosome's +-20 px shifts). On crowded small particles the low-passed centroid is
pulled toward neighbours: rms centroid offset 15.5 px at the start; final shifts median 10.2 / q90 37.4 / q99 52.9 px
vs abinitio2D 4.3 / 8.5 / 13.5 (box 100). Per-iteration increments are small (<=3.4 px, then <1), template centring
converges (rms ~1 px) -> particles never return. Blur, smear, spokes and many two-lobed/U classes are plausibly
neighbour-centred particles. Test: flex2D_c100_goldlp20_nopc (FLEX2D_PRECENTER=0, shared memory nparts=1 nthr=64).
(dc) 2026-09-30 PRE-CENTRING CONFIRMED from logs (no new run): rms offset uncrowded toys 6.1-6.6 px (true 4.1), PURE
NOISE 15.9 px, REAL hemoglobin 15.5 px -> real pre-centring is noise-level. Default now OFF (88b444673,
FLEX2D_PRECENTER=1 opt-in). FAST DEV TOY hbK15k (15k random subset of hbK via selection nran, ncls 60, shared memory
32 thr): ~200 s per run. pc on/off there: 0.858/NMI 0.180/rot 42.7 vs 0.869/0.186/40.0 (off no worse).
state_recall.py now counts only active particles (subset-safe). User: "this test loop is far too slow".
(dd) 2026-09-30 15:01-15:14 the hemoglobin NICE workspace was REDONE by its owner: 6_extract re-extracted at box 150
(was 100; stacks rewritten in place), 8_abinitio2D deleted, new 7_abinitio2D (ncls 200, center=yes, sigma_est=group,
refine=prob_snhc, mskdiam=0). All earlier real flex2D results used the box-100 particles (runs finished before the
rewrite). User chose the box-150 set: copy hemoglobin_flex2D/box150/hb150.simple (+ hb150_50k.simple, nran 50000);
runs flex_full200 (distributed 10x8, gold lp20, pre-centring off default, 200 cls) and s_flex_pc / s_flex_nopc / s_a2d
(50k, 20 cls, shared memory 32 thr). Old-copy projects (hb.simple in flex2D_* dirs) now point at rewritten stacks.
(de) 2026-09-30 box-150 real 50k subset (hemoglobin_flex2D/box150, 20 cls): a2d crisp 7.0-8.1 A; flex2D from scratch
10-15 A blobs with pre-centring on OR off (so pre-centring is not the whole story on real data); flex2D STARTED FROM
a2d labels+poses (cls_init=prev, hetero mode) DEGRADES them to 6.7-13 A, most blurred -> flex2D's refinement damages
good real classes. Isolation running: s_fix_gold / s_fix_merged (FLEX2D_FIX_LABELS=1, nstates=1, cls_init=prev).
Crowded toy hbL: pre-centring on/off same purity 0.801 (uncrowded 0.886), shift err 5.1 vs 4.0 px. 15k sweep: best
ncls 200 0.884, GOLD_LP=15 0.880; soft/robust/SH3 neutral; a2d 0.806.
(df) 2026-09-30 USER: flex2D started from abinitio2D's classes (cls_init=prev) IMPROVES them ("probably the right way to
go") -- I had misread the honest gold resolutions vs a2d's merged ones. Recipe direction: abinitio2D -> flex2D
cls_init=prev. Cold-start gap: real hemoglobin defocus median 0.72 um (q10 0.52, min 0.20), Cs 0.01 vs toys 1.68 um,
Cs 2.7 -> CTF only ~0.2-0.5 at 30-60 A (toy ~0.9): stage-1 is SNR-starved on real data. Realistic toy hbR_realctf
(hbK molecules/views/shifts, real defocus draw, Cs 0.01; synth.simple = 15k subset, synth_full.simple = 50k) being
generated. Seeded toy arms hbK15k/seed_het, seed_fix; real frozen-label arms box150/s_fix_gold, s_fix_merged.
(dg) 2026-09-30 SPOKES = merged-reference noise alignment, settled on real data: from the same a2d start (box150 50k,
20 cls) flex2D labels-FROZEN + MERGED templates reproduces abinitio2D's spokes/speckle and claims 5.6-7.0 A; frozen +
GOLD (shared band <20 A) removes them (honest 7.5-12 A); labels free + gold looks like frozen + gold. User: "the spokes
are a problem in the abinitio and not in the flex so its certainly promising". Seeded toy (hbK15k from s_a2d): hetero
0.822 / 3q 0.50 (a2d 0.806/0.49; flex scratch 0.87-0.88/0.74-0.79) -> seeding keeps alignment gain, loses separation
on the toy. Template-centring sign bug (noisy tie locked wrong sign -> 53-83 px shifts) fixed in 561785d1e. Full seeded
run box150/seeded_full200 (from NICE 7_abinitio2D, 200 cls, hetero, gold) launched on 561785d1e.
(dh) 2026-09-30 SEEDED FULL SET (box150/seeded_full200: flex2D cls_init=prev from NICE 7_abinitio2D, 200 cls, hetero,
gold lp20, 561785d1e, 54 min): top classes 4.2-4.6 A honest, cleaner (few spokes) than a2d's, similar diversity
(paired/crowded classes kept). Variance maps (hemoglobin_flex2D/varmaps.py, per class, members CTF-flipped + aligned):
at 8 A flat (noise), at 20 A interior 0.92 / edge 0.94 / neighbour zone 1.06 / far 0.82 of far-field -> neighbours
add ~15 % over the interior, no internal high-variance patches -> a per-pixel noise model would change weights only
~10-15 %. Realistic-optics toy hbR (Cs 0.01, real defocus, 15k, 60 cls): flex scratch 0.809 / 3q 0.50 vs a2d 0.778 /
0.36 -- the toy advantage shrinks from ~0.07 to ~0.03 purity with real optics (cold start less SNR at low res).
(di) 2026-09-30 10076 SEEDED RECIPE (10076_flex2D/seed_a200n: flex2D cls_init=prev from abinitio2D a200n, 200 cls,
nstates 4, gold lp20 default, 561785d1e): published-label score purity 0.444 / NMI 0.069 / recall C/D/E 0.02/0.08/0.05
vs the a200n seed 0.385 / 0.031 / 0 and flex2D FROM SCRATCH at 200 cls/30 A 0.492 / 0.085 -> seeding more than doubles
a2d's NMI but separates less than a cold start, same pattern as the matched-noise toy (0.822 vs 0.87-0.88).
Gold+no-precentring are defaults since a02f08924; guide section 14 (f54cc01ec).
   FINAL (2599 s, 88 cores): 0.444 / 0.070; eo_res median 6.50 A (gold, honest) vs a200n 7.29 A (merged); the a2d's
   bottom ~25 classes are honeycomb noise classes at 38-70 A, the seeded run has ~none. Montages ~/Downloads/e10076_seeded.
(dj) 2026-09-30 seeded mode turned OFF all exploration (neigh_frac 1, no power sampling, no extremal) -> a2d's view
   partition is only refined locally. 16faac7a8 adds FLEX2D_SEED_EXPLORE=1 (keeps them on for cls_init=prev) and fixes
   FLEX2D_EXTR (was overridden by FLEX2D_ROBUST's env value). Toy A/B hbK15k/seed_ctl2 vs seed_explore (snap 16faac7a8).
   PfCRT seeded run: pfcrt_from_scratch/seed_a2d7 (7_abinitio2D sieved 343k, 200 cls, mskdiam 180, nstates 4, a02f08924);
   a2d7 reference montage eo median 9.71 A, mostly junk dots/micelle discs/spokes.
(dk) 2026-09-30 PfCRT SEEDED (seed_a2d7, 343k sieved, 200 cls, 3503 s): deposited-particle separation (pfcrt-deposited-
   good-labels) AUC 0.888 vs a2d7 0.870, rec@5/10/20% 0.45/0.65/0.82 vs 0.38/0.59/0.79, top classes 46/31/27/26 % good vs
   36/22/19/17 % -> flex2D separates the real transporter views from junk better; looks softer (gold eo median 10.7 A vs
   a2d merged 9.7 A), spoke-free in the good classes. Evidence pruning removed only 8 particles (junk forms classes,
   is not an outlier). Montages ~/Downloads/pfcrt_seeded/*_bygood.png. Full-592k seeded run seed_a2d6 (300 cls) launched
   on e02e9dbe6. SEED EXPLORE (16faac7a8) on hbK15k: 0.816 / 3q 0.46 vs ctl 0.822 / 0.50 -> exploration is not the gap.
   e02e9dbe6: cls_init=prev with ncls > seed's classes fills the empty templates by latent splits after the init
   iteration (65 of 65 filled on hbK15k 55 -> 120); arms seed_x120 vs s_a2d120 vs cold s_ncls120 (0.876 / 3q 0.78).
(dl) 2026-09-30 10076 SEEDED EXPANSION (e02e9dbe6; seed_x400 = a200n 200 -> flex2D ncls 400, nstates 4, 56 min): purity
   0.519 / NMI 0.109, recall B/C/D/E 0.03/0.15/0.19/0.24 (72 classes >70 %), all 131,899 kept; abinitio2D at 400
   (a400n, distributed, 68 min) 0.402 / 0.041 / no class >70 %; seeded at 200 0.444 / 0.070; pruned cold start at 400
   (rprune400, 105.5k kept) 0.584 / 0.125. Expansion is the first seeded-mode lever that moves real state separation.
   Launched: seed_x800 (200 -> 800) and 4-state abinitio3D from seed_x400 (ab3d_seedx400) beside ab3d_seed200
   (references ab3d_a200n 0.418, ab3d_prune300 0.564). PfCRT full set see pfcrt-deposited-good-labels.
   seed_x800 (200 -> 800, 74 min): 0.564 / NMI 0.130, recall B/C/D/E 0.06/0.20/0.27/0.37 (196 classes >70 %), all kept
   = the cold-start hetero 800 run's recall (rhetero800 0.19/0.27/0.39, purity 0.623 after pruning 39.5k). Seeded +
   expansion reaches cold-start state separation with the cleaner seeded averages. PfCRT 300 -> 600 pair running
   (seed6_x600 vs a2d6_600).
(dm) 2026-10-01 3D INTERIM (interim3d.sh counter; ab3d_prune300 was 0.424 at ~115 and 0.557 at ~196): ab3d_seed200 at
   199 = 0.402 (C+E lump, like a2d's final 0.418); ab3d_seedx400 at 169 = 0.475 (13.6k-D state, C+E lumped); ab3d_seedx800
   at 142 = 0.347 chance (lagging). User asked about steerable PCA: already tried as alignment-free invariants (fails at
   real SNR); remaining use = CWF denoising for the cold start's alignment, low priority vs expansion split direction,
   class-level junk score on PfCRT (deposited-label instrument), merging expanded classes before 3D.
(dn) 2026-10-01 PER-PIXEL MODEL (synth2D/scripts/pixmodel.py, oracle state means inside mixed flex2D classes, 5-fold CV,
   nearest-mean with isotropic vs per-pixel diagonal metrics): toy hbK15k seed_x120 (34 classes, lp 15) iso 0.738 / pix
   (within-state var) 0.756 / pixtot (class var, label-free) 0.758 / LDA 0.805; REAL 10076 seed_x400 (60 classes, lp 20)
   iso 0.967 / pix 0.964 / pixtot 0.962 / LDA 0.948 -> no real-data gain. And 10076 states are ~perfectly separable
   inside classes once the state means are known: the bottleneck is FINDING the means (unsupervised), not the metric.
   STEERABLE PCA denoise (synth2D/scripts/spca_denoise.py): DC must be kept out of the whitening (npow[0] from a
   mean-subtracted background is ~0 -> 1e23 eigenvalues) and each m restricted to rings k >= m/pi; hbR toy run pending.
   PfCRT softness: at a common 8 A filter with ml_reg=no rebuilds (fairfilt/) flex2D classes are STILL smoother than
   a2d's -> part of the softness is real alignment precision; EMD-20806 projection FRC test (fairfilt/vsmap.sh) pending.
(do) 2026-10-01 STEERABLE-PCA DENOISING = NEGATIVE. Working denoiser synth2D/scripts/spca_fb.py (Fourier-Bessel basis on the
   pixel grid, exact LS coefficients + exact noise covariance Pinv Pinv^T, per-k generalised eigenproblem, MP edge,
   Wiener; the polar-grid version spca_denoise.py is LOSSY: bilinear polar<->Cartesian round trip cc 0.836, exact only on
   the axes -> cross artefact). hbR toy 50k: 1930 components kept (flat ~1.4 bulk = neighbour signal). Classification on
   the denoised stack (reimport ctf=flip): abinitio2D 0.716 / NMI 0.005 / rot 135 / shift 78 px (vs raw 0.778 / 0.067),
   flex2D 0.751 / 0.044 (vs raw 0.809 / 0.116). The global-covariance Wiener filter pulls particles toward the dominant
   (view/neighbour) modes and erases the weaker state differences. Do not retry as a preprocessing step.
   PfCRT vs EMD-20806 projections (fairfilt/vsmap.sh: reproject nspace 3000 + match_stacks lp 12 + frc_vs_ref.py): FRC top40
   a2d 0.128 / 0.015 / -0.007 vs flex 0.103 / 0.005 / 0.007 at 30-12 / 12-8 / 8-6 A -> INCONCLUSIVE (sharpened protein map vs
   micelle-dominated averages; nothing beyond 12 A for either). Launched seed_a2d7_glp12 (FLEX2D_GOLD_LP=12) to test whether a
   wider shared band sharpens the classes without losing the junk separation.
   seed_a2d7_glp12 (FLEX2D_GOLD_LP=12, 58 min): junk AUC 0.888 (= lp20), eo median 9.71 (lp20 10.71), common-8-A montage
   ~unchanged (no helix stripes) -> the shared band is not the lever. flex2D shifts are continuous (pftc_shsrch_grad), so
   pose quantisation is not it either; the softness = half-template alignment noise in the 12-8 A band (merged mode is
   sharp with spokes, box150 s_fix_merged). Larger classes would sharpen the half templates.
(dp) 2026-10-01 3D FINAL ab3d_seed200 (from seed_a200n, 7.0 h): purity 0.404 = a2d's 0.418, C+E lumped (16.7k C + 23.4k E)
   -> the 200-class seeded partition gives NO 3D gain. ab3d_seedx400 / x800 pending (interim 0.475 / 0.347).
   ab3d_seedx400 0.484 (clean D state 14.8k D; C+E lumped 12.1k C + 23.2k E), ab3d_seedx800 0.448. SAME-SUBSET rescore
   (labels/same_subset.sh, the 88,818 labelled particles ab3d_prune300 kept): prune300 0.564, a200n 0.444, seed200 0.429,
   seedx400 0.516, seedx800 0.469 -> expansion helps 3D (best at 400), pruning is ~0.03 of the gap, ~0.05 remains: only the
   cold-start partition splits C from E in 3D. 800 classes is worse than 400 for the class-restrained 3D sampling.
   TRAP: a ptcl dump regex r"state=(...)" matches mi_state first; parse with dict(findall)['state'].
(dq) 2026-10-01 SPLIT DIRECTION (791da4654: fill waits for class latents; FLEX2D_CLS_RANK=4 FLEX2D_SPLIT_CLS=1
   FLEX2D_SPLIT_DEFLATE=2): toy seed_x120_cls 0.834 / 3q 0.64 (plain 0.828 / 0.57; the stage-1 rebalance filled the empties
   first, with the same class-latent split), 10076 seed_x400_cls 0.516 / 0.106 (plain 0.519 / 0.109) -> class count matters,
   split direction does not. 3D replicates launched (ab3d_seedx400_r2, ab3d_prune300_r2) to see if the 0.516 vs 0.564
   same-subset gap is beyond abinitio3D run-to-run noise.
(dr) 2026-10-01 3D REPLICATE INTERIM (~iteration 190 of the interim counter): ab3d_seedx400_r2 0.549 (C 16.2k / E 20.8k / D 14.4k
   in separate states, all particles kept) vs ab3d_prune300_r2 0.554 (pruned) -> the run-1 gap (0.484 vs 0.564) was mostly
   abinitio3D run-to-run variation: n=1 per arm cannot rank 3D arms on 10076. Finals pending (same_subset.sh rescoring queued
   in the b1dlcwyn1 watcher).
(ds) 2026-10-01 PER-PIXEL ALIGNMENT on PfCRT = NEGATIVE (seed_a2d7/pixalign.py + pixalign/run.sh: 40 top-good + 40 random seeded
   classes, box 128 @ 2.43 A, gold halves each realigned to its own half mean, +-6 deg x +-3 px, 2 rounds, FFT-weighted score):
   eo FRC 30-12 A top uni 0.522 / pix 0.512, random 0.371 / 0.345; 12-8 A equal; pix better in 20/40 and 19/40. The 1/var weight
   map is nearly flat at 20 A (p90/p10 1.13-1.19): at single-particle SNR the per-pixel variance is ~all noise, micelle and
   neighbour variability add a few %. A principled per-pixel noise model cannot reweight; only an a-priori focus mask
   (protein / structural-variance mask) could. TRAP: per-micrograph stacks (3300) -> open-read-close, no memmap cache.
(dt) 2026-10-01 2D FOCUS MASK on PfCRT = POSITIVE (seed_a2d7/maskalign.py, same gold harness as pixalign): per-class masks from
   the class average alone (lp 20, sign-normalised; ENVELOPE = Otsu of the circle, fill + dilate 2 px, soft; CORE = Otsu inside
   the envelope; focus = core 1 / rim 0.15). Top-40 good classes, eo FRC in the CORE at 12-8 A: uni 0.486, env 0.541 (+0.063,
   33/40), focus 0.531 (32/40); core FRC0.143 8.1 -> 7.7/7.6 A; circle FRC0.143 15.5 -> 11.9 A. Random 40: +0.03 (25/40).
   The tight envelope (solvent removal) carries most of it; rim down-weighting costs rim FRC (env 0.304 -> 0.246).
   Montage ~/Downloads/pfcrt_focus/pfcrt_focus_masks.png. ab3d_prune300_r2 died at 16330 s on a BeeGFS 'Communication error'
   (sigma2_state.bin unreadable); its last interim 0.560 vs seeded 0.556 stands.
(du) 2026-10-01 BUG: flex2D's polish 'labels frozen' (and FLEX2D_FIX_LABELS) only stopped rebalancing/k-means/prices; the E-step
   class search (ctx%l_cls_search, default on) kept re-assigning the argmax class. So every hetero run's polish and the
   hemoglobin s_fix_* arms moved labels. 5b60817e2: FLEX2D_POLISH_FREEZE=1 stops the class search in frozen stages
   (implied by FLEX2D_ENVMASK); the log says which. FLEX2D_ENVMASK (a349c5d6a, automask2D of the class template rotated into
   the particle frame by rtsq(-e3), applied in prepimg4align before the circular mask; build_batch_particles2D cls_masks=)
   first A/B without the freeze: mask area 0.06 of the box, eo median unchanged 11.10 A, junk AUC 0.888 -> 0.857 because the
   class search ran on per-class-masked particles. Rerun with freeze: seed7_envmask/envbase (PfCRT), box150 hb_x400 / hb_x400_mask.
(dv) 2026-10-01 HEMOGLOBIN with the current recipe (box150/hb_x400, hb_x400_mask: NICE a2d 200 -> flex2D ncls 400, nstates 2,
   FLEX2D_POLISH_FREEZE=1, 5b60817e2, 36-40 min): top-60 ranked classes mostly clean single molecules, few paired, weak spokes
   vs abinitio2D's speckled/spoky/paired top 60; envelope mask (area 0.116 of the box) changes nothing visible, eo median 8.86 A
   both. Montages ~/Downloads/hb_improved/. PfCRT frozen A/B (seed7_envmask / seed7_envbase) pending.
(dw) 2026-10-01 PfCRT ENVMASK A/B with frozen polish (seed7_envbase / seed7_envmask, 5b60817e2): junk AUC 0.887 / 0.880, eo median
   11.10 both; final-class core FRC (core_frc.py, top-40 good classes, noreg even/odd) 12-8 A 0.892 / 0.875, 8-6 A 0.646 / 0.674;
   circle FRC0.143 12.4 / 11.3 A -> mixed, within run-to-run noise: labels agree only 61 % although the mask acts after labels
   freeze -> distributed flex2D is not run-reproducible (random sigma2 bootstrap subset; see flex-pca-ab-needs-simple-seed).
   The prototype's +0.06 core gain does not carry into flex2D. Envelope mask = not a default. Best hemoglobin run: box150/hb_x400.
(dx) 2026-10-01 TEMPLATE-FREE POSES (user idea: no class average in alignment). synth2D/scripts/sync_rot.py = angular
   synchronisation inside flex2D classes: FB coefficients, pairwise relative rotation via one FFT over k, leading eigenvector
   of H_ij = e^{i th_ij}. hbK15k seed_x120 (106 classes, 60 partners): rot spread 35.1 deg (lp 15) / 35.6 (lp 10) vs the
   run's template-aligned 33.8 -> no-average rotations are as good, GIVEN the run's shifts with sign -1 (centre = shift by
   +x,+y in fftfreq convention); wrong centring -> 120 deg (random 104). Next: shift-free version (rotation from the
   autocorrelation mod 180, Z2 flip sync, translation least squares). hb_x400 3D (ab3d_x400, 2 states 267k/237k, 1.4 h):
   both states best match hb_T (0.476 mirr / 0.508), corr below the old box-100 runs (0.54-0.69); no compositional state.
   SHIFT-FREE attempt 1 (sync_se2.py: rotation mod 180 from the per-particle AUTOCORRELATION in a padded FB basis, Z2 flip sync,
   translation least squares): rot spread 105.7 deg = random on hbK15k -> the autocorrelation is quadratic in the image,
   noise x noise swamps it at toy SNR; dead end. Attempt 2 (sync_se2b.py): joint pairwise rotation (5 deg grid) x shift (FFT,
   +-8 px) search, angular sync, then translation LS in the synchronised frame -- running.
(dy) 2026-10-02 PURE NOISE, template vs SYNC (realsync_v2.py on hbN_noise N_flex / N_goldlp20 partitions, 30+30 classes, gold halves):
   half-average power 30-8 A over the random-rotation average of the same particles: template poses 2.79-2.90 (merged) / 2.53-2.94
   (gold lp20), SYNC 0.84-0.89; even/odd FRC in the circle at 30-12 A: merged 0.75-0.77, gold 0.36-0.37, sync 0.12-0.15.
   -> template alignment (incl. gold-standard) invents ~3x noise power per half; angular+translation synchronisation (sparse random
   partner graph, spectral) does NOT: noise pair angles are loop-inconsistent, the eigenvector cannot lock in. The clean
   'no Einstein-from-noise' figure for a template-free pose engine. ab3d_seedx400_r2 stalled 11.5 h at iteration 70 after the
   cluster outage (workers lost) -> cancelled; its interim 0.556 vs prune300_r2 interim 0.560 stands.
   Toy solvers (sync_se2c.py, hbK15k seed_x120, curves computed once): template 33.8 deg / 1.27 px; spectral 34.5 / 1.28-1.29;
   weighted 34.4; spectral + ICM (discrete 5-deg coordinate ascent against partners) 33.9-34.0 / 1.26; 150 partners = 60 ->
   both at the ~34 deg floor of within-class view spread; pose error cannot separate them on this toy. Ground-truth FRC test
   (truthfrc/: class averages from template vs sync+ICM poses, matched to reprojections of hb_T/hb_3q/hb_pair, FRC vs truth)
   launched; real hemoglobin realsync (box150/realsync) still running.
(dz) 2026-10-02 GROUND-TRUTH FRC (synth2D/truthfrc: truthfrc_avgs.py + match_stacks vs reprojections of hb_T/hb_3q/hb_pair +
   truthfrc_score.py), hbK15k seed_x120, 106 classes, same particles: template-pose averages FRC vs truth 30-12 A 0.191 / 12-8 A
   0.013; SYNC+ICM averages 0.241 / 0.020; per class sync better in 76/106 (median +0.063) at 30-12 A, 56/106 at 12-8 A.
   Absolute FRCs low (phase-flipped, not CTF-amplitude corrected; mixed classes) but the comparison is matched. Template-free
   averages are closer to the truth -> consistent with templates injecting noise structure (pure-noise test).
   REAL hemoglobin gold-halves FRC (box150/realsync, hb_x400, 40 largest + 40 random): template core 30-12 0.944 / 12-8 0.719 vs sync
   0.922 / 0.634; circle 12-8 sync 0.230 vs template 0.044 -> NOT interpretable alone: the pure-noise baselines are template core
   0.90 / 0.66 (gold shares the <20 A band) vs sync 0.41 / 0.17. Real-data referee = FRC vs projections of the known tetramer
   (toy hb_T map is hemoglobin at the same 1.3 A/px): box150/truthfrc (80 largest classes, cap 1000, cropped to box 100) running.
   CROWDED toy (hbL_crowd/truthfrc, L_nopc partition, 194 classes): FRC vs truth template 0.179 / 0.003 vs sync 0.216 / 0.020 at
   30-12 / 12-8 A; sync better in 115/194 and 114/194 classes -> the template-free gain holds with neighbours in the box.
   REAL HEMOGLOBIN ground truth (box150/truthfrc, 80 largest hb_x400 classes, cap 1000, vs 2000 projections of the hb_T map):
   template 0.234 / 0.004 vs SYNC 0.169 / 0.011 at 30-12 / 12-8 A; sync better in only 27/80 at 30-12 -> the toy gain does NOT
   carry to real hemoglobin; the sync halves' 12-8 A agreement (0.23) matches no truth = shared artefact (neighbours / CTF rings).
   Suspect near-focus optics (median 0.72 um, Cs 0.01: single particles carry little 15-60 A contrast; pairs starve, templates
   do not). Next: realistic-CTF toy hbR (r_flex partition) to reproduce, then CTF-weighted pair scores / higher band / partners.
(ea) 2026-10-02 SYMMETRY TRAP: every hemoglobin abinitio3D of mine (ab3d_flex, ab3d_a2d, box150/ab3d_x400) ran pgrp=c1; the dataset
   owner's NICE 11_abinitio3D / 12_refine3D_states use pgrp=c2 (alpha2beta2 tetramer). User spotted it in the maps. Reruns:
   box150/ab3d_x400_c2 (hb_x400 partition) and ab3d_a2d_c2 (box-150 NICE a2d seed control), nstates 2, same scoring.
   Use pgrp=c2 for any hemoglobin 3D.
   REALISTIC-OPTICS toy (hbR_realctf/truthfrc_lp15/lp8, r_flex 55 classes): template 0.219 / 0.010 vs sync 0.242 / 0.038 (lp 15,
   32/55 better) and 0.245 / 0.032 (lp 8, 34/55) -> near-focus optics alone do NOT explain the real-hemoglobin loss. Next
   hypothesis: class size (real classes cap 1000 vs toy ~140-270: template SNR grows with n, 60-partner sync does not) -> real
   hemoglobin rerun with cap 150 (box150/truthfrc_n150).
   Real hemoglobin cap 150 (truthfrc_n150): template 0.230 vs sync 0.173 (17/80) -> class size is not it (template poses still came
   from flex2D's full-class alignment). Next hypothesis: NON-WHITE real noise (ice/carbon gradients dominate raw pair correlations,
   templates cancel them) -> truthfrc_avgs.py whiten=1 (radial PS whitening + 60 A high-pass for the pair comparisons only),
   box150/truthfrc_wh running.
   Whitened pairs (truthfrc_wh): sync 0.187 vs template 0.234 (28/80) -> noise colour is a small part. Next: whiten + matched
   filter (x|CTF| for the pair comparisons) + lp 8 (truthfrc_whm8, truthfrc_avgs.py args ... 8 80 1000 100 1 1).
(eb) 2026-10-02 HEMOGLOBIN 3D with C2 (box150/ab3d_x400_c2 vs ab3d_a2d_c2, nstates 2, best hand vs hb_T): flex x400 states 0.572
   (228k) / 0.421 (276k); a2d control 0.492 (242k) / 0.566 (262k). C1 x400 was 0.476 / 0.508. Best states tie (0.572 vs 0.566),
   flex's weaker state is weaker; every state matches hb_T best (no 3q state) -> 2 states split better/worse, not composition.
   Maps ~/Downloads/hb_3D_c2/. Mirror hand best in all four.
   Real hemoglobin sync variants vs template 0.234 (30-12 A, FRC vs tetramer): plain 0.169, cap150 0.173, whitened 0.187,
   whitened + matched filter (x|CTF|) + lp 8 0.135 (11/80) -> higher band hurts single-particle pairs on real data. Last try:
   VDM-style cycle-consistency edge weights (truthfrc_avgs.py robust=1; box150/truthfrc_whvdm). If it loses, close the
   real-data template-free branch: bias-free (noise test) and toy wins are real, real hemoglobin accuracy lags templates.
(ec) 2026-10-02 TEMPLATE-FREE BRANCH CLOSED for real data: VDM cycle-consistency edge weights (truthfrc_whvdm) 0.170 vs template
   0.234 (25/80) -- no better than whitening alone (0.187). Summary: sync is bias-free (noise power 0.86x vs ~2.8x) and wins on all
   three toys (clean +0.06, crowded +0.02, realistic optics +0.03 vs truth), but on real hemoglobin it loses in every variant
   (0.135-0.187 vs 0.234). Not explained by optics, class size, noise colour, band, matched filter or outlier edges; the toys miss
   something about real single particles. Scripts: synth2D/scripts/{sync_rot,sync_se2b,sync_se2c,realsync_v2,truthfrc_avgs,
   truthfrc_score}.py.
(ed) 2026-10-02 SPLICEOSOME 10180 46k (10180/flex2D_46k/{a2d50,flex_x150,cls3}; project = 4_abinitio3D copy, mskdiam 400):
   abinitio2D 50 cls 29 min; flex2D seeded x3 (ncls=150 nstates=4 FLEX2D_SHORTLIST=8 FLEX2D_POLISH_FREEZE=1, 77e76ac94) 23 min.
   Classes <=14 A: a2d 14 (39 % ptcls) vs flex 60 (48 %); <=20 A: 14 (39 %) vs 98 (69 %). Matched cls_split ncls=3 on a2d50
   (150 subclasses, flexcls snap, 4.4 min). Vs flex_pca 5-state labels (flexcls/10180/prep/ptcl3D_states5.txt): a2d 0.449/0.114,
   cls_split 0.481/0.121, flex2D 0.602/0.213 (majority/NMI). CAVEAT: those labels are VIEW-DOMINATED -- 50 view cones from the 3D
   poses alone score 0.623/0.238, so purity mostly measures view separation. Beyond-view readout (view_ctrl.py, I(K;T|V) minus a
   within-cone shuffle null) at 50/200/500 cones: a2d 0.072/0.057/0.043, cls_split 0.083/0.052/0.039, flex2D 0.171/0.110/0.078 nats
   -> flex2D carries ~2x the state info beyond view; cls_split adds nothing beyond its parents. Scripts in 10180/flex2D_46k/
   (score_cls.py, view_ctrl.py, mont.sh). Cluster numpy lacks default_rng (use RandomState).
(ee) 2026-10-03 READOUT + REPLICATES (8a1251a1e, flex2D_repro.f90 = nearest-class cross-half repro, cls2D nn/nncc/repro60/
   repro30/repro15/repro): 10180 46k flex_x150_r2 bands >60/60-30/30-15/20-10 = 0.951/0.893/0.535/0.146 (20-10 = only
   independent band, 1 % of classes > 0.5); cls_split-seeded flex2D identical 0.957/0.897/0.538/0.146. Scores vs states5:
   rep2 0.598/0.211, cls_split seed 0.608/0.215 (rep1 0.602/0.213). Partition agreement: rep1 vs rep2 NMI 0.554 (overlap 0.53),
   a2d-seed vs cls-seed 0.376, cls_split partition vs flex2D from it 0.252 (overlap 0.13) -> seed split irrelevant, re-sorting
   is the gain. CODE FACT: the terminal dense pass (cmd ~L1277-1310) rebuilds ctx with class search ON and relabels every
   particle at the final band -- FLEX2D_POLISH_FREEZE does not reach it; delivered labels are 8 A argmax, not label-phase labels.
   abinitio3D_cavgs with force_lp_range lpstart=30 never converges (EVEN/ODD ANGULAR DISTANCE 88.5); use lpstart=70.
   Default flex2D search objective = cluster2D euclid to the class template (FLEX2D_LATENT_EV=0); latent = splits/shortlist/corr.
   Manuscript draft: doc/manuscripts/flex2D_cls_split_2026_10_03/main.tex (+ figures/, make_figures_*.py), tectonic-compiled.
(ef) 2026-10-03 THREE-SET CAMPAIGN (/mnt/beegfs/elmlund/afan/flex3sets/<set>/{flex,cls}, script flex3sets/run_set.sh <dir> <src>
   <mskdiam> <flex|cls> [ncls]; flex = snap 8a1251a1e seeded nstates=4 SHORTLIST 8 POLISH_FREEZE; cls = flexcls snap ncls=3).
   rnasep: src rnasep_states/afm_cryoem.simple (932,496 active, only 112 populated of max label 498) -> flex500 (ncls=500;
   ncls must be >= max label), msk 200. trpm4c: src processed_test_datasets/trpm4_correct/20_abinitio2D/trpm4_correct.simple
   (162,326 active, 150 cls; the trpm4_sieve/8_extract project has the WRONG smpd) -> ncls 450, msk 250. e10028: src
   10028/flex2D_tests/ref.simple (83,626 active, 239 cls) -> ncls 600, msk 300. Also running: 10076 seed_new_x600
   (10076_flex2D/seed_new_x600, 600 cls from flexcls/10076/a2d_new) and NOT x900. Manuscript v3 = integrated maths, no
   proofs (main.tex; v1/v2 kept as main_v1.tex/main_v2.tex).
(eg) 2026-10-03 10076 DROP DIAGNOSED: seed_new_x600 (0.392/0.032) was seeded from flexcls/10076/a2d_new = protein-DARK stack.
   Same recipe from a200n (_neg stack): R1 old build 8a1251a1e 400 cls 0.523/0.109, R2 fixed build ca660d327 0.524/0.106
   (earlier seed_x400 0.519/0.109) -> build and terminal-pass fix are neutral; the seed's contrast was the whole drop. The
   label-free test flagged it: independent 20-10 A band median 0.09 (dark seed) vs 0.48-0.51 (good runs). Terminal-pass fix
   (ca660d327, keep label-phase classes under POLISH_FREEZE) also neutral on the spliceosome: 0.596/0.208, beyond-view
   0.169 nats @50 cones (was 0.602/0.213, 0.171). Label-band sweep running: 10076 dg_lp20/dg_lp45_fix400, spliceosome
   flex_x150_lp30/lp65 (FLEX2D_LPLABEL via sbatch --export).
(eh) 2026-10-03 RNase P (conformational; labels = 25_refine3D_states final states, NOT view-dominated: view-only 0.351/0.027,
   chance 0.304): a2d 112 cls 0.383/0.038, cls_split 336 0.388/0.037, flex2D 500 (flex3sets/rnasep/flex500, old build)
   0.427/0.054; beyond-view excess @50/200 cones a2d 0.132/0.143, cls 0.141/0.152, flex 0.212/0.227 nats (~1.6x).
   Readout indep band: rnasep 0.47, trpm4c 0.31, e10028 0.42, 10076(a200n) 0.48-0.51. TRPM4 old-build flex2D left 40/450
   classes with 1-4 ptcls (terminal relabel drain) -> rerun on ca660d327 (trpm4c/flexfix). cls_split bands: rnasep
   0.97/0.94/0.90, trpm4c 0.82/0.80/0.70, e10028 0.46/0.67/0.55; e10028 cls_split turned two sharp parents into blobs (bug?).
(ei) 2026-10-03 LABEL-BAND SWEEP + RULE: 10076 (a200n, 400, fix build) lp20 0.444/0.068, lp31 (d/9) 0.524/0.106, lp45 0.520/0.089;
   spliceosome (150) lp30 0.567/0.190 [beyond-view 0.150], lp44 0.596/0.208 [0.169], lp65 0.622/0.221 [0.185] (monotone).
   Within-class difference FRC (flex3sets/diff_frc.py on cls_split even/odd, best sibling pair, median over parents) falls
   below 0.5 at: spliceosome 36 A, 10076 14.5 A, TRPM4 13.4 A, RNase P 8.7 A (curve rises again at high res = shared-ref
   bias of cls_split's cluster2D round), 80S estimator broken (noisy first shell). Candidate rule: optimal label band ~2x
   that crossing (65/36, 31/14.5). PREDICTION written before the RNase P sweep finished: optimum ~17 A, so lp15/lp22 beat
   lp33 (runs flex3sets/rnasep/flexlp15, flexfix (22), flexlp33, all ca660d327).
(ej) 2026-10-03 TRPM4 tiny classes are NOT the terminal relabel: fixed-build rerun (trpm4c/flexfix) still has 51/450 classes
   with <=4 ptcls. Cause: stage-2 re-seeding of 69 emptied templates from the parts' worst-particle pools (junk), then drain;
   rebalancing needs the largest class >= 2x floor (578) but max is 569, so nothing refills them. Fix idea: fold classes
   below ~10 ptcls into their best neighbours at the end of the label phase.
(ek) 2026-10-03 RNase P sweep CONFIRMS the predicted rule: lp15 0.429/0.055 beyond-view 0.211/0.225 (46 pure), lp22 (d/9)
   0.425/0.052 0.206/0.222 (27 pure), lp33 (d/6) 0.383/0.028 0.134/0.149 (0 pure). Rule "label band ~ 2x the resolution where
   the within-class (cls_split sibling) difference FRC falls below 0.5": spliceosome 72 (best tested 65, monotone), 10076 29
   (best 31), RNase P 17 (best 15). A fixed mask fraction cannot work: d/6 best for spliceosome, worst for RNase P.
   Turnover check queued: spliceosome flex_x150_lp72 / lp95.
(el) 2026-10-03 spliceosome turnover: lp72 0.619/0.219 beyond-view 0.182, lp95 0.603/0.206 0.170 (lp65 0.185) -> peak 65-72 A,
   rule predicted 75 A (crossing 37.3 A with the interpolated estimator; 10076 15.4, RNase P 9.0, TRPM4 13.6, 80S 22.3).
   Figure doc/manuscripts/.../make_figure_band.py (fig_band.pdf).
(em) 2026-10-03 TRAP: 10028/flex2D_tests/ref.simple ptcl2D labels (240 classes) do NOT match its registered cavgs stack
   (cavgs_iter027, 300 classes) -> montage parents were wrong; regenerated seed cavgs with make_cavgs (flex3sets/e10028/a2d).
   Real cls_split anomaly left: one sharp 80S parent (res 5.9) -> all 3 subclasses ring-shaped (rotationally averaged), likely
   the cluster2D in-plane round randomising rotations for that class.
(en) 2026-10-03 NOT (805k, not_flex2D/x900, 900 cls, 8a1251a1e, 6.0 h): particles in classes <=5 A: a2d300 0.61 (137 cls) vs
   flex2D 0.68 (596 cls); readout indep band median 0.765 (70% > 0.5). Montage flex3sets/not_3way.png (cls_split from
   flexcls/not/x3, 300 parents). Paper updated with RNase P/TRPM4/80S montages, RNase P state table, label-band rule
   (fig_band, tab:band, eq:bandrule), 10076 seed-contrast note.
(eo) 2026-10-03 MATCHED-COUNT CONTROLS (a2d_controls/<tag>, abinitio2D at flex2D's ncls, snapshot ca660d327; score_ctrl.sh,
   res_dist.py, clus_ctrl.sh, post_ctrl_remote.sh dispatcher). SPLICEOSOME: abinitio2D 150 EQUALS flex2D 150 on the state
   labels (0.610/0.216 vs 0.602/0.213; beyond-view 0.179/0.116/0.089 vs 0.171/0.110/0.078 at 50/200/500 cones) and is
   SHARPER (54% of ptcls in classes <=10 A vs 30%; flex only wins <=20 A, 70% vs 63%). The paper's "2x abinitio2D" on the
   spliceosome was a 50-vs-150 class-count artefact. 10076 600: a2d 0.413/0.045 vs flex 0.559/0.122 (flex wins, compositional).
   RNase P 500: beyond-view a2d 0.177 vs flex 0.212 (1.2x, not 1.6x). EVERY set: a2d sharper (<=8 A frac + best class),
   flex fewer blurred (<=20 A). TRPM4: flex's 97% <=20 A may be junk absorbed (a2d isolates carbon edges). Empty classes
   must be dropped before cluster_cavgs (80S a2d600: 157 empty formed a blank cluster). Paper has tab:matched + revised
   abstract/discussion; NOT 900 + spliceosome ab3d of a2d150 (10180/flex2D_46k/ab3d_a2d150) pending at write time.
(ep) 2026-10-03 TILT-MARGINAL PROTOTYPE (FLEX2D_CLS_MARG=1, needs FLEX2D_CLS_RANK; 762906e0e; d1a20eff6 adds per-class pose-tangent
   deflation of Ucls): class evidence + 1/2 b'(A_c+I)^-1 b - 1/2 logdet(A_c+I) at shortlist + refined poses. hbK15k seed s_a2d ->120
   (synth2D/hbK15k/tilt, toy_arm.sh SNAP=): base 0.828/0.123 rec .93/.57/.01 rot 35; CLS_RANK=2 alone IDENTICAL to base; r2m 0.818 rot 48;
   r4m 0.812; r2md (deflated) 0.813 rot 46 -> marginalising the learned class basis HURTS (classes broaden in view 60->64 deg and the
   basis absorbs the state difference of mixed classes). A learned per-class subspace cannot tell tilt from state; needs an
   external tilt tangent (neighbour classes / rigidity). REFEREE (honest resolution): a2d_controls/referee.sh + referee_score.py
   (reproject map, match_stacks, FULL 360+mirror in-plane search -- match_stacks writes the matched ref UNROTATED), outputs in
   afan/referee/<set>. Spliceosome 150: a2d150 FRC30-12 0.47 / 12-8 0.15 vs flex150 0.21 / 0.03 (sharp classes equal 0.72 vs 0.67;
   flex's 10-20 A classes really blurred) -> a2d's even/odd res is NOT just inflation on real data. View diagnostic (diag_views.py):
   within-cluster flex classes 17 deg apart; at equal view separation flex pairs differ in state labels 30-40% vs a2d 5-9%.
   FINAL_MERGED (afan/flex_merged, run_flex_env.sh): sp e/o <=10 A 30% -> 43% (m1) / 42% (m2), labels unchanged; referee pending.
(eq) 2026-10-03 RESOLUTION GAP vs abinitio2D at matched count is NOT the final-polish templates (FINAL_MERGED: referee unchanged)
   and NOT frozen labels (sp_nofreeze: e/o <=10 A 30% unchanged; labels 0.610/0.221, beyond-view 0.181 = slightly better than
   frozen). It is POPULATION BALANCING: heterogeneity mode floors every class at 0.8 x mean pop every label stage -> flex
   pops 154-478 (sharp classes mean 388) vs a2d150 2-1004 (36 sharp classes mean 697 = 54% of ptcls; 103 blurred junk-sink
   classes mean 168). Arms sp_rb0 (FLEX2D_REBALANCE=0) / sp_rb035 in afan/flex_merged. TILT TOY (afan/tilt_toy, tilt_toy.py,
   Python, one view cone, 2 states): ANY learned subspace (per-class M1 or tied M2) collapses state separation (T/3q snr .02:
   M0 0.89, M1 0.61, M2 0.81, oracle-deflate 0.89, LDA 0.89); at 5 deg cones M0 is already at the LDA ceiling; T/R12 and
   T/R6 near-inseparable even supervised at snr <= 0.02 (LDA 0.55-0.62).
(er) 2026-10-03 TILT TOY VERDICT (tilt_toy_v3.py, T vs R12 conformational, snr 0.1, 1500/state, 4 views): within-class tilt
   only matters at ~20 deg cones (real flex2D classes: particle-weighted median member angle 14 deg ~ cone radius 20 deg):
   2-class mean-only 0.559 (splits by tilt) vs LDA 0.723; 8 classes (flex2D-like over-clustering) 0.655 = oracle tilt removal
   0.658; oracle + 8 classes 0.682; NEIGHBOUR-estimated tangents (2D-only, offset 0.5-1 x cone, cos 0.79-0.96 to oracle) +8
   classes 0.650 = no gain; second-order Taylor directions do not help. At 10 deg nothing helps (M0 0.735 ~ LDA 0.77).
   Compositional (T vs 3q): tilt never matters. -> do NOT build 2D tilt modelling into flex2D; only 3D tangents add ~+0.03.
   POPULATION: at matched class pop a2d and flex cavgs agree equally with the map (sp 400-700 ptcls: 0.545 vs 0.534); a2d's
   sharp edge = 14 classes of 700-1000 ptcls; rebalance off (sp_rb0) map FRC 0.214 = unchanged, pops still 129-619.
   Merged final polish: rn_m1 0.878 = flex 0.879, r80_m1 0.766 vs 0.763 -> no gain anywhere. Paper revised accordingly.
   Cluster-merge test (referee/sp/merge): pop-weighted merge of cluster_cavgs clusters LOWERS map FRC (flex 0.225, a2d 0.47->0.36):
   clusters span ~17 deg of view, so merging blurs; removed the 'merge same-view classes' suggestion from the paper.
(es) 2026-10-03 NOT 900 matched control (7.2 h): flex2D wins every column incl. resolution: classes <=10/<=20 A 622/764 vs a2d900 323/444;
   ptcls <=5 A 0.68 vs 0.62, <=8 0.70 vs 0.65, <=20 0.86 vs 0.78. Embedding figure: make_figure_embedding.py (data/emb_<set>_table.txt from
   a2d_controls/emb_extract.py): UMAP of the 16-D global latent; DRAW IN RANDOM ORDER (fixed order made the central island look single-state);
   one view family = one cluster_cavgs cluster, sub-islands by state, flex2D classes tile them. Spliceosome clear; RNase P states intermixed.
(et) 2026-10-03 PAPER RESTRUCTURE (doc/manuscripts/flex2D_cls_split_2026_10_03, backups main_v4_prerestructure / main_v5_pre_rnfirst):
   Results ordered by heterogeneity type, user chose RNase P FIRST ("10076 is pretty low resolution"): synthetic -> evaluation criteria
   (sec:res-eval) -> RNase P (+ fig_rn_variability: within-family eigenimages localise on the arms) -> 10076 (+ fig_embedding_10076:
   agreement maps flex 55.9% vs a2d 41.3%, a2d only gets E) -> large/junk/homogeneous (NOT first + fig_not_clustered matched + fig_embedding_not:
   cross-half repro marks 201 blurred classes 23% ptcls as indistinct = label-free junk ID) -> spliceosome -> six-set summary -> cls_split -> band.
   NOT within-family variability is tilt/shift stripes (no domain) -> NOT is not a heterogeneity showcase; variability figure kept out.
(eu) 2026-10-03 METHOD-NEUTRAL EMBEDDINGS (a2d_controls/embed_aligned.py, RNase P 40k ptcls: aligned images by each method's own poses,
   64 px, 15 A, PCA16; stacks are ctf=flip -> do NOT re-flip): 15-NN state 0.334 (a2d) / 0.350 (flex), view 0.327 / 0.377, state within
   cone excess +0.070 / +0.075 -> per-particle state info is the SAME; the 2D embeddings are mostly views. flex-aligned UMAP shows view
   islands, a2d-aligned is a blob: flex2D's split children inherit the parent's in-plane frame (classes of a view co-registered), a2d
   classes have independent frames. fig_neutral_rn.png (make_figure_neutral.py).
(ev) 2026-10-03 FLOW PASS (backup main_v6_preflow.tex): abstract/intro/discussion/procedures rewritten around "judge a partition by the
   heterogeneity it resolves"; RNase P section gained fig_neutral_rn (co-registered frames point); spliceosome resolution/states figures
   dropped; sp-clustered, TRPM4, 80S montages moved to Supplementary Figures S1-S3 (after refs); cls_split table drops the mislabelled
   TRPM4 row (six data sets); large-data section ordered NOT -> TRPM4 -> PfCRT -> 80S. User: fewer non-case-study figures; every
   case-study figure must derive a point.
(ew) 2026-10-03 IUCr FORMATTING (user: "put this formatting on the paper", IUCrJ rq5017 = Van et al. 2026 polar-Fourier paper): main.tex now
   A4 two-column Times (newtxtext/newtxmath), full-width title block with editor/keywords minipage beside the abstract, 'research papers'
   running heads (fancyhdr), numbered sentence-case headings (titlesec "1." "2.1."), captions "Figure n" newline + sentence case,
   author-year natbib with a generated IUCr reference list (make_refs_iucr.py -> refs_iucr.tex: no titles, italic abbreviated journal,
   bold volume, '&', single/two/et-al ordering, 2012a/b). All figures/tables are figure*/table*; long equations split (gathered/split).
   Placeholders for authors/editor/dates (no invented metadata). Backup main_v7_preiucr.tex. The IUCr PDF itself triggers a download
   dialog -- use the article's index.html instead.
(ex) 2026-10-04 REAL HEMOGLOBIN on d1a20eff6 (503,914 ptcls, 400 cls, nstates 4, shortlist 8, POLISH_FREEZE; box150):
   seeded hb_x400v2 (NICE a2d 200 -> 400) 38.7 min, repro(indep 20-10 A) median 0.64 / 85 % > 0.5, median class res 11.5 A;
   COLD START hb_x400cold (cls_init default = one particle per template) 38 min, repro 0.77 / 99 %, median res 13.0 A (92 vs 224
   classes <= 12 A) -> the cold start now WORKS on real data (old 50k cold start gave 10-15 A blobs) but seeding gives sharper averages.
   Matched abinitio2D 400 (a2d_controls/hb400, 50 min): 393 populated, median 8.5 A, sharpest 3.8 A (shared-reference FRC, inflated).
   DECISIVE: per-particle partitions do NOT reproduce on hemoglobin: NMI flex seeded vs a2d400 0.049, vs flex cold 0.049, a2d400 vs
   seed 0.033, cold vs seed 0.031, random 400x400 baseline 0.029 (only flex vs its own seed 0.16). Spliceosome for scale: two flex
   runs 0.55, flex vs a2d 0.25. Class averages differ reproducibly (they are real views), but which particle lands in which class is
   near-arbitrary -> hemoglobin is a noise-floor test, not a heterogeneity showcase. Embedding figure: figures/fig_embedding_hb.png
   (preview vs seed: fig_embedding_hbpre.png) via make_figure_embedding_not.py <maxfam> hb|hbpre (generalised; NOT output unchanged);
   emb_extract_not.py fixed for boxes not divisible by 4. Colour-bar labels were clipped in both embedding figures (fixed, paper rebuilt).
(ey) 2026-10-04 HEMOGLOBIN MAP FRC (referee/hb: hb_T padded to 150, 3000 projections, full in-plane search; referee_score.py now
   saves per-class curves to <out>.curves.npz and crops boxes not divisible by 4). Pop-weighted mean FRC per band, >60/60-30/30-12/12-8 A:
   flex seeded 0.449/0.400/0.415/0.051, abinitio2D 400 0.439/0.357/0.435/0.059, flex cold 0.311/0.315/0.314/0.009. Shell by shell flex
   leads 60-17 A (coarse shape, 0.50-0.54 vs 0.42-0.50), abinitio2D leads 15-11 A (0.45 vs 0.30), both ~0.05 beyond 11 A. USER said
   flex hemoglobin is "visibly so much better"; the map says mixed: flex's cleaner look = better coarse shape AND less 12-15 A detail.
   Reported class resolutions (a2d 8.5 A, flex 11.5 A medians) are both beyond what the map supports. Figures: fig_hb_montage.png
   (make_figure_hb_montage.py), fig_hb_mapfrc.png (make_figure_hb_mapfrc.py). Caveat: hb_T may come from these particles (NICE).
