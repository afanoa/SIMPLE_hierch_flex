---
name: flex-cls-split-state-2026-09-30
description: "cls_split pca_mode=flex (branch flex-cls-split, merged onto master after the flex refactor, on origin): suite 36/36; in-plane pass, FRC Wiener restoration, repro column; beats diffusion maps on synthetic and label-free band readouts on 10180 and the converged 10076"
metadata:
  node_type: memory
  type: project
  originSessionId: 03e8b343-4213-4208-890c-56346e7de4b2
  modified: 2026-09-30T21:01:59.142Z
---

Branch `flex-cls-split` (Mac ~/SIMPLE_flexcls, batch2 ~/SIMPLE_flexcls with snapshot binaries in
~/SIMPLE_flexcls/snap/bin used by all jobs). Plan + results: doc/implementation_notes/flex_cls_split_plan_2026_09_30.md.
Runs: /mnt/beegfs/elmlund/afan/flexcls/{t4_mix5050,t4_mix5,10076,10180,10028_ref}.

What it is: per 2D class, CTF+sigma2-weighted PPCA on the 30 A band (rank 8, 3 ALS restarts,
per-class seed), pose tangents + class mean as nuisance columns, divisive bisection by Ashman's D
into exactly ncls subclasses, CROSS-FITTED embedding (5 folds, held-out latents mapped through the bases: z = G^-1 U_full^H W U_fold z_fold; Procrustes on latents smeared rank-4 columns) with a
posterior-precision kernel for the weights, CTF-corrected weighted sub-averages. Distributed via the
cls_split master/worker parts; 46k particles at box 320 in ~80 s on 80 cores; 327k in 5 min.

Results (subclass majority purity / NMI):
- synthetic 50/50 (EMD-8440 vs 8445, abinitio2D 60 cls): flex 3-way 0.85/0.23 (143 of 180 pure),
  diffusion maps 0.56/0.01. Five-state ladder: flex 0.62/0.18 vs dmap 0.42/0.01.
- 10076 (published labels, 200 cls): parents 0.359; every flex variant 0.362-0.368; dmap 0.362.
- 10180 46k selection (34 cls, 5-way, vs flex_pca 3D states): parents 0.615/0.253; flex hard
  0.664/0.255, cross-fitted kernel 0.636/0.229; dmap 0.619/0.227.

Evening additions: placement latent = [z, log cross-fitted residual] (residual competes only on its own
axis in the bisection) -> junk test 1.0/1.0/1.0 at 15 % and 30 % junk; unit suite 27/27 in 4 s is the
development gate (user: "proper synthetic test suites"). 10180 calibration: cross-fit scatter /
posterior width 1.13-1.76 per component (1.0 = noise), in-sample 2.0-2.4.

2026-10-01: continuum ladder test added (32/32). Two levers built, measured, DELETED at the user's
request (plan note sections 11-12 keep the numbers): (a) full-band likelihood EM of the weights
(free mixture collapses on a continuum; mass-balanced version only shuffles, 10180 0.646 -> 0.620);
(b) same-view pooling of parents (9 mutual-NN pairs, composition verified; +0.034 vs +0.031
unpooled). Standing best, 46k 3-way vs 3D states: flex 0.646/0.247 (smallest subclass 46), diffusion
maps 0.616/0.237 (smallest 6). The user's rule: delete code that is measured not useful.

2026-10-01 later: label-free readout `repro` (cross-half reproducibility of subclass differences,
even/odd stacks written) is the validation that matters: diffusion maps reproduce only >60 A
(background), flex reproduces at 60-30 A on 10180 (0.74) and 10076 (0.45). It also caught a flex
defect (residual axis followed defocus on 10076, 0.97) fixed by detrending the residual against CTF
power (0.17). Scorers on batch2: flexcls/repro_eo.py, repro_bands.py, split_nature.py, defocus_ratio.py.

2026-10-01 evening: branch merged onto origin/master (be67a33a4, the flex architectural refactor)
with parity (36/36, 10180 0.648/0.249, bands identical) and pushed to origin as flex-cls-split.
Reuse from the refactored flex: gmm_state_weights, kmeans_latent_targets, spd_inv_dp (posterior
service); the rest of the refactored services take a build/reconstructor and are 3D-only. Added
since the morning: cross-half `repro` readout + even/odd stacks, residual detrended on CTF power,
FRC-based Wiener restoration (per half, then merged), Friedel mates of the h=0 column on write-back
(stripe fix), one in-plane pass per member vs the other half's sub-average (band to FRC 0.143),
SIMPLE_FLEXCLS_LABELS ablation (model vs restoration), SIMPLE_FLEXCLS_NOALIGN. Write-up:
doc/implementation_notes/completed/flex_cls_split_2026_10_01.md.

2026-10-01 night: the diffusion-map/kPCA cls_split path was DELETED at the user's request ("useless");
cls_split = flex only, ptcl2D only, ncls>=2 required, no pca_mode (git history before the removal
keeps it). Strategy 2106 -> ~1400 lines. cluster2D's in-plane mode (simple_private_exec
prg=cluster2D_distr refine=inpl) is callable on a split project but needs a canonical sigma2 state
path (or objfun=cc); experiment in flexcls/10180/inpl.

2026-10-02 early: the brute-force in-plane pass is gone; cls_split = split + one cluster2D round
(refine=inpl_local: 10-degree window, fixed-angle shift search object; startit=3; inpl_cont=no) +
restoration from the refined alignment, in one call (~4 min for 46k). Bands: 10180 0.88/0.62,
10076 0.81/0.79 (dmap 0.27/0.07 and 0.34/0.19). Three traps that cost a night: (a) cluster2D's
factory treats nparts=1 as distributed (worker env/logs lost) -> omit nparts for in-memory tests;
(b) the angle-optimising pftc_shsrch_grad re-selects the rotation whatever irot is passed -> use the
fixed-angle object (grad_shsrch_obj2) for a bounded refinement; (c) oris%get_label_inds reads the
class field whatever the label. Build: git hash moved to simple_gitinfo (SimpleGitHash.h); a commit
now recompiles one object, not 603.

Lessons that cost time: (1) unit noise weights collapse the basis (use sigma2 or the class residual
spectrum); (2) the largest-variance latent direction is a pose residual and the most bimodal one is
junk, hence divisive bisection; (3) the IN-SAMPLE posterior precision is overconfident (a basis fitted
on a member's own noise resolves that noise): posterior kernels collapse to empty tiles; cross-fitting
fixes the collapse; (4) pooling one basis across views is worse than per-class (60 vs 143 pure);
(5) 3_selection is a 34-class selection; 2_abinitio2D's compenv is qsys_name=local; (6) distributed
cls_split worker stdout is deleted by qsys_cleanup, per-class prints need a shared-memory class=N run.

User's verdict 2026-09-30: heterogeneity "more visible on the flex but not good enough, long way to
go". The limit on real data is per-particle SNR of the state difference within a view class, not the
model form; the cross-fit calibration (latent scatter / posterior width per component) is the per-class
readout of that, not yet inspected on real classes. Related: [[keep-test-loops-minimal]],
[[use-the-full-cpu-cap-for-tests]], [[flex2d-findings-2026-09-21]].

**2026-10-02 update.** Branch flex-cls-split at a04218870 (on master be67a33a4), local + batch2 remote only.
Core footprint reduced to the two optional args on transform_ptcls: the 10-degree `inpl_local` window and the
cluster2D sigma2 bootstrap were measured and removed (stock `refine=inpl` equals the window on the same split,
0.870/0.631 vs 0.877/0.623; gate 0.874/0.632; suite 34/34). The commander runs calc_pspec itself when the
project's sigma2 state does not resolve (projects copied out of their abinitio2D dir register it by a bare name).
Seven real sets run (dirs `flexcls/{trpm4,flipqr,10028,not,e13553}/x3`, generic launcher `flexcls/runds.sh
<dir> <proj> <args>`, montage `montds.sh`, scorers `split_amp.py`/`split_amp2.py` = effect size + balance).
Effect size ranks sets the way montages read (ribosome 2.87, spliceosome 1.06, not 0.57, 13553 0.33);
reproducibility alone does not. 13553 (one 550 GB stack) is 10x slow in split/restore from per-class random
member reads; cluster2D streams in index order. Fix = read each worker's range once in index order (not done).
10028 is slow in every phase incl. stock cluster2D (its data path). Snap binaries = a04218870's code.

**2026-10-05.** Renamed cls_split -> cls_expansion (program, commander, strategy, model `simple_flex_cls_expansion`, suite key
`flex_cls_expansion`, outputs `cls_expansion_*`); branch flex-cls-split at e8d1561a8. Model audited to conventions (explicit
`use only`, named constants, no umbrella import). Minimal tester = 5 tests / 15 checks (local only, not pushed). One note:
`doc/implementation_notes/completed/flex_cls_expansion_2026_10_01.md` with the plan as appendix. Gate 0.872/0.631, suite 15/15.
TRAP: an incremental batch2 build after reordering optional dummies gave a gate of 0.742/0.359 with identical labels (strategy
object compiled against the old .mod); a full rebuild (cmake rerun) restored 0.87/0.63. After interface changes, rebuild fully.
batch2 scripts renamed (run_arm/runds prg=cls_expansion, unit.sh suite, scorers accept both prefixes). gitinfo build change
parked (not for the push); push plan = squash onto origin/master (38 commits ahead, cluster2D->refine2D rename) minus tester/gitinfo.

**PUSHED 2026-10-05:** origin/master 08a46838b -> 41ec883e4 "Flex_pca version of cls expansion, trashing diffusion map
version" (one squashed commit, 15 files, +3220/-1668; tester and gitinfo build change NOT pushed, stay on local branch
flex-cls-split). Rebase resolutions: master's CMake flag block and deletions (matcher source reader, stream p03,
denoise_project commander) win; UI lives in simple_ui_refine2D.f90; commander uses commander_refine2D / prg=refine2D,
inpl_cont dropped (internal flag on master); strategy uses qsys_declare_part_finished. Verified by master's compile_clean.sh
and compile_debug.sh on a batch2 compute node (worktrees ~/SIMPLE_push_rel, ~/SIMPLE_push_dbg): fast gate 14/14 both,
28-29 s of a 30 s budget. Local branch cls-expansion-push = the pushed commit.
