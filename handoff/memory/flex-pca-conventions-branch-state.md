---
name: flex-pca-conventions-branch-state
description: "2026-09-07/08: Steps A-E, E2, F4, G1, G2, G4, G5 gated and pushed on box branch flex-pca-conventions (tip a21fac77e); the box working tree additionally carries G6 + a ~3 800-line dead-code cut (34 env switches) that is UNCOMMITTED and UNBUILT (no builds allowed while elmlundho runs); F1-F3, G3, G7 remain; the four-arm 400-particle gate is ~2 min"
metadata:
  node_type: memory
  type: project
  originSessionId: 7e382bce-970c-4cdb-bc32-4a1a1c1801a9
  modified: 2026-09-07T22:30:00.000Z
---

Box `~/SIMPLE-h0901-int2`, branch `flex-pca-conventions` from `auto4-rhofix`
(5c2442dcd = the validated recipe + delivery fixes). Local repo has remote `box`
(`atifao@fwl-c142339.ncifcrf.gov:SIMPLE-h0901-int2`); `git fetch box
flex-pca-conventions`. Commits: a741d8f12 (A-C), 8f6c6a1fb (D), 351647f3f (E), 553e3ad08 (F4), 3f0d01d80 (G1), 89fb652c0 (G4), G5 (em_iter split);
`doc/policies/flex_pca_policy.md` is the execution contract; the plan doc's
status header lists the deviations.

Shape now: commander -> `create_flex_pca_strategy(cline)` -> shmem/master/worker
strategies (`strategies/parallelization/simple_flex_pca_strategy.f90`); the domain
receives a `flex_pca_rounds` (leaf `simple_flex_pca_rounds.f90`) and never queries a
module-level role; workers get `pindfile=` lists and `which_iter/maxits/nfits/stage`
in job_descr; `PCA_STAGE_POLISH`=7 is the joint fit; `simple_flex_pca_distr.f90` is
gone. Worker entry: `run_flex_pca_worker` -> `probe_worker_pass` /
`embed_worker_pass` / states accumulation.

**Gate harness** (box): `~/gate_run.sh <outdir>` (four arms s1/s4/p1/p4 on 10028 rows
1001-1400, n_probe_iters=2, ~2 min total), `~/gate_cmp.py <base> <new>` (eigenvalues,
coordinates, embedding magnitudes), `~/gate_step.sh <tag>` = regenerate args + cmake +
make install + gate + compare vs `gate_base`. Never `pgrep -f` a pattern that your own
ssh command line contains (it killed my shell twice); anchor with `^/bin/bash ...`.

**How to apply:** continue with F1-F3 (probe_subspace_iteration onto the paired phases
with fits(1:1) + fit_estep_device_pass; master_execute owning the loop; embed/states
accumulate-reduce-finalize splits), the parts dissolution, G2 (projected_latent_model
survivors into latent_ops), G6 (em_pose split), G3 (GPU dead code, CUDA build) on this
branch, one gated commit per step. Over the cap now: em_iter 2585, gpu 2505, model 1870,
em_pose 1620. Merge target is the user's call. See
[[flex-pca-polish-destroys-island]] for the recipe the branch carries.

**Update 2026-09-07 (evening):** remote tip is now cbcd762bd (Step E2, parts module dissolved); the box checkout
`~/SIMPLE-h0901-int2` is at a21fac77e (Step G2, projected_latent_model merged into latent_ops) with `em_pose.f90` dirty and a new
`em_polar.f90` (G6 in progress). `box/h0901-int2` has 5 commits NOT on this branch (SIMPLE_COV_FOCUS_MSK d9eb70d34, Sinkhorn states,
config-of-record default). origin/master (124b277f4) has diverged 62 commits since the common base dc2eb9a74; a merge-tree dry run
conflicts in six files (sigma2 canonical refactor vs the moved flex files). Branch plan for the priors work is in
[[flex-pca-real-space-priors-proposal]].

**Update 2026-09-08 (dead-code cut, no builds):** the box working tree on top of a21fac77e now
holds G6 (`em_pose` deleted, `em_polar` 400 lines) PLUS the approved "delete outright" cut:
EM_POLAR former + parity diagnostic + DEFLATE_POSE + pose perturbation; ONLINE windows/blend
(`fit_window_walk`, `full_pass`, ob_lo/ob_hi -> 1..npp); XFSC march/stop/scaffold/anneal
(`xfsc_ctx_t` keeps ridge + paired writer, record layout unchanged, `XFSC_CRIT` constant);
GMM deconv/orient (+ `views` chain); ARCHIVE_BASIS (it000 stamps stay, gated on
MERGE_PRUNE/PAIRED_MERGE because the merge reads them), DUMP_ACC, ZHALF, *_STRIDE, OSTRAT,
RKW, AWEIGHT, DSIGNAL, CGSOLVE, EMBED_RECAL (+recalibrate_gamma), ZSHRINK (now unconditional),
ZWHITEN (+whiten_latent_frame), MASK_IMAGES/CONTRAST/SPLIT_CHECKER overrides. git diff: 19
files, +202/-3986; flex dir 20 239 lines. Scripts on the box: ~/cutParity.py, cutEmPolar_*.py,
cutOnline.py, cutXfsc.py, cutGmm.py, cutInstr.py, cutStale.py (all idempotent-guarded;
fblock.py is the block remover). NOT syntax-checked, NOT gated: the first action when the box
is free is `~/gate_step.sh deadcode` vs `gate_base`, then commit, then apply F1a/cache patches.
Production driver still exports SIMPLE_COV_ONLINE=0 (now a no-op; harmless).

**Update 2026-09-08 (squash + rebase):** everything (47 branch commits + G6 + the dead-code cut +
README v2) is now ONE commit 5e04aefcc on box branch `flex-pca-squash`, parent origin/master
124b277f4 (46 files, +11 917/-9 306); fetched locally as `box/flex-pca-squash`. Conflict
resolution: our side for all flex hunks; master's canonical `simple_sigma2_files.f90` taken with
our group-count dispatch ported into its legacy path; master's `ensure_canonical_sigma_state`
called in the commander after `strategy%initialize` on the master process only;
`load_sigma2_groups` takes `build%spproj`. Backups: tag `flex-pca-conventions-pre-squash` and
branch `flex-pca-conventions` (tip 6b4eae16a WIP commit) untouched. Not built, not gated: the
first action when the box is free is `~/gate_step.sh` on `flex-pca-squash` (note: gate_base was
made on the pre-master tree; a bit-identical result vs gate_base still needs the master-side
changes to be numerically neutral).

**Update 2026-09-08 (afternoon):** the dead-code cut had removed a LIVE else-branch (the polar/Cartesian
E-step former call before `fit_estep_solve_stats` on the single-fit CPU path) together with the parity
block: every particle then subtracted an unallocated mean plane -> silent OOB in release -> joint-fit
worker segfaults. Found by bisection (pre-cut passes, cut-only fails) + `-fcheck=all` Debug tree; fixed
(former dispatch restored) and gated bit-identical: squash tip 44ed91143 (branch flex-pca-squash), which
also carries: legacy sigma default (canonical bootstrap fails on 10028: no eo split), master thread
restore after the sigma hook, master budget = nparts*nthr (no 64 cap), esig built via `esig%new` in the
unit-sigma fallback. Branch flex-pca-dev (9b0446e9d, worktree ~/wt_dev) adds: `umap=yes|no` typed
parameter (default no; production driver passes umap=yes), non-worker label delivery, merge-walk
weight lookup, numlen to qenv%new -- gated bit-identical. Production ribosome paired stage now 1312 s
(was 2335). **cache=yes is still NOT consumed** (`l_pcache=.false.` in em_iter); enabling it on the
400-particle gate gives 43% eigenvalue differences while per-particle residual energies agree to 0.1%,
so it needs the 20k harness to judge, not the gate. Debug-tree caveat: `-fcheck=all` reports a false
sigma-table bounds error in the mean-stage prep unless an unrelated print is present.

**Update 2026-09-08 (evening):** work continues on box branch `flex-pca-resident` (worktree `~/wt_res`,
build `~/wt_res/build`; local worktree in this session's scratchpad `wt_resident`; pushed via the
scratch ref `flex-pca-resident-in` because the branch is checked out on the box). Commits on top of
flex-pca-dev 7ffebe43d: c9a9a3b7c cache reads in the E-step (queue5 tree), d284f913a/7c07e9e0a
resident planes, 1e46ddb7e taper_width (shared image edit), d0feeaec8 mean+embed on the cache read
path, 2dcb67ba3 residency only on the cropped grid, 0387b433c continuum kernels. Gated bit-identical
through 2dcb67ba3 (raw arms); the taper-fixed cache arm and the v7_kernel states rerun were queued
(queue9b/queue10) when this was written. See [[flex-pca-cache-taper-mismatch]],
[[flex-pca-continuum-kernel-fix]], [[flex-pca-speed-levers-2026-09-08]].
