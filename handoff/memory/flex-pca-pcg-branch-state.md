---
name: flex-pca-pcg-branch-state
description: "2026-09-11: rec_backend=pcg for flex_pca is implemented on branch flex-pca-pcg (Mac worktree /Users/atifao/SIMPLE_pcg, batch2 worktree ~/SIMPLE_pcg built from the staging ref flex-pca-pcg-in); state maps on reconstructor_pcg + the coupled M-step on pair Gram kernels; Hans's note doc/implementation_notes/flex_pca_envelope_support.md is the document being followed; PfCRT battery tools in ~/pcg_tools on batch2"
metadata:
  type: project
---

**The document:** Hans's `doc/implementation_notes/flex_pca_envelope_support.md` (origin/master d8ae11d15, "flex with mask
note for Afan"): step 0 smpd_target working box, step 1 rec_backend=gridding|pcg at BOTH sites (state maps via sigma2/w on
reconstructor_pcg; the M-step on pair-weighted Gram kernels with the per-voxel solve as block-Jacobi), step 2 pcg_mskfile as the
solve support. User order (2026-09-11): step 1 parity with gridding first, then the mask, then step 0.

**Branch:** `flex-pca-pcg` = snapshot of the Mac's dirty consolidated flex tree (commit 8f8f5ba5a, byte-identical to batch2's
~/SIMPLE_flex the user runs) + merge of origin/master d8ae11d15 (two trivial conflicts) + my commits. Mac worktree
`/Users/atifao/SIMPLE_pcg` (the main checkout /Users/atifao/SIMPLE stays untouched, still dirty on flex-pca-resident). batch2 has no
git module but git works; ~/SIMPLE is a clone on master, `~/SIMPLE_pcg` a worktree on flex-pca-pcg. Push from the Mac to
`batch2` remote (ssh://atifao@batch2.ncifcrf.gov/home/atifao/SIMPLE) as `+flex-pca-pcg:refs/heads/flex-pca-pcg-in` (the branch
itself is checked out there), then `~/SIMPLE_pcg/build_inc.sh` (waits for make, ff-merges flex-pca-pcg-in, cmake ., make -j24
install) via nohup; build_driver.log + build.log. A fresh configure+build takes ~10 min at -j24 on the login node.
Upstream cutover gotcha: `params%l_sigma_canonical` is gone (canonical is the only sigma store) and
`ensure_sigma2_for_iteration` lost an argument.

**Code (all `rec_backend=pcg`, defaults maxits_pcg=20 rtol=1e-3 from the commander):**
- `src/main/flex/simple_flex_pca_rec3D_pcg.f90`: states per (state,half) on reconstructor_pcg; sigma2/w with a 1e-3 weight floor;
  relative lambda 1e-3; sphere support min(mskdiam/2smpd_rec, box_rec/2-6-1); workers write `flex_pca_pcgraw_part<NN>_<SS>_<e|o>.bin`
  under the part dir, master add_raw_accum in part order; cold solve max(maxits_pcg,5); delivery = window*u with the same eo-FSC
  filter policy as gridding, no zero_background/mask; `flex_pca_state_pcg.txt` table. Dispatch at the top of
  reconstruct_flex_weighted_states; `write_state` moved to util as `flex_pca_write_state`.
- `src/main/flex/simple_flex_pca_pcg.f90` (`flex_pcg_t`), RESTRUCTURED 2026-09-11 (commit after 3ded47e91) to the production
  conventions after the first design failed on real data: (1) the masked block-Jacobi start was garbage (GATE: |pcg0| 1e7x |grid|,
  corr 0) because the per-voxel divide is unbounded where rho~0 and any real-space op (the support) leaks spectrum there ->
  shell-relative floor FLEX_PCG_RHO_FLOOR_FRAC=1e-2 (Hans's RHO_FLOOR_FRAC) in a module-local floored coupled solve
  (apply_precond on nimg, packed rho of the insert + its lbound); (2) the 1x KB-deposited RHS is inconsistent (10-20%) with the exact
  Toeplitz Gram -> the RHS is now deposited at DOUBLED coordinates on the 2x lattice from the same native-subset samples
  (accumulate_rhs: data_scales x pf2 x conj(T)y, complex (ncomp, lattice) accumulators racc_e/o, packed rpk_e/o, fold_rhs,
  finalize_rhs = ifft_2N, crop, /dep_c, Nyquist-ball band limit, support); (3) CG variable = physical u, operator P Pi T Pi P
  (bare kernel, band limit on both sides so it stays symmetric), preconditioner P G_floored P (Jacobi inverse of T, NO envelope
  brackets), warm start = the gridding solution of the 1x numerators /E (budget 0 puts it straight back = gridding, bit-level
  modulo fft round trip), put_back writes E*u so the untouched tail (/E, FSC-Wiener, bp, soft mask) ships u. Part format 8/9
  (kernel block then rhs block), merge rotates rpk with rotate_rhs_packed_add (Y'=Y.R). Hooks: probe_fit_t members (l_pcg, pcg,
  kacc/kpk/racc/rpk e/o, mg_kpe/o, mg_rpe/o); fit_iter_begin, fit_batch_insert, fit_iter_reduce, pairmerge, fit_iter_finish
  (solve(Y, rho, rpk, maxits, rtol, out, tag)). Memory: rhs adds 0.54 GB/process at box 64 K=16 (kernels 2.3 GB accumulators).
  Evidence of the failure (runs 27 grid vs 28 pcg0 vs 29 pcg, old code): -2logL/N at it 2 7.7 vs 46.8 vs 34.6, eigenvalues 98% off,
  20 CG iterations never converged (update 0.3-0.6).
- Self-test (production/tests/simple_test_flex_pcg) after the restructure: (A) operator vs exact Gram 1e-7 on-lattice, 1-2% off-
  lattice; (B) 2x rhs deposit vs exact adjoint <1%; (C) 48-slice synthetic solve: 12-13 its to resid 5e-4, 0.4% error in the inner
  half band (the outer ball is undetermined by 48 slices: judge recovery in-band, never full-band).
- CONDITIONING (synthetic sweep + real data, 2026-09-11 ~06:30): the unregularised masked LS is only PSD; noise-free it converges
  in 8-13 its, with noise the CG chases the null space (err x78) and stalls; a relative Tikhonov lam = lam_rel x low-band mean
  diagonal rho (Hans's PCG_LAMBDA analogue; module default 1e-3, env SIMPLE_COV_PCG_LAMBDA) fixes it: on PfCRT lam_rel=1e-2 takes
  the M-step residual 1.85 -> 3.5e-2 in 20 its (1.04 without). The cross-FSC ridge (SIMPLE_COV_XFSC_REG=1|2) alone did not change
  the CG at all (too small where it stalls) but improves gridding's -2logL (gridx). Vanishing curvature now returns the iterate
  (no THROW). Arms 38 gridx2 / 39 pcgl2 (1e-2) / 40 pcgl1 (1e-1) running; step 2 (pcg_mskfile) implemented in commit after
  a3c33b7a3: flex_pcg_support_volume (read_and_crop, min/max/out-of-range log, clamp, floor 0.05, occupancy native vs target),
  flex_env_init/flex_window_apply/flex_pcg_install_window replace every sphere site (mstep tails, dfl_basis, basis init, probe
  load, merge gate, merged solve) and set_mask_volume at box_rec in the states path; launcher arms pcgm/stpcgm use
  9_refine3D_auto/automask3D_state01.mrc.
- ARBITER RESULT (-2logL/N at probe it 5, halves A/B; lower is better): grid -0.06/0.52; gridx (REG=1) -1.30/-0.41; gridx2
  (REG=2) -0.71/0.16; pcg0 (gate) = grid to 7 digits; pcg (lam 0, 20 its unconverged) -0.43/0.62; pcgx -0.05/0.98; pcgl2
  (lam 1e-2, converged resid 3e-2) 1.81/2.75; pcgl1 (lam 1e-1, resid 1e-2) 0.95/2.00. So the CONVERGED Tikhonov-LS M-step is
  ~2 units worse than gridding+FSC-Wiener while the unconverged one is on par: the exact coupled LS fixed point is not the
  better estimator with these priors; the FSC-shaped prior inside the operator (P_tau, the note's "second version") is the next
  thing, and budget 0 (gate) is the safe default. State maps: populated state corr 0.997 pcg vs grid.
- STEP 2 on PfCRT (runs 41 pcgm = paired pcg + pcg_mskfile + lam 1e-2 + REG=2, 42 stpcgm = states with the envelope): the
  automask3D_state01.mrc occupies 2.4% of the box = 33% of the mskdiam sphere; Fourier crop to box 64 rings (min -0.08, max 1.06,
  50% of voxels outside [0,1] before clamping, all solvent ripple) and the 0.05 floor restores the occupancy (32.5% of sphere);
  no ringing at box 300. The envelope-windowed INITIAL basis starts EM at -2logL/N 22.9 vs -3.1 for every sphere arm: PfCRT's
  micelle lies outside the protein envelope and the constrained basis cannot model it (the note's 5 anticipated exactly this
  test). The registry bound pcg_mskfile as 'B' (binary) and rejected MRC masks (fixed 883fe7ea5, 'M','S'); upstream refine3D
  would hit the same. Status note: doc/implementation_notes/flex_pca_envelope_support_status.md (untracked until committed).
- Envelope M-step trajectory (run 48, lambda 1e-2, REG=2): -2logL/N 22.9/14.2/9.1/6.8/4.5 vs the sphere at the same settings
  (pcgl2: 9.8/7.5/4.7/1.8): the protein automask is the wrong support for PfCRT (micelle outside); try nu_envmask3D_state01.mrc
  or a dilated automask next. A class(image) dummy in flex_window_apply SEGFAULTED at -O3 (gfortran 16) on the first
  uninstrumented call site (M-step tail, then the merge tail once the first was instrumented): typed entry points
  flex_window_apply (image) / flex_window_apply_rec (reconstructor) replaced it (4e25c4d98); never pass array elements of
  type(image) to a class(image) dummy in flex at -O3. CONFIRMED: run 49 (job 53187278, snapshot 4e25c4d98) passed both the
  M-step and the merge tails (identical -2logL trajectory) and went on to the state stage. HEAD of flex-pca-pcg at hand-off:
  the doc commit after 4e25c4d98 (status note doc/implementation_notes/flex_pca_envelope_support_status.md); branch pushed to
  batch2 (flex-pca-pcg-in, ff-merged into ~/SIMPLE_pcg). Merge into master is the user's call.
- Real-data status (2026-09-11 ~05:00): restructured pcg0 (run 32) tracks gridding within ~10% in -2logL/N (8.73 vs 7.72 at it 2)
  -- the residual gap was the hard support P applied to the budget-0 solution: the tail's half-set FSC is taken on UNMASKED
  halves and P inflates it (weaker Wiener shrinkage). Fixed 4b2662009: budget 0 ships the gridding solution untouched; with
  iterations the solution lives on P by construction, so the pcg arm's FSC is the masked one (a systematic difference to note
  against the -2logL arbiter). CG on real data converges SLOWLY: init rel resid 1.8 (warm start far from LS), 1.04 after 20 its;
  diagnostics corr(b,Bx0)/scale + SIMPLE_COV_PCG_VERBOSE=1 per-iteration log added to find out whether it is conditioning or a
  convention mismatch. State maps (stgrid 25 vs stpcg 30, Hans's reconstructor at box 300, 20 its, resid 3.5e-2): populated
  states corr 0.97/0.98 with identical eo-FSC (8.2/4.2 A); the two ~45-particle states are noise on PCG (corr 0.12, FSC 310/44 A
  vs 21/8 A gridding) -> PCG needs the ML prior / lambda for near-empty states; 3 min per (state,half) solve vs 4 min total gridding.
- Self-test `simple_test_flex_pcg` (production/tests): PASSES since 933bdf32d against an EXACT nonuniform-DFT Gram (kernel within
  0.6-1.2% at box 32/64; the direct KB-interpolation Gram is itself 18-22% off, so never use it as the reference). Two bugs found
  by numeric probes (delta/const->fft, kpk/tker/khat spike location): (1) alloc_accum allocated the pair accumulator with LOGICAL
  lower bounds while the scatter (assumed-shape dummy, 1-based) and fold_accum (allocatable dummy, keeps bounds) disagreed ->
  fold read past the end (mass landed at the Nyquist corner, alternating-sign output); now 1-based everywhere. (2) fold_accum
  stopped at h=lims3(1,2), missing the 2x lattice's Nyquist plane (its column is the wrapped -Nyquist one). kscale=padf**3=8
  confirmed numerically (khat(0)=16 for two mates at the origin).
- Worker sigma failure after the upstream canonical cutover ("nonpositive or nonfinite sigma2", model.f90 load_and_validate_sigma):
  esig%new allocates sigma2_noise on params%fromp:top but a distributed worker's pinds are the master's pindfile partition, so rows
  outside its fromp/top were never filled. Fixed 2026-09-11 by widening fromp/top to 1..noris around load_sigma2_groups (the unit
  fallback already did that). Both backends failed identically before the fix; unrelated to PCG.

**PfCRT on batch2:** /mnt/beegfs/elmlund/pfcrt, project 9_refine3D_auto/pfcrt.simple (box 300, 1.035 A, 16905 ptcls, mskdiam 160,
automask3D_state01.mrc present for step 2), 14_flex_pca = the user's paired baseline (embedding for states-only resumes). SLURM
csbdevel; `script=yes tag=<t>` writes `flex_pca_script_<t>` / `FLEX_PCA_OUTPUT_<t>` with NO env exports (add SIMPLE_PATH/PATH and
the recipe after the `cd`, mem 128000). `~/pcg_tools/launch_pfcrt.sh [states|full|all]` generates+chains (after:+1) stgrid/stpcg
(infile resumes) and grid/pcg0/pcg (paired recipe); `~/pcg_tools/pcg_cmp.py A B` compares eigenvalues, |z|, state maps, eo-FSC.
Launcher runs jobs from `~/SIMPLE_pcg/snap_<commit>` (bin+scripts copy) so rebuilds never replace binaries under running
processes; `build_inc.sh` removes build.log first so a waiter can poll for "make rc=". Job log ~/pcg_tools/pfcrt_jobs.txt.
Runs: 25_flex_pca = stgrid (53186689, done), 26 = stpcg first try (state workers 2-4 failed: rounds%nparts() is 1 on a worker ->
fixed b19175dac, params%nparts), 27 = grid (53186699, done, the control), 28/29 = pcg0/pcg on the OLD M-step (cancelled, keep as
the failure record), 30 = stpcg (53186790, snapshot b19175dac; state solves 20 its at box 300 ~3 min each, resid 3.5e-2 at
maxits), 31 = pcg0 GATE diagnostic (cancelled). Earlier 25-28 from the failed sigma attempt were deleted before relaunch.
Launcher takes arm names (`launch_pfcrt.sh stpcg`, `pcg0 pcg`); 2.6 GB raw part files per (state,half,part) at box 300.
See [[flex-pca-real-space-priors-proposal]], [[flex-pca-conventions-branch-state]].
