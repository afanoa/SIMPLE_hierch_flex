---
name: flex-pcg-pool-memory-fix-2026-10-16
description: "REVERTED - fix was downscaling instead (refine3D_states_cavgs flex at the lpstart box). flex_pca PCG (flex_pcg_t) memory at box 420: per-thread 840^3 pool sized to ALL threads + kt copy of the 55-pair khat -> ~750 GB for refine3D_states_cavgs' paired fits; fix = pool width min(nthr,ncomp) + chunked fold (no kt); A/B pcgmem_ab"
metadata:
  node_type: memory
  type: project
  originSessionId: 950a20f7-af3a-4db4-93f0-bde563ac38a8
  modified: 2026-10-16T15:20:39.676Z
---

2026-10-16, NOT refine3D_states_cavgs (7955 state-3 averages, box 420): its flex_pca initialisation runs in the master,
shared memory, box_crop=420, rec_backend=pcg. At nthr=16 it OOMed a 64 GB job; at nthr=32 on cn001 (500G request) RSS
reached 735-775 GB (not killed). User: "is there a memory leak" -> no: `simple_flex_pca_pcg.f90`
- ensure_pool sized the per-thread pool to omp_get_max_threads(): per slot wpool + tp_emb (boxpd^3 real) + tp_acc
  (cdim complex) = ~7.5 GB at boxpd=840 -> ~250 GB per solver at 32 threads;
- khat = npairs(55 at ncomp 10) x cdim real = 65 GB, plus finalize's staging copy kt (another 65 GB);
- op_cq/op_hq (ncomp x cdim complex) ~47 GB per solver; two solvers live (paired fits A/B).
Production PCG (simple_reconstructor_pcg) uses one wimg with threaded FFTW, one Khat per state/half, colouring
instead of per-thread grids, master cap 32 -> a few GB per state.
Fix (uncommitted, ~/SIMPLE_m1015 + scratch worktree, snap/b57538b7b_vol1_pcgmem): pool width min(nthr,ncomp); every
pooled q-loop num_threads(nthr_pool); finalize folds pairs in chunks of the pool width (slot = chunk position) and
scatters each chunk into khat, no kt. Expected ~190 GB per solver (~430 GB run) instead of ~375; fold up to ~3x slower.
Numerics by construction identical (same per-pair/per-component ops, single-threaded FFTs, same real()*KSCALE).
A/B job 53905436 in /mnt/beegfs/elmlund/afan/toyexp/pcgmem_ab (old1/old2/new, flex_pca box_crop=128 nthr=16 on the
15,910-row work project copy; compare.py). Result: all 84 state maps byte-identical old1/old2/new; embedding +
21 txt files differ but old1 vs old2 differ identically (max rel 5.0e-3) = flex's own run-to-run nondeterminism at
nthr=16, NOT the change. RSS 33.3 -> 32.8 GB at box 128; fold total 205-212 -> 239 s.
Then user: "do the production pattern": l_tfft mode (+473/-53 lines): above FLEX_PCG_POOL_BUDGET_GB=16 (pool width x
slot bytes) no pool; transforms one at a time on the shared wimg/nimg (threaded FFTW plans, nthr_glob) with OpenMP
array passes (cmat_view via get_rmat_ptr_padded + c_f_pointer; fill_rmat/fill_cmat/masked_copy/bandlimit_tfft;
finalize_tfft stages nb=ncomp kernels then scatters into khat). SIMPLE_FLEX_PCG_POOL_GB overrides (0 forces it).
Built snap/b57538b7b_vol1_tfft. Tests: 53905813 (box128 forced, vs old1/old2), 53905814 (box420 default, n_probe_iters=1
rec_states=no, nthr=32, 450G) in pcgmem_ab/{tfft,box420_tfft}. box420 tfft: fold 75-80 s (old 88-147), 274 GB at
fit A, 381 GB peak with BOTH fits' solvers alive (A keeps khat+op arrays after its M-step). box128 forced tfft: first
PCG MSTEP line identical to old to all printed digits, but fold 284 s vs 8-12 s (small-lattice threaded-FFT overhead;
small boxes use the pool by default). User: "the master should fit in 256 GB" -> tfft2 (+513/-53 total): solve()
releases khat/op_cq/op_hq/pc_cq at its end in tfft mode (rebuilt by the next finalize anyway) + in-place coupled
multiply into op_cq (no op_hq). Built snap/b57538b7b_vol1_tfft2; tests 53905940 (box420, rec_states=yes, NORMAL node
--mem=240G) and 53905939 (box128 forced numerics) in pcgmem_ab/{box420_tfft2,tfft2}.
Related: [[box-oom-from-paired-pcg-solve]], [[flex-pca-ab-testing-traps]], [[flex-pca-exactness-verification]].
User review 2026-10-16: comments cut to 2 lines ([[max-two-comment-lines-per-change]]); SIMPLE_FLEX_PCG_POOL_GB env switch
REMOVED (user: "the environment switch should be removed"); direct cmat_view buffer access + dual pool/threaded paths
accepted ("the other two are fine"); FLEX_PCG_POOL_BUDGET_GB=16 is my arbitrary switch-over (pool <= 16 GB), crossover
untested. Current build snap/b57538b7b_vol1_tfft3 (+484/-57 over b57538b7b, 4 files).
OUTCOME 2026-10-16: FFTW is linked as libfftw3f_threads (pthreads; cmake/Dependencies.cmake prefers it over _omp);
threaded flex PCG at box 128 = 284 s fold vs 10 s pool; LD_PRELOAD libfftw3f_omp -> 7.7 s (cause confirmed);
OMP_WAIT_POLICY=passive no help. User declined the FFTW backend change and asked for ONE memory-safe path, then
suggested downscaling: ALL flex PCG changes REVERTED. Fix instead: refine3D_states_cavgs's policy l_flex_native_box
(box_crop=420, hael f83f12687) -> l_flex_lpstart_box: flex gets lpstages_setlims(box,1,smpd,lpstart,lpstart) box_crop
(box 162 @ 3.37 A for NOT, the march's first-stage factor; march itself runs at the lpstop box ~272). flex state maps
are delivered by reconstruct3D at native sampling regardless (project_gateway), so the controller check holds.
Build snap/b57538b7b_vol1_lpbox (vol1 + lpbox only, 4 files +36/-9).
Then (user): solve3D_cavgs maps stay vol_cavg (old convention: cavg maps are vol_cavg since 2020; particle chain never
relabels maps - solve3D state=N / selection + refine3D_auto vol1). vol1 option and selection carry_state_maps DROPPED;
refine3D_states_cavgs gets state=N like solve3D state=N: cline state taken+deleted early, restrict_to_state
deactivates averages (cls2D) and particles of other cls3D states (labels kept, so map N matches), consensus = vol N
else vol_cavg N; no state= -> master vol(1). Patch (worktree, not built): 4 files +58/-10 incl. the lpstart-box change.
Chain: 4_solve3D_cavgs -> refine3D_states_cavgs state=3 nstates=3 (no selection step).
FINAL design (user 2026-10-16): selection route + explicit vol1 ("to avoid any confusion"); state=N and
carry_state_maps dropped. Patch = refine3D_states_cavgs vol1 (+UI, how2) + flex at the lpstart box (controller) +
selection state= for cls2D/cls3D rows (l_require_updated only for ptcl segments; class rows have no updatecnt, so
selection oritype=cls3D state=N used to deselect everything silently). 5 files +37/-10; the running job 53906124
(snap b57538b7b_vol1_lpbox) is this minus the selection fix. Chain: 4_solve3D_cavgs -> selection oritype=cls3D state=3
-> refine3D_states_cavgs nstates=3 vol1=4_solve3D_cavgs/rec_final_state03.mrc.
