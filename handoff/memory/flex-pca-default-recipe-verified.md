---
name: flex-pca-default-recipe-verified
description: "2026-09-12: the flex_pca cleanup is DONE and verified on 10028 -- paired+merge is the workflow, 70 of 84 env switches deleted (~1900 lines), rec_states is a CLI parameter, min_neff is now an occupancy floor; a BARE command line with no SIMPLE_COV_* at all reproduces the validated recipe (0.803 vs its 0.804 self score)"
metadata:
  type: project
  originSessionId: 7e382bce-970c-4cdb-bc32-4a1a1c1801a9
  modified: 2026-09-12T17:57:47.902Z
---

**Run:** `flex_h0901/default_b64_lp16` on the box (10028, 105,247 particles, cryoSPARC poses).
Two invocations, since the paired driver still returns after the deconvolution:
1. fit+embed, env = `PAIRED=1 PAIRED_MERGE=1 PROBE_ONLY=1 MERGE_POLISH=1 XFSC_REG=1 POLAR_ESTEP=1
   EM_DEFLATE=4 DEFLATE_BG=1 EM_MIX=16 PROBE_CONV=999999 SOLVENT_MSK=0 GPU=0 STATE_REC=0`,
   cli = `neigs=20 n_probe_iters=6 box_crop=64 box_rec=360 lp=16 mskdiam=320 npreimages=8
   state_axis=0 min_neff=2000 nparts=10 nthr=8 cache=yes umap=yes`. 1024 s with a warm plane cache.
2. states: same env with `STATE_REC=1`, `infile=flex_pca_embedding_deconv.bin`, on a FRESH copy of
   the original project (a delivered project carries state labels and is refused). 870 s.
NOT set (and provably inert): MERGE_PRUNE (the CV rank cut only ever logged "nothing pruned"),
RELPRIOR (the paired path forces the plain prior), the solvent/focus masks, SIGMA_RECALC, COMPOSE.

**Result: the defaults ARE the validated recipe.** kNN purity 0.804 vs the repeat origin's labels
(repeat's own self score 0.804), 0.803 in reverse (self 0.804), 0.767 vs the original sphere-origin
labels (repeat scores 0.768); per-cluster purities identical to 2 dp; island d' 3.58 / purity 0.71
(references 3.66 and 3.56, both 0.71); K=16. Differences are below reduce-order noise.
Build rc=0, so the marching removal (ce59f2062) compiles as well as syntax-checks.
Per-shell mean scale flat 0.93-1.03 for shells 2-20 (only DC and shell 1 deviate) -> inputs correct.

**Open defect found by this run:** the state stage shipped 8 states with occupancies
3874 / 2455 / 20946 / 49576 / **232** / **284** / 23492 / 4388. States 5 and 6 are reconstruction
artifacts at box 360 (state 6 is striped and has the wrong envelope; state 5 is a uniform amplitude
offset with a halo); the other six are clean ribosomes with structured, localised difference maps.
`min_neff=2000` did NOT prevent them: that floor sets kernel bandwidth, not whether a state is
reconstructed. Fix before deleting the state-stage alternatives: a hard occupancy floor that drops
or merges a sub-threshold state. Also state 4 alone holds 47% of the data.

**Artifacts:** montages `state_maps.png`, `state_diffs.png` (projection montage tool lives at
`~/state_montage.py` on the box, numpy+matplotlib only, no scipy). Volumes pulled to the Mac at
`~/Downloads/flex_paired_20260908/default/maps/` (16 x 186 MB + kernel tables).
Plan document: `doc/refactoring_notes/flex_pca_cleanup_plan.md` in the flex-delivery worktree.
Related: [[flex-pca-basis-composition]], [[flex-pca-marching-erodes-structure]], [[box-shared-with-elmlundho]].

**Cleanup executed 2026-09-12 (branch flex-delivery, worktree ~/SIMPLE_wt_delivery), 10 commits
d234f83df..b0ffe0e37:** environment switches went from 84 live reads to 14; ~1,900 lines deleted;
simple_flex_pca_ica.f90 removed. Paired + merge + one polish iteration is now the workflow, not an
arm. Deleted as measured-not-to-work: ICA, the external-basis probe, FINCH/external targets, the
smoothness rotation (pcrot), the Gaussian-shrink fallback, sigma2 recalibration, the focus and
solvent masks on the basis, the three compose variants, the inert CV rank cut, the AUTO-K
instrument, the two-job pcafit harness, the losing state-weight arms (STDZ/KMEANS/RELPATH/
STATE_HARD/NORMW/GMM off), raw-plane residency. Baked in: polar E-step, mixture width 16, 4-shell
deflation + background + dilation, cross-fit FSC ridge, full iteration budget, plain prior,
per-half contrast, pooled per-shell delivery kernel, GMM island/continuum factorisation.
STILL LIVE (deliberate): the 7 GPU/PREP_CHECK switches (cannot be tested on the CPU box -- do them
on the Blackwell box), SIMPLE_COV_MERGE/_MERGE_R, and the survivors that should become CLI
parameters: STATE_REC -> rec_states, COMPOSE -> compose_dirs, COMPOSE_CUT + CUT_SNR, HP.
**Build trap hit:** deleting a source file breaks the box build until `cmake ..` is re-run in the
build dir (the CMake glob has no CONFIGURE_DEPENDS) -- add it to every queue script that rebuilds.
**Editing trap hit:** index-based (`s.index`) block cuts in a long Fortran file spliced the wrong
region into simple_flex_pca_weights.f90 (190 spurious insertions, caught by `git diff --stat`);
use unique exact anchors, and check `git diff --stat` after every scripted edit.
Plan + status: doc/refactoring_notes/flex_pca_cleanup_plan.md on that branch.

**FINAL VERIFICATION (15:41, queue63/64/65 on the box, dirs default_/cleaned_/bare_b64_lp16):**
| run | environment | kNN vs the others' labels |
| pre-cleanup default | full SIMPLE_COV_* set | 0.804 (self) |
| cleaned tree | same env | 0.804 |
| bare cleaned tree | NO SIMPLE_COV_* at all | 0.803 |
Calibration a=2.509 in all three; K=16; loglik/particle -91.7784 to 5 dp; island d' 3.55-3.58,
purity 0.71; merge worst-axis match 0.46/0.17 identical. The defaults are in the code.
rec_states exercised BOTH ways on the command line: =yes 743 s with 6 maps, =no 48 s with 0 maps
(15x cheaper for embedding-only work). The occupancy floor dropped the same 2 states (232, 284
particles, 0.49% of the data) in every run.
**Two build traps, both hit and both now in the queue scripts:** (1) `cmake ..` must run before
`make` after deleting a source file (glob has no CONFIGURE_DEPENDS); (2) a NEW CLI key is written
by the args generator to build/lib64/simple/simple_args.f90 while CMake compiles
build/lib/simple/simple_args.f90 -- copy lib64 -> lib between cmake and make or the key is
silently rejected ("not allowed") at runtime.
**Masks (for the PCG arms, not yet run):** the NU envelope size has TWO knobs, both on nu_filt3D:
`amsklp` = the length scale in Angstrom at which the envelope is DEFINED (lp_smooth, default 8;
this is the "cutoff lp" that sets how much it encloses) and `nu_msk_sig` = the cut in MADs above
the solvent null (default 3.0). The filter-bank resolutions printed in the log are a diagnostic,
NOT an envelope input. At the defaults the 10028 envelope is 6.4% of the box (7.2e6 A^3, ~5.3 MDa
vs a 2.5 MDa ribosome); dilated 15 A + 15 A edge it is 10.6%. Files:
flex_h0901/solvent_envelope/{outvol_nu_envmask.mrc,solvent_mask_360.mrc}, also on the Mac at
~/Downloads/flex_paired_20260908/solvent_mask/. flex_pca on branch flex-pca-pcg takes
rec_backend=pcg, pcg_mskfile=<mrc> (read_and_crop'd to box_crop automatically; absent = the
mskdiam sphere), maxits_pcg, rtol and smpd_target (which derives box_crop). batch2
(fsitgl-xfer03p, 36 cores, 377 GB) has that branch built and mounts the same beegfs data.
