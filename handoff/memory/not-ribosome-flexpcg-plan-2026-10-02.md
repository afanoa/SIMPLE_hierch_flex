---
name: not-ribosome-flexpcg-plan-2026-10-02
description: not (80S + Not module) flex PCG with a NU-envelope solvent mask, following the 10028 recipe that worked; why refine3D_states' internal flex failed here
metadata:
  type: project
---

User 2026-10-02: the not refine3D_states states show heavy blurring (far more states) and the "junk" states are not
believable (good particles; states 4-6 had joint overlap 0.004-0.12 = references never stabilised). Asked: flex on the
ORIGINAL consensus (not/4_abinitio3D, 452,386 ptcls, 3.16 A FSC0.143, cFAR 0.21, mskdiam 400, box 420 @ 1.3) with PCG +
the CORRECT solvent mask, get states without alignment, then hard/soft refine them.
"Correct mask" = the NU envelope from `prg=nu_filt3D vol1=<odd_unfil> vol2=<even_unfil> smpd mskdiam` (writes into
1_nu_filt3D/ under cwd: outvol_nu_envmask.mrc), as for 10028 (flex_h0901/solvent_envelope, 6.4% of box). NOT refine3D's
automask3D_state01 (that is the 20-A envelope-FSC mask). Job not_envmask 53556078 -> /mnt/beegfs/elmlund/afan/not_flexpcg/solvent_envelope.
10028 success (flex_verify_20260915/R_P_pcg -> R_Ppcg_gridst_ax0, ~/flex_verify/flex_verify_arm.sh ribo): standalone flex_pca,
cryoDRGN/cryoSPARC-refined poses, rec_backend=pcg pcg_mskfile=NU envmask (support 7.3% at box 64), neigs=20 n_probe_iters=6
box_crop=64 lp=16 box_rec=360 mskdiam=320 npreimages=8 state_axis=0 min_neff=2000, NO population floor; deconv K=16 ->
6 macro-clusters -> 8 states of 6.4k-40k, state maps 4.05-6.18 A straight from flex labels (no re-sorting).
Contrast with the failed not runs: refine3D_states' internal flex forces rec_backend=gridding (no mask possible),
applies the 10% floor + random fill, then re-sorts by image likelihood. Lesson: standalone flex + envelope + no floor,
then refine the delivered states.
Mask result: master nu_filt3D needs `nu_envmsk=yes` (param default 'no'; the 10028 mask predates that) and `mkdir=no`
to write in cwd. not envelope = solvent_envelope/outvol_nu_envmask.mrc: 2 components found, 1 kept, 5.8% of box 420,
9.5e6 A^3 (~2x an 80S+Not, same ratio as 10028), r99 158 A, fully inside mskdiam 400. Fit job not_flexpcg 53556199
(`not_flexpcg/fit.sh`, flexref5, vol1 = 4_abinitio3D/rec_final_state01, rec_backend=pcg pcg_mskfile=envelope neigs=20
n_probe_iters=6 lp=16 npreimages=8 state_axis=0 min_neff=2000 min_state_frac=0 mskdiam=400 sigma_est=global rec_states=yes,
default box_crop).
FIT RESULT (2.5 h): 8 states 120k/129k/54k/30k/22k/45k/28k/23k, no floor/random fill; flexref5 delivers state maps as
Epanechnikov kernel-regression "DIFFMAP NYSTROM PRE-IMAGE" maps -> blurry eo-FSC 9.1/13.0/8.5/6.6/9.4/6.0/7.1/11.1 A
(10028's hard-label maps were 4-6 A); weight files = Gaussian-classifier posteriors on the hard labels (20 comps).
Next (user chose BOTH): hard-selected refine3D_auto all 8 (`not_flexpcg/auto_hard.sh`, hard_state0S, nparts=4, jobs
notH_st1-8 53558542-49) + weighted refine3D_auto (`auto_w.sh`, w_state0S, all 452k merged, flex_weights=yes flex_state=S,
nparts=6) batch 1 = states 6,4 (notW_st6/4 53558550-51); remaining weighted states 7,1,3,5,8,2 to start as cores free.
Spliceosome weighted refine3D_auto FINAL: 7.35/7.66/7.06/6.88/7.06 vs hard-selected 7.45/9.21/7.55/7.35/7.25 (all better).
STOPPED 2026-10-02 ~15:20 (user): "the consensus volume is extremely biased, so we will fix the sampling first" (not has one
strongly preferred view; consensus cFAR 0.21). Cancelled notH_*, notW_*, default r2 53555834 + workers. Before the stop,
hard-selected refine3D_auto finished for states 5/7/8 only: 11.87 / 9.93 / 9.58 A (FSC0.143) = no better than the
smeared flex maps. User judged the flex PCG states distinct but extremely smeared; the UMAP states panel showed states
merging different islands (state 1 = right island + part of upper continent). Basis itself reproducible (20 comps,
principal cos 0.84-0.99). Next for not: fix the orientation sampling / view bias, then redo states.
