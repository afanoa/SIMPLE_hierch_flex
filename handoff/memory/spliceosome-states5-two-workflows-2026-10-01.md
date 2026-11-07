---
name: spliceosome-states5-two-workflows-2026-10-01
description: "2026-10-01 10180 46k 5-state comparison: A = refine3D_states flex=yes, B = flex_pca pcg+solvent mask; each then refine3D_auto per state; dirs, scripts, job ids"
metadata:
  node_type: memory
  type: project
  originSessionId: b8bc4431-3743-4b17-aafd-6d31e354ea19
  modified: 2026-10-01T17:48:19.555Z
---

Root `/mnt/beegfs/elmlund/afan/10180/states5_20261001/` (master build ~/SIMPLE f9bbf72ce). Both workflows start from a copy
of `spliceosome_10180/4_abinitio3D/spliceosome_10180.simple` (46,783 ptcls, poses match recvol_state01; see
[[spliceosome-project-is-4-abinitio3d]]), sigma_est=global ([[sigma2-group-fails-on-selected-subsets]]), mskdiam=400 pgrp=c1.

- `A_refine3D_states/` (job s5A 53534239): refine3D_states flex=yes nstates=5, default min_state_frac 0.1, nparts=10 nthr=8.
  Per-state start maps read as `recvol_state0S.mrc` in the workflow dir (naming NOT verified before launch).
- `B_flexpca_pcgmsk/` (job s5B 53534240): flex_pca rec_backend=pcg pcg_mskfile=allptcl_realign/mask/automask3D.mrc
  npreimages=5, nparts=10 nthr=8. Start maps `flex_pca_state_00S.mrc`.
- Step 2 (inside `pipeline.sh`, MODE=A|B): for each state map found, sbatch `auto_state.sh` (W,S,VS env) ->
  `auto_state0S/`: copy p.simple, `selection oritype=ptcl3D state=S mkdir=no` (keeps state S, relabels it 1),
  refine3D_auto vol1=VS nparts=4 nthr=8. Logs: pipeline.<jobid>.out, step1.log, auto_state0S/run.log.
Recipe for per-group refine3D_auto came from `~/rnasep_mbauto.sbatch`.

**Outcome (2026-10-01 ~14:50):** B per-state refine3D_auto 8.91/9.71/7.77/7.55/7.25 A (states 1-5; st2 = 3,952 ptcls,
st1 degraded 7.77->8.91 over 5 iterations) -- no better than the 46k abinitio3D consensus itself (CORRECTION: the consensus final map is 7.15 A FSC0.143 /
9.21 FSC0.5, cFAR 0.36, per 4_abinitio3D/RESOLUTION_FINAL_STATE01; 8.91 was an intermediate cropped-stage number). User stopped both workflows (s5A cancelled mid refine3D_states): "we need a better
starting point for multi state sorting". Next: `10180/consensus46k_refine3D_auto/` (job c46k_auto 53538090,
refine3D_auto on all 46,783 from the abinitio3D project + map, same settings), then re-split into states from it.
Fallback lever if the consensus stays ~9 A: the 327k realigned set (allptcl_realign).

c46k_auto progress 15:42: FSC0.143 5.85 (reg pass, band-capped) -> 5.78 -> 6.25 A; cFAR 0.34 -> 0.30 -> 0.31 -> 0.21 (anisotropy worsening).
c46k_auto DONE 16:21 (1 h 30): final 5.91 A FSC0.143 / 7.66 A FSC0.5, cFAR 0.24 (was 0.36 at start) -> resolution better, anisotropy worse. Refined project: 10180/consensus46k_refine3D_auto/p.simple.

**Rerun from the refined consensus (user: map "looks fine", proceed), 16:40:** root `10180/states5_c46k_20261001/`, same
A/B design; P0 = consensus46k_refine3D_auto/p.simple, B vol1 = its rec_final_state01.mrc, B pcg_mskfile = its
automask3D_state01.mrc. Jobs c46A 53542126, c46B 53542127; per-state jobs named c46A_stN / c46B_stN.
refine3D_states final per-state maps are rec_final_stateNN.mrc (simple_final_rec.f90:193); A step 2 now reads those.
