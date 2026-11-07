---
name: spliceosome-mest-refine3dauto-2026-10-02
description: 10180 46k M-estimation (soft flex weights in refine3D_states) + per-state refine3D_auto from the refined consensus; master flex_pca writes ONE-HOT weight files, the flexref5 recipe gives soft ones
metadata:
  type: project
---

User 2026-10-02: the PCG flex states are the better set; M-estimation (refine3D_states flex=no flex_weights=yes = each
state's reference reconstructed from ALL particles weighted by the flex posterior) followed by refine3D_auto has never been
tried. Measured first: master f9bbf72ce flex_pca (states5_c46k B, npreimages=5, rec_states default) wrote
flex_weights_state_NNN.bin that are ONE-HOT (max posterior 1.000 for all 46,783) = hard labels, useless for M-estimation.
The RNase P S4 fit (rnasep_split/S4/F_lp8, ~/rnasep_split.sbatch, ~/SIMPLE_flexref5, `npreimages=3 min_state_frac=0
rec_states=no` + pcg_mskfile) is genuinely soft: max posterior p10 0.53 / median 0.83 / p90 0.99. S4mest result
(~/rnasep_S4mest.sbatch, ~/SIMPLE_flexref6): soft 4.04/4.04/4.11 A vs hard 4.34/4.11/4.50 A.
Weight file layout: 512-byte header | 112-byte real64 scalars | float32 w(nrows) | int32 hard(nrows).
Pipeline `/mnt/beegfs/elmlund/afan/10180/mest_c46k_20261002/`: fit.sh (flexref5 flex_pca on consensus46k_refine3D_auto
p.simple, vol1 = its rec_final_state01, pcg_mskfile = its automask3D_state01, npreimages=5 min_state_frac=0 rec_states=no,
prints WEIGHTS flatness, sbatches both arms) -> arm.sh ARM=soft|hard (flexref6 refine3D_states flex=no flex_weights=yes|no
lpstart=8 lpstop=4.5 mkdir=no, then sbatches auto.sh per state) -> auto.sh (flexref6 selection state=S + refine3D_auto
nparts=4). All --mail-type=NONE. Fit job mestF 53549775.

**RESULT (2026-10-02 ~09:00):** fit weights median max 0.77 (24,452 of 46,783 weighted). refine3D_states final FSC0.143
soft 8.63/8.36/8.24/8.11/8.36 vs hard 8.91/10.07/9.06/9.06/8.77 (soft better by 0.3-1.7 A, as on RNase P). After per-state
refine3D_auto (hard selection of each arm's final labels) the gap vanishes: soft 7.45/9.21/7.55/7.35/7.25 vs hard
7.35/9.71/7.88/7.35/7.45. Deposited 10180 map is 7.2 A (user); our refined 46k consensus 5.91 A; locres core 8-10 A,
periphery 15-20 A = flexibility, not pose failure (user: particles are good; we refined from scratch).
User 2026-10-02: the 327k realigned set loses the conformational signal (only discrete/compositional states, same as
cryoDRGN on 10180) -> stay on the 46k subset for conformations. refine3D_auto must take flex weights: on ~/SIMPLE_flexref6
it already propagates them (cline_boot = cline keeps flex_weights/flex_state; main refine3D gets cline; simple_final_rec
passes both; builder loads column flex_state when nstates=1). Master has NO flex_weights/flex_state params at all (port
needed later). Test: mest_c46k_20261002/autow.sh -> soft/autow_state0S (selection states=1,2,3,4,5 merges to one state,
refine3D_auto vol1=soft/rec_final_state0S flex_weights=yes flex_state=S nparts=4), jobs mestw_st1-5 53555880-84; compare
vs hard-selected soft->auto 7.45/9.21/7.55/7.35/7.25 A.
