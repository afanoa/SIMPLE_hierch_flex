---
name: spliceosome-realign-deselected-unlocks-basis
description: "2026-09-14: re-aligning the DESELECTED particles back in (refine3D maxits=1 greedy volrec=no against the clean subset volume) took EMPIAR-10180 flex_pca from 2 of 20 reproducible components to 12 of 20 -- data volume was the binding constraint, not any tuning knob"
metadata:
  type: project
---

**The user's recipe (their idea, and it worked).** Compositional heterogeneity had been removed by a
2D selection, so flex_pca was fitting only the 46,783 selected particles. Re-align the deselected
ones against the CLEAN subset volume to recover poses for all of them, then fit flex_pca over the
full set while keeping the clean volume as the reference:

```
# A = the 2D project (ALL particles active), C = the abinitio3D (clean 46k volume)
simple_exec prg=refine3D projfile=A/p.simple vol1=C/recvol_state01.mrc \
  maxits=1 refine=greedy volrec=no lp=8 trs=10 \
  pgrp=c1 mskdiam=500 filt_mode=nonuniform nparts=20 nthr=16 mkdir=no
```

**`volrec=no` is CORRECT and deliberate** -- do NOT reconstruct from all particles. flex_pca models
deviations FROM vol1, so the reference must stay the clean 46k reconstruction; averaging the
deselected (partly junk) particles into it blurs the reference and mutes the very signal wanted.
I argued for `volrec=yes` and was wrong.

**RESULT (box_crop 64, lp 20, neigs 20, rec_backend=pcg + automask solvent envelope, 4 probe iters):**

| | prior 46,783 ptcls | after realign, 327,490 ptcls |
| components >=0.7 both halves | **2 of 20** | **12 of 20** |
| components >=0.85 both | 1 | **5** |
| mean min\|cos\| | 0.452 | **0.738** |
| spectrum flatness max/min | 3.68 | 4.68 |
| top components | one >0.85 | 1.00/1.00, 0.99/0.99, 0.97/0.95 |

All THREE prior runs (8_flex_pca_paired_br320, 9_..._v2, 9_..._v3) give identically 2/20 at 0.452,
so that was a stable ceiling, not an unlucky draw. Rank 20 on 46,783 particles is ~2,340 particles
per component -- the basis simply could not be determined. **The binding constraint was DATA VOLUME,
not box/band/rank/iterations/masks**, none of which ever moved cnga1 past 8/10 either.

**Parameter traps hit:** `filt_mode=nu` is INVALID -- the accepted set is
`none|uniform|fsc|nonuniform|nonuniform_lpset` (`simple_parameters_phases.f90:858`); abinitio3D had
used `nonuniform`. `pgrp` is REQUIRED by refine3D (no default) -- use the value the reference volume
was refined with (c1 here).

**Paths:** work dir `/mnt/beegfs/elmlund/afan/10180/allptcl_realign` (p.simple = copy of
`spliceosome_10180/2_abinitio2D/spliceosome_10180.simple`, all 327,490 active);
vol1 `spliceosome_10180/4_abinitio3D/recvol_state01.mrc`; solvent mask
`allptcl_realign/mask/automask3D.mrc` (from `prg=automask` on rec_final_state01_pproc, PCG hard
support 0.041). Geometry: box 320, smpd 1.699. Baseline runs are the `9_flex_pca_paired_br320*` dirs.
States were NOT reconstructed (STATE_REC=0) -- that is the open next step.
Related: [[afm-cryoem-1m-run-is-io-bound]], [[flex-pca-em-iterations-overfit-cross-half]].

**CAVEAT 2026-10-04:** these 327k poses were never globally searched (see [[refine3d-shc-virgin-field-trap]]). The 46k good particles are 112° off their consensus poses, so the "12 of 20" basis may reflect pose error rather than conformation. A proper realignment of the 281k against the 46k consensus is running in `10180/autostates_c46k/emb327k/align281k`.
