---
name: agnostic-order-protein-mask
description: "2026-09-27: symmetry/axis-agnostic protein support from half-map order works -- band 1.35-1.8 x FSC0.143, Otsu on local e/o corr, membrane normal = least-spread axis of the disordered set, belt slab from its profile; cnga1 normal exact (0.0 deg), IoU 0.984 vs the hand mask; PfCRT (C1) gives a flat belt at 43.9 deg tilt"
metadata:
  node_type: memory
  type: project
  originSessionId: 2590d10f-7aff-4aff-b2cc-d9d85f962ad3
  modified: 2026-09-28T03:43:04.056Z
---

Script `agnostic_mask.py` (+ `pngout.py`, numpy-only PNG writer -- batch2 and the Mac have no usable matplotlib) in
`/mnt/beegfs/elmlund/afan/cnga1/flex_region/`; outputs `agnostic/{cnga1,pfcrt}_protein_agnostic.mrc` (native box) and
`*_agnostic_frame.png`. Inputs: the refinement's `recvol_state01_{even,odd}_unfil.mrc`, `automask3D_state01.mrc`,
`RESOLUTION_FINAL_STATE01`. Window sigma = the band's coarse edge; belt = disordered voxels inside the slab where the
disordered-count profile along the normal is >= 25 % of its peak; planarity lambda_min/lambda_mid < 0.5 gates it
(soluble protein -> automask untouched). ~1 min on 8 cores for both.

cnga1: band 5.95-4.46 A, Otsu 0.505, disordered 35.6 % of automask, normal (1,0,0) = 0.0 deg, centre (-10.0, 0.1, 0.1) A,
planarity 0.125, slab 56 A, belt 34.4 %, **IoU 0.984** vs labels_prot1_rest2 (the hand mask behind v2_R1_prot).
PfCRT (`/mnt/beegfs/elmlund/pfcrt/9_refine3D_auto`, box 300 @ 1.035, FSC 3.93): band 7.07-5.31 A, Otsu 0.353,
disordered 50.7 %, normal tilted 43.9 deg from box z (C1 refinement, arbitrary frame), planarity 0.225, slab 52 A,
belt 43.3 %, protein 56.7 %; slices show a clean belt ring and keep the extra domain above the bundle as protein.
IoU only 0.428 vs the old `flex_focus_protein_native.mrc` (arbitrary 100 kDa density threshold, see
[[flex-pca-complement-block-and-step0]]): agnostic-only 28.6 %, reference-only 9.0 % of the automask.
The CNGA1 component montage (before/after/belt block) is `figs/cnga1_pcs_side_top.png`; the belt block shows
textbook m=1 (shift) and m=2 (oval) belt deformations. Related: [[cnga1-micelle-leak-is-negligible]],
[[cnga1-angular-harmonic-label-weak]] (its axis attempt failed on the centre; this one does not).

**ROBUSTNESS STUDY 2026-09-29 (mask_v2.py, scored against the docked models).** The current rule (band 1.8-1.35 x
FSC0.143, Otsu, principal axes of the disordered piece) is best or joint best: normal error 0.0 deg cnga1 / 4.0 deg PfCRT
(vs the mean CA(i)->CA(i+4) helix direction), model kept 99.9 / 99.7 %. Robust-slab (RANSAC) normal: 3.0 / 6.5 deg --
no better. Data-chosen band (max Otsu separation): same on cnga1, WORSE on PfCRT (finer band, transporter 98.5 %, 0.9 % of
model into the belt). 3-band vote: marginal gain (PfCRT model in belt 0.24 -> 0.09 %). Noise-calibrated threshold
(99th pct of a phase-randomised null, 0.12-0.14) FAILS: the inner belt has real shared signal, mask swallows 78-82 % of
the automask. **Use a relative split (Otsu), not significance.** Domain split: erosion (<= 11 A) finds no neck; cross-
section profile has no 2x constriction; density-saddle watershed splits cnga1 into its 4 PROTOMERS and PfCRT wrongly
(Fab piece takes 69 % of the transporter) -- the Fab binds flat on the extracellular face, no saddle at this resolution.
What works for the helix region: protein mask x (belt slab +/- 10 A margin, fixed a priori) = tm_slab.py
(`agnostic_v2/<ds>_tm_slab_mask.mrc`): PfCRT transporter CA 97.1 % inside, Fab CA 13.8 % (its Fv tip binds at 27-55 A);
cnga1 keeps 66 % of CA (cytoplasmic domains out by design). 6UKJ models only the Fab Fv (~1780 atoms): the Fab constant
domains are ordered, UNMODELLED density -- that, not detergent, is the "40 % of mask far from model" on PfCRT.
Runs: pfcrt/flex_soft/p_auto (automask control), p_block (labels, belt reference), p_tm (helix region); cnga1 c_tm.
**PfCRT arms scored (score_pfcrt2.py, 2026-09-29):** automask p_auto: 57.5 % of component variance on the belt, leak
R^2 13.2 % (like cnga1's 14 %); order mask p_hard: belt 0.3 %, leak 0.14 %, transporter/Fab/shell 25/14/61 %; helix-
region p_tm (slab +10 A): belt 0.7 %, leak 0.06 %, transporter 45 % / Fab 1.5 %, components stay inside the TM region;
p_block (protein + rest blocks): protein reproducible 3/6 vs p_hard 0/6 (cnga1 s_soft1 7/10 vs s_hard 2/10) -- an
explicit nuisance block improves protein-component reproducibility at nparts=1. p_block's rest block is 82 % solvent.
**CNGA1 R1 state SNR (state_snr.py):** 15-8 A band, TM protein signal fraction 0.86-0.93 (SNR 6-14) for pc1-4, belt and
outside at noise; <15 A band has SNR 2-5 even OUTSIDE the automask -> a global low-res covariate (scale/contrast/ice/
defocus mix) differs between states; do not read <15 A state differences as protein motion without regressing it out.
Selection-bias null running: compose route, random protein-confined directions (flex_region/null_rand) vs R1 real pc3.
**CNGA1 helix-region run c_tm** (flex_region/c_tm, cnga1_tm_slab_mask): leak 0.59 % (s_hard 0.94 %), all 10 columns on
protein, subspace shares its 6 leading directions with s_hard (cancorr 0.89..0.69), 1/10 reproducible at nparts=1.
Overnight deliverables: ~/Downloads/cnga1_overnight_20260929/ (PNGs + ctm_states/ box-160 state maps).
