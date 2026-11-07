---
name: rnasep-s4-four-arm-support-test-2026-10-08
description: "2026-10-08 S4 (238k, lp 8, HP 15 band-sort, SIMPLE_bandsort build): regular automask vs 3-region multibody vs NU envelope vs locres-weighted support; multibody basis far more reproducible (18/18, lead |cos| 0.92 at it 4 vs 0.59-0.76), NU supports add nothing, every arm is a continuum (disc UMAP, flat spectra, GMM cuts are equal-mass slices)"
metadata:
  type: project
---

Runs `/mnt/beegfs/elmlund/afan/rnasep_autostates/S4/arm_{regular,multibody,numask,nulocres}` (script
`~/autostates/as_arms.sbatch`, env ARM; ~35-45 min each, 4 in parallel on csbdevel). Same build (SIMPLE_bandsort), lp=8,
npreimages=3, SIMPLE_COV_EMBED_HP=15, umap=yes; multibody = pcg_mskfile=rnasep_regions/G10cov_labels_B3.mrc neigs=6
rec_states=regions DECONV=0; numask = nu_envmask3D_state01.mrc (Jaccard 0.82 with the automask); nulocres = automask x
clip((locres-4)/8) from recvol_state01_nu_locres.mrc (`rnasep_hires/S4/nu_locres_weight_state01.mrc`; locres median 8 A
inside the automask, 78% worse than 5 A). Figures: `doc/implementation_notes/s4_arms/*_umap_states.jpg`; maps copied to
`~/Downloads/s4_arms_20261008/`.

| arm | pairs | lead |cos| it1 -> it4 | S(t) it1 -> it4 | GMM |
|---|---|---|---|---|
| regular | 10/10 | 0.90 -> 0.59 | 0.82 -> 0.43 | mode 162k + two tails 38k (support 0.16) |
| multibody | 18/18 | 0.98 -> 0.92 | 0.92 -> 0.69 | 1 macro-cluster; region 3: 5/6 cols >= 0.85, region 2: 1/6, region 1: 1/6 (3/6 >= 0.7) |
| numask | 10/10 | 0.90 -> 0.60 | 0.80 -> 0.46 | K=2 equal-mass cut |
| nulocres | 10/10 | 0.93 -> 0.76 | 0.71 -> 0.52 | K=2 equal-mass cut |

**Reading:** the region-masked (multibody) support is the only one that gives a stable basis on RNase P; NU-derived
supports are no better than the automask. S4 is a continuum in every arm (UMAP = uniform disc, merged eigenvalues flat
1.3-2.0k non-monotone, GMM boundaries are straight equal-mass cuts): no island inside S4 at whole-molecule or region level
before the gap test on region axes. Probe iterations erode the basis in all arms (S(t) falls), as on cnga1/TRPM4.

**How to apply:** SUPERSEDED 2026-10-14: the user does NOT want the region/multibody basis for RNase P ([[no-region-basis-for-rnasep]]); do not use or test it. Keep: do not spend runs on NU masks; read S4 by ladders
and nulls, not by GMM counts; keep n_probe_iters low. Open: ladder + spectrum-matched null readout per region axis for the
3.3 decision, and anchoring S1/S2/S4 to the deposited open/closed models (plan
`doc/implementation_notes/rnasep_ensemble_plan_vs_degenhardt_2026_10_08.md`). Related:
[[rnasep-periphery-ladder-2026-09-29]], [[flex-band-sort-removes-selection-bias]], [[flex-placement-not-latent-loses-small-states]].

**Ladder readout (same day, `S4/ladders/{reg_q2,mb_q13,mb_q1,null1_q1,null2_q1}`, 12 bins, HP 15, scorer
`~/autostates/as_score.py G10cov_labels_B3 automask`):** end-to-end SNR real/null in the unbiased <20 A band: multibody q13
in region 3 2.21/0.17 (energy beyond bias 0.63, localised, KEPT); multibody q1 in region 1 1.49/0.45 (0.42); regular q2
1.26/0.45 in r1, 1.00/0.17 in r3, 0.80/0.28 outside (0.30-0.40, a mixed global mode). In the 15-8 A sort band region 1 is
at bias level (76 vs 69) and region 3 72 vs 31. Adjacent-bin SNR is at null level for every axis (0.1-0.4), so 12 bins are
NOT resolved on 238k particles in the unbiased band; a region-3 movie is end-to-end plus ~2 intermediates.
**Low-rank placement (arm_mb_n1 / arm_mb_n2: 1 or 2 columns per region, npreimages 8, deconv on, HP 15):** 3/3 and 6/6
columns matched (lead |cos| 0.94 / 0.93, S(t) 0.88 / 0.74); per-axis deconvolution signal fraction 0.70-0.74 (regular arm
0.75-0.80), i.e. the latent is a genuine broad continuum, not per-particle noise; GMM gave 8 states of 20-68k with
frac_in_support 0.09-0.15 (slices of one Gaussian). The mb_n1 UMAP shows filaments and spikes: an ARTEFACT of the built-in Fortran UMAP on 3-column input (umap-learn on the
same coordinates and on 3/6/18-D Gaussians gives a disc, scratch control 2026-10-08); z1 vs z2 is one Gaussian. Full-band variants (arm_mb_n1_full, arm_mb_n2_full, HP unset) were run
because HP 15 removes the band where a 30 A domain motion has its per-particle signal.
**Full-band low-rank arms:** 3/3 and 6/6 columns at |cos| 0.98-0.99, signal fraction 0.82-1.02; GMM = core 182k + halo 55k; the halo is a CONTRAST tail (maps scaled 1.03-1.04 vs consensus, difference evenly spread over regions, halo-core e/o 0.13, NMI view/defocus 0.005 vs null 0.0001), i.e. the nuisance the HP-15 band-sort removes, not a conformation.
**S1/S2 chains (`~/autostates/as_mbchain.sbatch`, `rnasep_autostates/S{1,2}/{arm_multibody,ladders,null_r*,leaders.txt}`):** same fit
quality (18/18, lead 0.93/0.90; leaders S1 q6/q7/q15, S2 q3/q7/q13, 97-98% in-region). EB end-to-end <20 A: S1 r1 0.29, r2 0.35,
r3 <0; S2 r1 0.14, r2 0.36, r3 0.29; S4 r1 0.42, r3 0.63. No axis passes 0.5 on more than one state; adjacent bins at null
everywhere. TRAP: the chain's `[ -e null_r$b/... ] ||` guard SKIPPED building nulls because `rnasep_autostates/S{1,2}/null_r*` already
existed from the 2026-10-03 per-state auto-region runs (different label numbering), so "null_r1" lives in r3 (S1) / r2 (S2).
Pairing was done by measured energy share against the G10 labels (fcrop'd to box 72). Build fresh nulls into a new dir name
whenever the label set changes; always measure where a null lives before pairing.
