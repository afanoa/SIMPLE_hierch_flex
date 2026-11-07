---
name: rnasep-toy-ah-apo-holo-design
description: "the apo/holo x open/closed RNase P toy (\"ah\", 2026-10-16): 10 apo + 8 holo deposited conformers (EMD-74272..74279) docked into the apo frame (cc 0.95-0.98, same hand), scale-normalised, holo closed = 102/105/106 by the contact score; fast toy path (reused calibration, perturbed true poses, 100k) and the four board tests (plain, oracle site, truth-free site, tied two-layer)"
metadata:
  node_type: memory
  type: project
  originSessionId: 25d47be4-7d1d-48d8-abdf-97b1c3c946ab
  modified: 2026-10-16T15:12:31.838Z
---

Scripts (toy55, mirrored on batch2): `holo_maps_rnasep4.sbatch` (download + process_maps + dock_volpair into the apo
consensus frame -> maps160/holo_docked, kept OUT of maps160/docked so the real-data identification readouts stay
apo-only), `normalize_ah_maps.py` (least-squares scale of every docked map onto the apo consensus core -> maps160/ah_norm;
scale factors came out 0.95-1.10, so no brightness confound; the mean holo-minus-apo excess is NOT just rnpA: the
largest component is 85k A^3 (rnpA ~16k) = RNA rearrangement + protein, written as maps160/holo/rnpa_mask.mrc, 112k
A^3 dilated), `run_toyah_sim2.sbatch` (env AHTAG HFRAC FASTPOSES FASTERR BFAC SNR N: FASTPOSES=1 = true poses
perturbed by angerr 3 instead of the hour of refine3D, BFAC=0 SNR=0.0012 = the toy10c calibration, N 100k; CSV
rnasep_conformers_$AHTAG.csv, truth sim$AHTAG/gt_particles.txt with conformer >= 101 = holo, closed_ids.txt),
`run_toyah_boards2.sbatch` (A plain 7x7 Fortran board, B oracle focal 1x2 inside rnpa_mask, C derive_excess.py
candidates + focal 1x2, D Python tied two-layer board: --nlayer 2 --tie_layers --tie_support 0.5), `run_toyah_pose.sbatch`
(the Fortran board's local pose step at family scale vs the truth), `score_ah.py` (per-cell holo/closed fractions,
NMI, 2x2 best cells, layer-mass AUC). Holo closed/open by the P12/P10.1 contact score of the docked map: 105 (+1.40,
arm 0.99), 106 (+0.29), 102 (+0.23) closed; 101/103/104/107/108 open (3 of 8 ~ the paper's 43 %).
klem.py additions 2026-10-16: --tie_support (positive compact component, support = cross-half reproducible positive
part from round --tie_support_after, --tie_ncomp components) / --tie_support_mask FILE; the balanced E-step is
layer-aware with nlayer > 1 (cells balanced within each layer, the layer split free).
Related: [[rnasep-degenhardt-paper-closed-state-and-pipeline]], [[klem-f692-fortran-chain]], [[rnasep-toy55-design-2026-10-15]].
