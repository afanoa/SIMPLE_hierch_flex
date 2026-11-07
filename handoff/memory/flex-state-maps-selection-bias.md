---
name: flex-state-maps-selection-bias
description: "2026-09-29: on-axis state maps are inflated by selection bias; against a SPECTRUM-MATCHED random-direction null (TM 15-8 A SNR 3.8) cnga1 R1 pc1/pc2 exceed it 2.5-3.6x (mostly real), pc3/pc4 only ~1.6x (~1/3 real, ~2/3 bias); a flat-spectrum null (SNR 18) overstates the bias"
metadata:
  type: project
---

Setup (`flex_region/`): null = 10 random smooth 8 A-band directions confined to R1's protein support (make_null.py ->
`null_rand/` with R1's probe meta), embedded through the compose route (`SIMPLE_COV_COMPOSE=null_rand`,
`SIMPLE_COV_COMPOSE_CUT=0` -- WITH the cut the run stops: "signal subspace 1 of 10 ... fewer than 2 columns survived");
states along random direction 1 (states_null2_pc1). Control through the identical route: R1's real basis
(`states_cR1_pc3`) reproduces the original R1 pc3 readout exactly (TM 15-8 A SNR 6.30 vs 6.31). Readout = state_snr.py
(D = S5-S1 vs N = ((E5-O5)-(E1-O1))/2).

| direction | TM <15 A | TM 15-8 A | belt <15 A | outside <15 A |
| random | 26 | **18** | 1.6 | 0.15 |
| R1 pc1-4 | 143-349 | 6-14 | 7-51 | 2.2-5.4 |

**Why:** sorting particles by their own images along a direction aligns their noise with that direction; averaging ~48k
of them turns it into a clean difference with the direction's shape. Reconstruction noise (the half-map null) does not
see this. **How to apply:** never cite state-map fine detail (or "above half-map noise") as evidence of heterogeneity;
evidence is basis-level: cross-half component reproducibility and the compose-cut population-variance test (real 10/10
vs random 1/10). A spectrum-matched null (null_spec/, states_null3_pc1) was queued because the flat-spectrum random
directions have more 15-8 A power than the real components. Unbiased state validation needs sorting and reconstruction
from independent data per particle (e.g. frame halves). Related: [[flex-pca-state-reproducibility-gap-test]] (same
circularity), [[cnga1-micelle-leak-is-negligible]].

**SPECTRUM-MATCHED NULL (make_null_spec.py -> null_spec/, states_null3_pc1):** random directions shaped to R1's mean radial
spectrum (power share 0.756/0.233/0.011 vs R1 0.746/0.243/0.011; direction 1 overlaps R1's subspace 3.8 % -> null slightly
inflated, comparison conservative). SNR TM <15 A 103, TM 15-8 A **3.8**, belt <15 A 4.9, outside <15 A 1.9. vs R1: TM 15-8 A
pc1 13.6, pc2 9.7, pc3 6.3, pc4 6.0; outside <15 A pc1/pc2 5.4, pc3 2.2, pc4 3.4. Energy beyond bias ~ 1 - (null+1)/(real+1):
pc1 ~0.67, pc2 ~0.55, pc3 ~0.34, pc4 ~0.31. The "global low-res term outside the automask" is mostly GENERIC to sorting
on coarse directions (contrast/defocus mix follows), only pc1/pc2 exceed it. Always use a spectrum-matched null.
**HELIX-REGION RUN c_tm (2026-09-29 03:5x):** null4 = TM-slab-confined random directions shaped to c_tm's spectrum
(null_spec_tm/, states_null4_pc1; power 0.703/0.285/0.012 vs c_tm 0.662/0.325/0.014): TM 15-8 A SNR 4.7, TM <15 A 155,
belt <15 A 6.0, outside 2.05. c_tm states TM 15-8 A: pc1 13.5, pc2 9.6, pc3 10.0, pc4 7.1 -> energy beyond bias ~0.61 /
0.46 / 0.48 / 0.30; TM <15 A 176-219 (~null: coarse TM differences mostly generic); belt <15 A 13-37 (2-6x null: micelle
density covaries with the TM states). Shapes at 8 A: pc1 radial periphery-vs-centre, pc2 4-fold alternating (twist-like),
pc3 4-fold centre-vs-ring (iris-like), pc4 dipole (shift-like, treat with caution).
