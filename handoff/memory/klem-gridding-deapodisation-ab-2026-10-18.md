---
name: klem-gridding-deapodisation-ab-2026-10-18
description: Kaiser-Bessel gridding + deapodisation in the Fortran board (cfg key `gridding 1`, fboard GRIDDING=1): A/B on 50k synthetic parent boards doubles cell cross-half reproducibility (0.16 -> 0.37) and improves site purity (tethered2 best 0.42 -> 0.29), neutral on arms (0.55 -> 0.62), but 10x slower (3.5 -> 37 min) as implemented
metadata:
  type: project
---

Implementation (port simple_klem_core.f90: klem_set_gridding / klem_gridding_on / klem_deapod_envelope /
klem_kbinterp; M-step 3x3x3 KB stencil (KBWINSZ 1.5, KBALPHA 2) with Friedel mates; E-step KB interpolation instead of
trilinear; deapodisation = real-space division by the discrete stencil's envelope squared (M-step + E-step), capped at
4, applied through apply_focus_d with G0 = 0 before the focus; cfg key `gridding` in simple_klem_distr.f90 (master
and worker) and simple_commanders_klem_board.f90; fboard.sh env GRIDDING). Built in ~/SIMPLE_klem2 2026-10-18 10:54
(clean), frozen at rnasep_klem/gridding_ab/klem2_bin. Default OFF: every chain run is unchanged.
A/B (gridding_ab.sbatch, 50k random subset, 7x7+1, random init, 8 bands 20->8 A, 16 threads): tethered2 reproducibility
0.157 -> 0.376, argmax-cell site purity median 0.68 -> 0.58 (best 0.42 -> 0.29); arms2 0.173 -> 0.366, arm std ratio
median 0.82 -> 0.83 (best 0.55 -> 0.62). Wall 3.5 min -> 37 min (27-point stencil in the M-step per band coefficient;
needs a precomputed-weights / per-particle vectorised insertion before production use).
Related: [[klem-kernel-latent-em-2026-10-15]], [[klem-focus-step-was-the-board-floor]].

**Mode 2 = trilinear 2x2x2 gridding + its envelope (w 1/4 1/2 1/4):** reproducibility 0.42 (tethered2) / 0.48 (arms2),
site purity median 0.51 / best 0.28 (best of nearest / KB / linear), arms neutral (0.86 / 0.63), wall 7 min = 2x nearest
(KB 10x). Recommendation: `gridding 2` as the chain default candidate; RNase P DEV chain f692g launched with GRIDDING=2
to compare the landscape with the dev baseline.
**User decision 2026-10-18 12:xx: "i dont really want to switch gridding"** - keep it OFF (default); the f692g dev chain
was cancelled. Do not propose enabling it again unless asked.
