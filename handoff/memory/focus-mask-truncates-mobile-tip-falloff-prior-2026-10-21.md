---
name: focus-mask-truncates-mobile-tip-falloff-prior-2026-10-21
description: "TRAP: variation focus masks thresholded at 20 % of the peak (arm_mask.py) leave a mobile region's blurred tip outside, and pcg_focus (mask = continuous solve window P) freezes it -> states with a SHORT mobile region; fix = falloff prior (core + exp(-d^2/2 s^2), arm_mask2.py --falloff 6) AND plain (unfocused) delivered maps"
metadata:
  node_type: memory
  type: project
  originSessionId: 25d47be4-7d1d-48d8-abdf-97b1c3c946ab
  modified: 2026-10-21T19:44:59.914Z
---

User 2026-10-21: "the mobile region seems shorter than the ones in the paper, which signals masking". Confirmed with radial_extent.py
(routes/toy55): closed-arm density relative to the deposited closed map at 95-100 / 100-105 / 105-110 A: unit-7 state 6 (focused
refine, cheap-children mask) 0.70 / 0.41 / 0.21, but a plain reconstruct3D of the site route's top 2.5k 1.16 / 1.05 / 0.89 (the data
carry the tip; the 232 A sphere covers it). Mask coverage there: cheap-children 20 %/3 px 0.60 / 0.26 / 0.07, f692st root flex-state
mask 0.19 / 0.03 / 0.00; f692rf node masks ~1.0 (their states are blends, not truncated).
pcg_focus: rec3D_pcg_strategy set_pcg_solve_support -> set_window_volume(mask), (P H P) u = P b with P = the mask VALUES (only values
<= -min are floored) -> a soft edge acts as a falloff prior. Gaussian-blurring a small core dilutes the core (0.4 peak) - use a distance
falloff: arm_mask2.py --frac 0.2 --dil 3 --falloff 6 -> tip weights 1.00 / 0.95 / 0.80 at 21.6 % of the support (hard 10 %/6 px:
0.96 / 0.83 / 0.59 at 30 %). Significance-only cores (negative-tail null) cover 88-100 % of the sphere (global reproducible differences).
**Why:** focused states silently lose the moving part's extremity; the user judges maps against the paper's arm length.
Even WITH the falloff mask the focused stage-B MAPS still damp the tip (factorial session, f692st leaf r1_s5: focused 0.58 / 0.39 / 0.23
at 95-110 A, plain reconstruct3D of the same particles at the node poses 0.87 / 0.80 / 0.67; r1's falloff mask covers 0.45 / 0.28 / 0.16):
focus is fine for refinement, never for delivered maps.
**How to apply:** rflex13 builds node masks with the falloff prior (MFRAC 0.2 MDIL 3 MFALL 6); check radial coverage of any new focus
mask on the moving region before refining; DELIVER plain reconstructions (routes/toy55/rflex_final2.sbatch: one reconstruct3D per node,
one state per leaf, at the node's stage-B poses; the tree-weighted maps and score embedding take those half-maps);
see [[minority-soft-em-fixed-poses-no-rescue-2026-10-21]].
