---
name: rnasep-alternation-verdict-2026-10-14
description: overnight 2026-10-13/14 alternation (Pi-step sub-classes vs state maps, pose refinement, reconstruction) on 692k RNase P: coarse band rises (0.36 -> 0.43-0.48, from the pose refinement), every finer band falls below the tracked start; soft worst; state maps separate sub-classes by a median margin of 0.04 so no assignment is stable; tracked start remains the best fine-band partition
metadata:
  type: project
---

Five variants (`rnasep_coupled/alternation_v1..v5`, `run_alternation.sbatch`, metrics in
`doc/implementation_notes/coupled_2d3d_math_formulation_2026_10_13.md` section 13): flat vs 30-12 A score, beta 0.1 vs
0, tracked vs latent k-means start, hard vs soft. All hard variants at iteration 1: 117-30 A up from 0.36 to
0.43-0.48, 30-20 A down from 0.58 to 0.42-0.54, 20-12 A down from 0.38 to 0.22-0.32, 12-8 A down from 0.23 to
0.14-0.17; later iterations keep losing the fine bands (v1 it3 0.17/0.10). The coarse gain is the pose refinement
(main vs no-refinement 0.44 vs 0.34), not the re-assignment. The latent k-means start (v4) is the best it1 set and
drained one of its 15 states at the first Pi-step (zero sub-classes -> reconstruct3D "No active particles"). The v5
responsibility overlap matrix shows neither chain nor blocks: every state shares mass with 6-12 others because the
sub-class deviation correlates 0.22 with the best state map and the runner-up trails by 0.04.

**Why:** the state maps (pooled over the sphere, 10 A) carry only the smooth part of what distinguishes the tracked
slices; the fine structure cannot be scored, so each re-assignment moves sub-classes along the continuum toward the
nearest coarse state and reconstruction averages the fine part away.

**How to apply:** treat the alternation as a coarse-band pose polisher, not a state finder; keep the tracked
whole-sphere states as the fine-band partition; for a continuous description use the synchronised 2D-latent
coordinate ([[rnasep-sync-od-coordinate-2026-10-14]]). Related: [[mac-lsqr-vstep-proxy-invents-soft-gain]],
[[rnasep-basis-span-finding-2026-10-13]].
