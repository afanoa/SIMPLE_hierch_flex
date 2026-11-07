---
name: simsites-build-traps-2026-10-21
description: "TRAPS for RNase P synthetic sets: maps160/variants/domain_mask.mrc IS the arm's swept volume (not an independent site); simulate_particles soft-masks the VOLUME at mskdiam, so toy sims at mskdiam 200 clipped the closed-arm tip (use 232); simsites = rigid conf49 core + 3 disjoint sites"
metadata:
  node_type: memory
  type: project
  originSessionId: 1ca67e82-4d15-4725-b95d-0ee5719be9f4
  modified: 2026-10-21T22:04:27.309Z
---

Found 2026-10-21 building `simsites` (factorial site test, `doc/implementation_notes/factorial_site_layer_results_2026_10_21.md`):
- `maps160/variants/domain_mask.mrc` (make_variants.py "mobile domain" = largest high-variance component across the 10 conformers) is the
  closed arm's swept volume: 86 % of arm_readout's closed-arm region and 49 % of the open-arm region inside it, centroid 3 A from the
  closed arm. A "deleted domain" site from it and an "arm open/closed" site are the SAME place. simsites site A = a compact lobe of the
  conformer-49 core instead (ball r 17 A, 6.0 % mass, 23 A from the other sites; make_sites_maps.py).
- `simulate_particles` applies `vol%mask3D_soft(msk)` to the volume before projecting (simple_commanders_sim.f90): the toy10c / sim55
  runs at mskdiam 200 (radius 100 A) cut the closed arm, which reaches ~108-114 A. simsites uses mskdiam 232.
- The conf62 blob is exactly conf62 - conf46 x mask200^2 (make_variants multiplies by the soft mask twice); it sits outside conformer
  49's density (3.4 % of its mass).
**How to apply:** check any new site mask against arm_readout's regions before using it as an independent site; simulate RNase P at
mskdiam >= 232. Related: [[rnasep-klem-mask-clips-closed-arm]], [[rnasep-toy55-design-2026-10-15]].
