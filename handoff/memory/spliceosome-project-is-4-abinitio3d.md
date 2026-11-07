---
name: spliceosome-project-is-4-abinitio3d
description: the 10180 spliceosome 46k flex_pca runs MUST use 4_abinitio3D/spliceosome_10180.simple (poses match recvol_state01); 3_selection carries the IMPORTED star poses and mismatches vol1
metadata:
  type: project
---

2026-09-16: every 10180 flex_pca run since 09-15 (sel46k_*, resume_*, verification S_* arms) copied
`spliceosome_10180/3_selection/spliceosome_10180.simple` (46,783 selected) with
`vol1=4_abinitio3D/recvol_state01.mrc`. 3_selection's ptcl3D poses are the ones imported from
`consensus_10180_simple.star` (== 1_import_particles), NOT the abinitio3D poses the volume was built
from. Pose-frame mismatch -> per-shell mean-scale ratio 6.4/4.0/2.8/-1.2/... (must be flat near 1),
PC1 became "compositional", state maps 25-42 A. The 09-06..08 runs with conformational PC1
(8_/9_flex_pca_paired_br320_v3) used `4_abinitio3D/spliceosome_10180.simple` (same 46,783 particles,
matching poses; mean scale 0.86-1.08).

**How to apply:** for 10180 always copy the 4_abinitio3D project; check the mean-scale table in the
first 100 log lines before believing any basis ([[flex-pca-mean-scale-is-the-first-diagnostic]]).
The sigma2 selection issue ([[sigma2-group-fails-on-selected-subsets]]) is the same on both projects.
See [[flex-pca-fork-map-2026-09-15]].
