---
name: ole-relion-state-maps-2026-11-05
description: "OLE 12707 RELION per-state maps: batch2 ole_12707/relion; main set = SIMPLE consensus poses + CTF refine (cpose_w), 3.4-6.6 A; same frame as SIMPLE; tree differences present at ~half the tree-map contrast"
metadata:
  node_type: memory
  type: project
  originSessionId: a25e6f6d-c94c-4ba6-92f7-c1fb7553af83
  modified: 2026-11-05T06:25:23.096Z
---

2026-11-05, for the heterogeneity slide deck (parent session OLE RNA slides). Work dir batch2
/mnt/beegfs/elmlund/afan/ole_12707/relion (README.md, states_table.txt, occ_table.txt). 80,598 particles of the 31 6_hetero_tree
leaves; stacks are white; the RELION frame equals the SIMPLE frame (no mirror, unlike RNase P [[rnasep-relion-maps-are-mirror-of-tree-frame-2026-11-04]]).
CTF refine (2 optics groups): consensus 3.18 -> 3.15 A. No polishing (no movies downloaded).
Display maps (FSC-weighted, B = 0, unmasked, box 360, 0.954 A/px, C2): PostProcess/cpose_w_<leaf>/postprocess.mrc (A, delivered at the
top of ~/Downloads as ole_relion_<leaf>_<res>A.mrc), PostProcess/cons_w, branch_w_r1, branch_w_r2 (consensus 3.15, r1 3.18, r2 3.40 A).
B = local refine, worse everywhere ([[relion-small-state-local-refine-blurs]]): ~/Downloads/ole_relion_local_refine/.
Slide-region occupancy (occ_table.txt): cap and bridge track the tree maps across states (r = 0.99, half-map noise 0.01-0.03), with
about half the tree maps' spread. In A the r2 cap sits at about the contour level (0.8-1.2), so r1 vs r2 looks weaker than in the
tree maps. Related: [[ole-12707-setup-2026-10-23]].
