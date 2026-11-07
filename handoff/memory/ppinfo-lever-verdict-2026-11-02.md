---
name: ppinfo-lever-verdict-2026-11-02
description: "VERDICT on the per-particle information levers at the arm site (band, mask, pose) for the 6.5 % rare state - none the node can pull crosses the bar; band closed by information, core mask = recalibration that gives qtoynull a kept pole, the pose IS the lever but needs ~1-degree rare-particle poses (node has 3.8 median) and stage B is circular; Python port of site_em validated against the Fortran"
metadata:
  node_type: memory
  type: project
  originSessionId: ed42031e-40d4-419c-88a2-6d87bcb479d4
  modified: 2026-11-02T22:18:48.449Z
---

2026-11-02, session tree-ppinfo-lever (worktree suspicious-dewdney-bbfa4e on c5176bddc; scripts/hetero_tree_proto_ppinfo/, note
doc/implementation_notes/tree_ppinfo_lever_2026_11.md; batch2 math_reworks/ppinfo/ runs/<config>_<site>, ~25 core-h):
- ppinfo_site.py = numpy port of site_em + exec_site_rank's acceptance, validated: hold-out curves of the saved posteriors match the
  Fortran point by point (0.149 vs 0.147), EM baseline AUC .768/.792 vs .754/.780, first E-step to 3 digits; n_eff 134 vs 116.
- Ceilings at 6.5 % (box 48, 15 A, site mask): perfect classifier r .35-.37 (passes); perfect references AUC .880, r .13 (fails =
  the EM's own); needed AUC ~.98 (3x per-particle information). Band to 10 A: AUC unchanged, wider hold-out dilutes. Core mask:
  AUC down, z up for nulls too (qtoynull kept pole, 17 seeds). Oracle mask AUC .927, no kept pole. Whole-molecule E-step AUC .946
  but z 2.7 under the site mask.
- TRUE poses: EM kept pole at both seeds (AUC .89-.90, z 3.5-3.7, seed z 3.3); 1-deg jitter just short (z 2.95); 2-deg random
  jitter breaks the EM (45/55 reproducible partitions). Input poses err 3.8 deg median / 9.8 mean; stage B improves the rare
  particles only at 8 % (circular).
**Why:** the coordinator asked whether any per-particle lever lifts the 6.5 % endpoint over the bar with nulls pole-free.
**How to apply:** the lever line is closed for the tree as it stands; the only open direction is a soft pose marginal / local pose
refinement at the site (20-30x the E-step), worth at most 6.5 % -> ~5 %. Never random-jitter poses as a pose-error model.
Related: [[site-climb-lever-no-go-2026-11-02]], [[rare-state-detection-limit-2026-10-29]], [[math-reworks-verdicts-2026-11-02]].
