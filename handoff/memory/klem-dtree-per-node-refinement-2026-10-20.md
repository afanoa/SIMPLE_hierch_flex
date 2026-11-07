---
name: klem-dtree-per-node-refinement-2026-10-20
description: "user 2026-10-20: refine poses at every tree node (\"what's the point of the tree if we didn't refine at each node\"); dtree_v13 + node_refine.sbatch = f692e family recipe per node and half; frequency march kept; timings; cheap variant under test"
metadata:
  type: project
---

Context: the gold-standard trees (dtree_v10-v12) ran every node at the CONSENSUS poses; f692e had refined each family
(refine3D_states, k-means K=10 on the child board, local 20->8 A, maxits 6, ml_reg=no). User: per-node refinement,
same recipe as Oct 18, keep the frequency march, run in parallel; expects the heterogeneity model to work better with
refined poses and refinement can be cheap.
- toy55/dtree_v13.sbatch: from depth REFDEPTH (1), after the node's boards, node_refine.sbatch per half (k-means REFK
  states, 4-10 by size, on the node board's posterior means -> refine3D_states pose_policy=local lpstart 20/14/12 by
  depth -> lpstop 8, maxits MAXITS 6, mskdiam 220, REFINE_EXTRA ml_reg=no); the node's split test (splitbench_node6
  POSEA/POSEB) and its children's boards (dtree_board5 PROJ) use the refined poses; RESUME_ROOT reuses a finished
  root (root stays at consensus poses, like f692e's parent board). refine3D_states nstates=3 crashes, so K >= 4.
  tree_finalize4 (merge_leaf_poses.py: each leaf takes its deepest refinement) + state_maps_readout.sbatch (bench maps,
  frame check 0.95 identity, closed-arm readout validated on sixc8: closed state 0.55 vs 0.12-0.36).
- Timings: real depth-1 nodes (65-163k per half, 10 states, 8-10 parts) 28-29 min each regardless of size, poses moved
  2.4-3.2 deg median (p90 5-6.5); synthetic depth 1 (6-16k, 4 states) 10-13 min, ~2 deg; synthetic depth 2 ~7 min,
  median 0 deg (p90 2.9): a second refinement barely moves poses. Cost = per-state master work, not particles.
- Runs (2026-10-20 ~00:00-00:55): real gs2 (vs gs1 unrefined), synthetic sixc8r and sixc2r (vs sixc8/sixc2), cheap
  variant sixc8c (REFK=4 MAXITS=5 REFINE_EXTRA=filt_mode=none); f692e validation also with f692e family poses
  (f692e_half_poses.py, f692e_validate2.sbatch, xem_prot_famposes / xem_noprot_famposes). Chains: refined_chain.sh.
- RESULT sixc8r vs sixc8 (same root, 700k six + closed 18 at 8 %): final states 6 vs 5 (47 now separate from 49,
  0.54/0.81), NMI 0.325/0.347 vs 0.292/0.265, closed state 18 purity 0.53/0.47 vs 0.47/0.38, closed-arm 0.60 vs 0.55;
  the leaves are the same (NMI 0.226/0.210 both), so the gain is in the final pass (refined poses in E-step + maps).
- POSE ACCURACY vs truth (pose_truth_check.py; sim all.txt e1-e3, true shifts are 0 because simulate_particles with
  oritab never applies sherr): consensus 4.1 deg median rotation error (p90 7.6), merged refined 3.8 (p90 7.3); 38 % of
  particles better, 28 % worse; shift |error| 1.46 px both. Refinement moves poses ~2 deg but gains only 0.3 deg.
- POSE CEILING (tree_finalize5.sbatch POSEA/POSEB: the SAME sixc8r leaves through the final pass): consensus 6 states
  NMI 0.30/0.30, closed 18 purity 0.49/0.40, closed-arm 0.56; refined 6, 0.33/0.35, 0.53/0.47, 0.60; TRUE poses 5,
  0.45/0.42, 0.73/0.65, 0.74. Poses ARE a major limit; per-node refinement recovers ~1/4 of the gap. Next: oracle-label
  refinement (oracle_refine.sbatch, true conformers as states, same recipe) separates reference purity from noise.
- CHEAP vs FULL (sixc8c REFK 4, MAXITS 5, filt_mode=none): 6 states, NMI 0.31/0.32, closed 0.51/0.45, arm 0.60 = about
  half the full recipe's gain (0.33/0.35). 2 % closed (sixc2r, refined): 6 states, NMI 0.33/0.33 (unrefined 0.31/0.32),
  still NO closed state (18 <= 5 % of any state); leaves identical to the unrefined tree again.
- REAL gs1 (unrefined, consensus poses) final: 28 leaves -> 14 after merges -> 9 states after the cross-half EM; no closed
  map (max closed-arm 0.56, 4 intermediate); known-closed proxy (f692c zoom 1:3) enriched ~3x in branch r.2 (3.9 %,
  leaves to 6.8 %) and in final state 2 (4.7 % / 5.0 %, 55 % of all proxies, 61k particles) -> enriched, not isolated.
  NOTE: closed_proxy_track per-node lines start with two spaces; don't filter "^  [A-Za-z]" (banner filter) on them.
- ORACLE / POLISH pose accuracy (half A; half B lost to the corrupted project, rerun 02:15): oracle (true conformers as
  refine3D_states labels, all 49k) 3.53 deg median, p90 6.8; POLISH (the final pass's 6 states as labels, all 49k) 3.53
  = oracle; per-node merged 3.84; consensus 4.13; shift error unchanged (1.47) in all. -> refining the WHOLE half with
  the final states beats per-node refinement (per-node states have 1.5-4k particles) and reaches the label ceiling;
  the rest of the gap to truth is noise/recipe, not labels. Truth poses in the final pass: NMI 0.45 vs 0.33.
- REFINE3D_STATES RE-CLASSIFICATION is the separating step (user's manual-workflow observation): from the TRUE labels
  (NMI 1.0) it ends at 0.52 (25 % moved); seeded with the final pass's 6 states (NMI 0.33) it ends at 0.46 (half A,
  34 % moved) = the true-pose final-pass ceiling (0.45). It converges to its own ~0.5 partition whatever the seed.
  Per half -> gold-standard; K fixed by the seed. Next: final pass with these states as candidates (tree_finalize6
  CANDA/CANDB, run polishc). TRAP hit: a failed oracle_refine left poses_full.txt -> symlink to the INPUT table; the
  rerun wrote THROUGH it and overwrote lat_55_700k_B/ptcl3D_r0.txt (restored 02:45 from final_xem_cons/poses_B.txt,
  md5 verified); oracle_refine now rm -f's its outputs first (node_refine* still have the hazard on reruns).
- REAL gs2 (per-node refinement, v13) final: 12 -> 8 states after the EM (3 drained, two strongly half-asymmetric);
  known-closed proxy: 70 % / 67 % of all proxies in state 2 (5.4 % / 5.8 % of 68k) vs gs1 55 % (4.7 / 5.0 %); best
  leaves r.2.1.2.2 7.5 % / 7.4 % and r.2.2 children 7.7 / 6.9 % (gs1 max 6.8) -> refinement concentrates the closed
  particles slightly more, still not isolated; max closed-arm 0.55; one depth-2 refinement lost (r.3.4 A, corrupted project).
Related: [[f692e-validation-cross-half-2026-10-19]], [[klem-state-graph-umap-and-real-halves-2026-10-19]].
