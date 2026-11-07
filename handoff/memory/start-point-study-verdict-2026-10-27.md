---
name: start-point-study-verdict-2026-10-27
description: "Start-point study closed (545 core-h): annealing (state_anneal) and majority consensus fail; 'keep the most-confident stage-B restart' picks the truth-best basin (rank corr .88-.93) but is NOT integrable until the split test rejects the quality null"
metadata:
  type: project
---

2026-10-27, session start-point-sensitivity (branch claude/start-point-sensitivity-cb0d4a, uncommitted; cluster
hetero_tree_regression_20261026/start_point_v1; write-up doc/implementation_notes/planned/hetero_tree_start_point_sensitivity_2026_10_27.md there).
- Most-confident restart (mean max state posterior of the STAGE-B partition, prob_state_evidence.txt): rank corr with truth
  r1211 .88 (23 restarts), CAHRA .93 (22), qtoy2 .89 (13); top restart always in the truth-best basin; ties inside a basin;
  fails on stage-A partitions (.12) and annealed restarts (.16). Fails the null: quality splits on the one-conformer null are
  MORE confident (.86-.90) than real data.
- state_anneal (opt-in refine3D key): does not cross r1211's basin, costs the good start 7-18 points precision -> don't commit
  (diagnostic only). Majority-vote consensus picks the majority (wrong) basin. Not non-convergence (9-it stage B converges at 4).
- Start-sensitive nodes: RNase P r1211/r12111, OLE r111/r1111; barely: qtoy2, CAHRA, OLE r11111.
**How to apply:** revisit restart selection only after the split test has a quality-aware null
([[split-test-accepts-quality-split-on-null-2026-10-27]]); cost is R x node refinement.
- Pass count (r1211 subtree, 20:40): fixed site-pole start closed P/R .75/.53 at 1+1 vs .98/.98 at 3+3; plain trees (site_layer=no)
  1+1 scatters closed over 3 leaves (best .691/.334), 3+3 keeps one leaf .474/.914 but cannot isolate it without the closed-pole start.
