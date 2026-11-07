---
name: joint-gate-c7ac3b0bd-2026-11-02
description: "joint gate replay of c7ac3b0bd (step-1 cleanup + subspace site test + SVD discovery + draw-free split test) = math_reworks/gates/head_*: three nulls pass (one leaf), r131 closed corner kept on both seeds, qtoy2 NMI 0.773 (10 leaves); the bc arm (+ seed-anomaly fix + P/N >= 2 effect floor) 0.812 with 8 leaves = the 8 true states; OLE 8 leaves = the accepted df arm's tree (bc identical but 95 particles at r3); wall time unchanged"
metadata:
  type: project
---

2026-11-02. build_head (= c7ac3b0bd, built by the tree-lever-floor session on batch2 at hetero_tree_20261028/math_reworks/build_head,
NOT by the user) and build_bc (head + the two uncommitted lever-floor changes); readouts gates/{head,bc}_readout.txt and
{head,bc}_s1_readout.txt (closed_nodes2 views in gates/s1_readout and gates/bc_s1_readout).
- nulls qtoynull / olenull / micnull: one leaf (r_g1_2) in both arms = PASS.
- r131 seed 2: root 25,397 / 18,937; corner = leaf r1_s1 (node r1 state 1, 1.00 closed / 0.49 open), both arms identical trees.
- r131 seed 1: root 23,449 / 20,885; corner = leaf r2_s2 (node r2 state 2, 0.97 / 0.53), both arms identical.
- qtoy2: head 10 leaves NMI 0.773 (r3 and r5 split within a conformer; combo 28 F1 .64, 101 .71); bc 8 leaves NMI 0.812
  (floor merged r3 and r5; every combo F1 .73-.95); previous df arm 9 leaves 0.742; quality-axis NMI 0.000 (genuine).
- OLE: head and bc 8 leaves (4 root groups x 2) = the df / sub4 arms' tree to a few hundred particles (base had 10); bc = head but 95 at r3.
- OLE landscape (landscape_score on the base root embedding, gates/ole_headbc.txt): span ends identical in every arm (a1 .48-.49, a2 .56, a3 .49), OT W1 head/bc vs base 0.002 (sub4/df 0.003), head vs bc 0.000; partition coverage vs sub4/df 0.44-0.54 = the null's level, so the leaf partitions differ at the null level while the landscape is the same.
- Wall minutes head vs df: nulls 25-35 vs 28-38, qtoy2 51 vs 49, r131 55-58 vs 56-59.
**How to apply:** c7ac3b0bd is gate-clean on 6/6 fixtures; the bc changes (seed fix + floor) improve qtoy2 and are neutral elsewhere, committed
as cb3196f0a / 577ee25fc on claude/tree-lever-floor-bccc29 (with D step 1 fe02cffb8 and E d062ba9f2; HEAD fast-gated as build_lf). Related:
[[math-reworks-verdicts-2026-11-02]], [[site-climb-lever-no-go-2026-11-02]], [[rare-state-detection-limit-2026-10-29]].
