---
name: inode-prototype-2026-11-06
description: "integrated hetero_tree node prototype (mixture of PPCA with latent poses, cross-validated halves) for the paper; paths, fixtures, first toy verdicts"
metadata:
  node_type: memory
  type: project
  originSessionId: dc17548d-1b94-4fe8-ad02-93f812b3f2a7
  modified: 2026-11-06T20:44:04.971Z
---

User 2026-11-06: "Put the integrated node model in this paper: prototype it, then run the current tree as the baseline it is
compared against." Prototype = batch2 /mnt/beegfs/elmlund/afan/paper_inode (scripts/inode_em.py v1, inode_em2.py = + --env
--ctfabs --retest and a nested held-out test; make_fix.py; score_tree.py; export_node.py for real nodes; logs/, runs/, trees/).
Built on the toyJ EM (~/toyJ/toyj_em.py, frozen copy of toyj_core.py): K states x r components, z-marginal (k, pose) posterior,
rigid tangents projected off U_k, half h scored against half 1-h's maps (swap) / own (indep) / one model (merged = stage-B analog).

Fixtures from toy3 (t3s0.02, 39,990 ptcls, +-15 deg hinge, 31 levels) and toyPS t3null: f50 (<= -6 / >= +6 deg, 12,900 each),
f15 (2,276 / 12,900), cont (26k all levels), null (26k). Baseline = build_mu (production 72e5aaba0 source) root-only trees
(trees/<fix>, keys maxdepth 0, site_layer/site_states yes, mskdiam 200, minpop 500).

BASELINE production tree (root only): f15 2 leaves NMI .967, f50 .978, cont split at the median, null 1 leaf (PASS); root stage-B
pose scatter 1.85 deg. The production node BEATS the bare prototype on classification.
First verdicts (seed 1): classification integrated ~= stage-B analog (f15 NMI .91-.93 vs .92-.94) > fixed poses .86, oracle 1.0;
integrated keeps pose scatter at input (1.48-1.6 deg) where greedy stage B inflates it (1.77); NO arm reduces the within-state
conformation-coupled pose bias (~1.4 deg median, 2.4 max). HELD-OUT two-vs-one likelihood is NOT a calibrated split test: clean
design (indep training, per-half refits, rival nested in the K-state span; runs/rt3_*) gives f15 1.87 (minority particles -0.89,
majority +2.35), cont 1.78, null 0.58 nats/ptcl; swap training inflates (null 0.70). Keep the tree's draw-free test for decisions.
RNase P r1 (tree seeds, ctfabs, contrast, box 64 lp 10, whole-image scoring; runs/r1_*): closed recall / purity of the best state:
integrated .757 / .041, stage-B analog .630 / .044, fixed .594 / .038 vs the PRODUCTION tree's r1 refinement .889 / .063 (r13).
Seeds alone: .51 / .013 (closed mostly in the rest seed). FOCUSED round (inode_em6.py, --focus node mask, --stages 15:3/10:3 —
NEVER a comma in ARGS: sbatch --export cut the first round's 15:3,10:3 and dropped the seeds, [[sbatch-export-splits-at-commas]]):
the cut round (lp 15 only) already gave r1 closed .685 / .084 in BOTH arms (focus, not integration, lifts purity; no rest collapse).
VERDICT corrected focused round (v6: tree seeds, node-mask focus, lp 15:3/10:3, contrast on r1; runs/*_int6, *_stb6): r1 closed share /
purity integrated .831 / .062, stage-B analog .849 / .058, production tree .889 / .063 -> the prototype reproduces the production node
but INTEGRATION ADDS NOTHING over the stage-B analog; toys f15 NMI .929 vs .911, f50 .936 vs .948 (tree .967 / .978); null: stage-B
analog empties the 2nd state, integrated keeps 673; pose scatter 1.79 vs 1.91 deg. Focus is the lever, not integration. All arms collapse the 135k rest seed to 0.3-2k; posteriors soft (>0.99 for
2-5 %). Likely gap: the tree focuses references on the node mask and refines at lp 15 -> 10, box 64 -> 88; the prototype does not.

**Why:** the paper's integration claim must survive the null and beat the production node, not a strawman.
**How to apply:** read runs/rt_* (retests) and logs/inode_score.* (baseline trees) before claiming any gain; frame-free pose metrics only
(state maps absorb a coherent pose bias as a frame rotation, [[pose-state-coupling-toy-verdict]]). Next: real node RNase P r13
(46,991 ptcls, f692_lf_s1) via export_node.py with --ctfabs 1 (phase-flipped stacks, [[rnasep-klem-boards-need-abs-ctf]]) and
stkind/indstk lookup ([[rnasep-692k-project-rows-are-a-permutation-of-stack-order]]). Related: [[hetero-tree-paper-plan-2026-11-06]].
