---
name: mest-formulations-verdict-2026-11-06
description: "VERDICT 2026-11-06: why M-estimation adds nothing to the hetero_tree (one-hot or null-level posteriors; leaf-level weights = frequency-blind pooling that tree_wiener already does per shell); hierarchical/depth tree priors, predictive soft EM and Student-t IRLS all tested in Python on qtoy2 / r131 / nulls / 692k: no state-signal gain, no port"
metadata:
  node_type: memory
  type: project
  originSessionId: 7da7e826-c848-444f-a66e-ad78379bb739
  modified: 2026-11-07T00:00:09.816Z
---

Task from the paper coordinator (2026-11-06): explain why soft per-particle state weights (M-estimation) add so little and test principled
ML formulations (hierarchical tree prior, soft EM under it, Student-t). Python only (scripts/hetero_tree_mest in worktree mest-formulations,
cluster copy math_reworks/mest; runs math_reworks/mest/runs; note side-session-notes mest_formulations_status.md).

Mechanism (mu gate trees' stage-B log-evidence, runs/post/post.txt):
- Toy nodes with real differences: posteriors one-hot (qtoy2 root 77 % >= .99, median gap 11 nats) so M-est = hard labels (purity .675 vs
  .671, soft never purer on any node); the errors are confident (28 % of the >.999 assignments wrong, median 0.000 mass on a right state):
  no reweighting or tempering of the stated posterior can move them.
- Real nodes (r131) and homogeneous nodes: evidence gaps at the null level (median 1.8 vs 1.4-1.6 nats; closed rows' confidence .839 vs
  .825): posteriors soft but uninformative, so soft weights blend states by noise = the 692k closed-arm dilution (anc .91 -> .81).
- Algebra: weights depending only on (own leaf, target leaf) make the weighted map a frequency-blind linear pooling of leaf maps
  (dist_0.4: FSC 6.8 -> 5.0 A, dev .133 -> .121); tree_wiener is the per-shell version (4.6 A, dev .153), hence always better.
  The spliceosome gain came because refine3D_states had no such pooling and flex's weights were genuinely soft (median max .77).

Formulations (all vs the current hard labels + tree_wiener):
- Hierarchical tree prior = tree-covariance Wiener (tree2). Added the tree's own topology, depth-tied and sequential-MAP variants: the
  leaf covariance does not follow depth or topology (R^2 .003-.16 on qtoy2 and 692k); qtoy2 truth dev <= emp2; r131 identical maps and
  corner; 692k nominal FSC 4.73 -> 4.26 A by pooling, closed readouts unchanged to +-.02.
- Soft EM at frozen poses: predictive E-step (reference noise inflation 0.1 % of sigma2) identical to plain; all chains decay below
  production labels; nulls stay 50/50. f2r's latent-span references: same two regimes.
- Student-t IRLS: nu -> Gaussian on qtoy2 and nulls (weights +-1 %, junk AUROC .52 at 10 % junk); nu 30 on r131 with no enrichment change.

**How to apply:** do not revisit soft/M-est weighting, tree priors on the final maps, or frozen-pose soft EM for the paper's "integration"
gap; the only optional port is a depth-tied design in simple_tree_wiener for higher nominal leaf resolution on large trees, which adds no
state information. Related: [[tree-verify-components-and-soft-vs-hard-2026-11-06]], [[flex-refinement-integration-verdicts-2026-11-06]],
[[soft-multistate-fixed-pose-verdict-2026-10-30]], [[spliceosome-mest-refine3dauto-2026-10-02]].
