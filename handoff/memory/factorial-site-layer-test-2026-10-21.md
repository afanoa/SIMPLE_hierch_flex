---
name: factorial-site-layer-test-2026-10-21
description: "factorial site layer, DONE 2026-10-22: built INTO the tree (tree_sites.py post-pass; user likes the hierarchical tree): per node discovery -> focused flex per site + two-tailed site EM at the parent's shared poses, pooled up by a ranking gate; RNase P closed corner found automatically; closed pole reproduces across runs / halves; site EM mask-robust; discovery must not merge sibling factors and must rank by variance share"
metadata:
  node_type: memory
  type: project
  originSessionId: 1ca67e82-4d15-4725-b95d-0ee5719be9f4
  modified: 2026-10-22T08:58:14.420Z
---

Results doc: `doc/implementation_notes/factorial_site_layer_results_2026_10_21.md` (worktree laughing-carson-4e8ca1, uncommitted; 45
scripts in its `scripts/cls_expansion_proto/toy55/`, synced to batch2 routes/toy55; outputs rnasep_klem/factorial).
- Built: `tree_sites.py`/`.sbatch` = site layer inside a finished tree (read-only on the tree): per node discover_sites2 (local <= 30 %),
  per site a focused flex (pcg_mskfile = site mask) + 2-class site EM from both tails of the site flex z1 at the node's n.simple poses;
  pool to the parent while the parent EM ranks the node's particles at Spearman >= 0.5 (RNase P closed site 0.07-0.10 stays in r1;
  IgM arms 0.51-0.91 and the simsites arm 0.66 pool).
- Results: simsites tree+sites domain 0.44 / arm 0.64 eta2, 3/3 reachable combos (tree 0.36 / 0.63, 2); RNase P r1 closed pole
  1.07 / 0.49 at top-2.5k automatically; IgM every arm found, per-site flex both angles at the oracle level, but the product of many
  site labels collapses to the leaves -> output must stay factorial.
- TRAPS fixed: Dice dedup merged sibling FACTORS (lost r1's closed factor 1f1 and IgM arms) -> discover_sites2 keeps them apart;
  ranking by reproducibility picks the wrong factor (port runs) -> rank by the share of the node's reproducible variance.
- Overnight (coordinating session, plan of record `doc/implementation_notes/state_tree_reproducibility_overnight_2026_10_22.md` in
  /Users/atifao/SIMPLE): pole SET {closed, open A, open B} reproduces across 3 runs at a shared mask (Spearman 0.96-1.00); masks from raw
  latent-axis tails do NOT; the site EM is mask-robust (Dice >= 0.89 -> Spearman >= 0.997) and discovery on recB re-derives 1f1 in
  every run (Dice 0.92-0.94); disjoint halves reproduce sites (Dice 0.94-0.98) and poles (0.91-1.02 of the half-map noise ceiling);
  fflex inputs partly, cheap children never give the closed site.
- TRAP: top-N map readouts biased high at small N; compare at matched N.
**How to apply:** keep the tree as the backbone; site masks from node-level discovery on refined states, all local sites fitted,
poles from two-tailed starts, classify at the parent's shared poses, polish poses per class, factorial output. Related:
[[klem-integrated-hierarchical-model-vision]], [[tree-leaves-not-reproducible-across-runs-2026-10-21]], [[no-region-basis-for-rnasep]].
