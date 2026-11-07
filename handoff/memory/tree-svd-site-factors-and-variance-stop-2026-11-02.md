---
name: tree-svd-site-factors-and-variance-stop-2026-11-02
description: "VERDICTS 2026-11-02 (audit 6.2 / stop rule, branch claude/tree-mask-discovery-stop): site discovery by SVD factors (share >= 0.2, cross-half Fisher z >= 3) committed as default 413c0ce39: nulls unchanged, OLE span unchanged with 16 sites not 22, qtoy2 NMI 0.796 vs 0.764 under variant D, r131 corner neutral over five split-test draws (seed 2: 2/5 vs base 2/5, T_df 0.126 vs 0.117; seed 1: 0/5 vs the rework's 1/5, T_df 0.122 vs 0.121); diffuse-variance stop fires only on null roots (33-46 % vs <= 23 % real), saves ~1 %, kept as a log line only"
metadata:
  type: project
---
2026-11-02, session tree-mask-discovery-stop (worktree tree-mask-discovery-stop-7503a4), note
doc/implementation_notes/hetero_tree_mask_discovery_stop_2026_11.md, data hetero_tree_20261028/mask_disc on batch2 (gates/<arm>_<ds>,
draws/, builds/{base,m23,m23d,final}).
- Item 2 (discover_sites): a component splits into factor sites by the kept factors of a weighted SVD of the even half's deviations
  (Gram K x K), kept at eigenvalue share >= 0.2 and odd-half reproduction z >= 3 over the 15 A voxel count. At r131's pixel size
  n_eff is 2-10, so factors are rarely kept there (10f1/10f2 and 3f1/3f2 fold); qtoy2 keeps its six. Root site order under variant D:
  site 1 first at both r131 seeds, as the reference. Under D the qtoy2 root partition is the reference's to 20 particles; the NMI gain
  is r1's seeding (no kept pole, fflex states seed) plus an r5 draw (over-split of combo 101).
- The r131 corner nodes' production split test is one or two draws in five in EVERY arm (base seed 2 r1 2/5, rework seed 1 r2 1/5,
  m23d 2/5 and 0/5); the draw-free statistic (rework commit 9ee5d5c72, not default) passes all of them at z >= 6. Whether a given
  tree shows the corner at ~20k-particle nodes is the draw, not the site rule.
- Item 3: the stage-B variance level set holding 0.7 covers 33-46 % of the molecule on nulls, <= 23 % on real splits; a leaf rule at
  30 % skips only the split job (0.15-0.72 core-h per null root) and never fires on real nodes. Logged at every split step instead
  (variance_spread); the stop is the one-line frac > 0.3 branch if wanted.
**How to apply:** judge r131 md1 corner verdicts by several split-test draws (math_reworks/build_df + SIMPLE_SPLIT_SEED, draws.sbatch)
or the draw-free T_df, never by one tree; do not propose variance-diffuseness stops for the 692k deep lineages (they are localized).
Related: [[tree-mask-share-core-no-go-2026-11-02]], [[split-test-single-draw-variance-2026-11-02]],
[[site-strength-ranking-decides-r131-corner-2026-11-01]], [[tree-design-feature-freeze-2026-10-25]].
