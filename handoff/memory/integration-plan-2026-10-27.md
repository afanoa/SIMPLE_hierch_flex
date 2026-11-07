---
name: integration-plan-2026-10-27
description: "User 2026-10-27: integrate ONE production build only after the embedding fix, iteration count and crop timing are confirmed; then combine sequential null + equivalence merge-up + geometry reuse + those three"
metadata:
  type: project
---

User, 2026-10-27 evening: "Wait until the embedding fix, the iteration count and the crop timing are confirmed, then integrate
everything together."
- Ready now: sequential shuffle null (in integration tree), equivalence merge-up (adaptive branch, env SIMPLE_MERGE_RULE=equiv
  test build ~/SIMPLE_smr2), geometry reuse merge_geom_reuse (porting branch claude/tree-port-optim-de0b3e).
- Pending confirmations: embedding fix (flex_embed_fitted_contrast patch; OLE full trees amp_mstep_v1/trees/ole_{base,emb}_s{1,2}_full),
  iteration count ITA=ITB 2 or 3 (amp_mstep_v1 v_it22x trees + start_point_v1 it3 trees), crop timing (porting triple
  ole_nosl_{base,opt,full}).
- Base: committed tip 8c1b66f31; leave the uncommitted design layer out. NOT for production: SIMPLE_SPLIT_CONTRAST=yes (costs
  real RNase P states). Related: [[focus-three-threads-2026-10-27]], [[equivalence-mergeup-verdict-2026-10-27]].
- 2026-10-28 ~00:15: embedding fix CONFIRMED on OLE (no cost: leaves 31/30, 24/22; small <2k 9/10, 10/6; ARI to 6_ .131/.113,
  .108/.074; all OLE trees agree only at particle ARI 0.05-0.16, no 6_ leaf recovered at Jaccard >= 0.5). Crops confirmed (1,039 ->
  374 GB per OLE tree, 0.68x wall). recipe3 = recipe1 (fixm source, ITA=ITB=3) + intopt.diff + embed fix, built on scratch
  /scratch/cluster_scratch/atifao/recipe3 (beegfs Default pool 99.9 % full). Next storage lever: leftover per-part partial
  maps (~half of a cropped tree's bytes).
