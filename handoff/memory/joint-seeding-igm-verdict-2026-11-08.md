---
name: joint-seeding-igm-verdict-2026-11-08
description: "VERDICT: joint seeding (latent mixture over tied sites' axes) at the IgM roots: K* 6 on the 5 % null under every Gaussian null (per-particle noise, bin-wise, f2r residual/contrast); groups carry no arm; only the split test's effect floor keeps one group; split runs need ../stageB or the floor is skipped"
metadata:
  type: project
---

2026-11-08, site-chain session for the coordinator (notes: batch2 `hetero_tree_20261028/handoff/site-chain/joint-igm/STATUS.md`).
- Joint seeding = node-levers' prototype: a Gaussian mixture in the latent axes the tied sites' poles use (z1, z3 at every IgM root),
  K by cross-half held-out LL against a null, groups' fixed-pose maps, the root's own `state_split`.
- Three nulls, none null-safe on IgM: A (Gaussian signal + each particle's a Q_i^-1), B (Gaussian per noise-level bin of log det
  a Q_i^-1), C (the f2r null: one Gaussian, deviation scaled by s_i from residual and contrast). K* 6 on the 5 % null and on the
  10 % null forced onto z1, z3; the K 2 gain (.066-.099) is as large on nulls as on fixes. On IgM, residual and contrast explain
  1-7 % of the latent radius, so null C is nearly homoscedastic.
- The groups carry no arm (NMI .000), so joint seeding recovers nothing at the IgM roots (stock root split .034 / .001; chain .022
  on the 5 % fix).
- The tree's split test returns one group only through its effect floor (P / N < 2): the shuffle test passes 350-1,403-particle
  children and the locality guard keeps them.

**Why:** the effect floor reads `../stageB/recvol_stateNN_{even,odd}_unfil.mrc`; the prototype ran `state_split` without it, so
its floor was skipped and its qtoynull PASS is what these roots would deliver.
**How to apply:** any off-tree `state_split` run must put the tested groups' half maps under `../stageB`; do not trust a Gaussian
held-out-gain null for K on IgM-like latents. Related: [[effect-floor-pn2-after-guard-2026-11-02]], [[site-chain-verdict-2026-11-08]].
