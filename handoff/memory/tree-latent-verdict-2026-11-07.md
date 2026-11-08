---
name: tree-latent-verdict-2026-11-07
description: "VERDICT: shared sparse-PPCA site basis + XD tree + sat5 fission splits works on fixtures (IgM10 .322 with K=3 vs lever 3 .338; qtoyflex/qtoy2 beat the tree; nulls 1 leaf) but on RNase P only finds misfit; port on hold until a per-particle misfit model"
metadata:
  node_type: memory
  type: project
  originSessionId: 1218bb2b-a21a-4d25-9e48-c40ddd89ceff
  modified: 2026-11-08T05:19:30.142Z
---

The tree-latent session (2026-11-07/08) built the cohesive tree in latent space. Write-ups are on batch2 in
`handoff/tree-latent/` (VERDICT.md, STATUS.md addenda 1-3, PORT_DESIGN.md); runs are in `math_reworks/tree_latent/`.

The model: one site basis shared tree-wide, per-particle deconvolution (posterior precision A), K=2 plus K=3 XD proposals, and a
split accepted only on a sat5 DISCRETE call (pairwise for K=3), in the forest form (one tree per site), with Sidak bars.

| Fixture | Latent tree | Comparators |
| --- | --- | --- |
| IgM 10 % fix | five-arm product .322 with K=3 (.258 with K=2); a2 resolved (.298), a5 half (.216) | hetero tree .042, lever 3 .338 |
| qtoyflex | closed F1 .558, holo .713 | hetero tree .39-.47 / .62 |
| qtoy2 | closed .91, holo .96; conformer 1 \| 2 pair still merged | pp_on .891 / .944 |
| Nulls | one leaf everywhere; closest K=3 null side margin +.04 at k 2 | |
| RNase P | r131: one leaf, gate fails. 692k: the accepted split tracks likelihood / amplitude, closed-depleted child | |

M-estimation: hard maps by default. Soft dilutes contrast. Unmixing via the model's own confusion matrix wins on IgM but trades
bias for noise on small leaves.

**Why:** the user wants "a cohesive model over the entire tree". The prototype shows that the model works on fixtures and that
real data is blocked by per-particle misfit, not by the tree ([[real-data-fission-needs-misfit-model-2026-11-08]]).

**How to apply:** do not port it to Fortran until a per-particle band-wise noise model exists (rnasep-misfit session). Keep
memberships on y1 ([[fission-select-members-on-y1-labels-may-read-y]]). Related:
[[data-fission-merge-test-verdict-2026-11-07]].
