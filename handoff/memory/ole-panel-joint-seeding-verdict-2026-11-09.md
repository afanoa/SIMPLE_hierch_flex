---
name: ole-panel-joint-seeding-verdict-2026-11-09
description: "VERDICT: OLE on merge_master2 + October keys (4 maxdepth-1 trees): the 4-way root diverges in part (split .65-.77, XD K 4/6/5/5); depth-1 nodes keep no site pole and seed deterministically (same-rows replays 1.00), so depth-1 divergence is inherited from the root; null-safe joint seeding does not help (root .34-.67 < single-site .49-.95; depth 1 refined .76 = single-site .77); TRAP: coordinates file numbers particles 1..n at sub-nodes"
metadata:
  node_type: memory
  type: project
  originSessionId: 6f195b76-225b-466c-a432-7fc8360148e8
  modified: 2026-11-09T06:25:39.914Z
---

Side session sparse-ppca2, parallel slate ([[parallel-slate-2026-11-08]]). Notes: `hetero_tree_20261028/handoff/ole-panel/OLE_PANEL.md`
(shared notes `ole_panel.md`); runs and scripts in `hetero_tree_20261028/ole_panel/`. Total 244 core-h.

- **Trees:** merge_master2, ole_r14 keys, maxdepth 1, nodepar 1, KEEPMAPS, seeds 1-4. About 2 h 27 min and 42 core-h per tree.
- **Root:** splits 4 ways in every seed but reproduces only in part.
  - Split ARI .65 for s1 against the rest, .73-.77 among s2-s4. Seeds .49-.53 vs .86-.95.
  - XD K flips 4/6/5/5 (gaps <= .0007), which reshapes the factor sites and s1's site_states pole seeds.
  - Agreeing seeds (.95) lose a fifth in stage A/B.
  - s1 and s4 have identical fflex coordinates (r 1.00) yet seeds .53, through the XD K flip alone.
- **Depth 1:** no node keeps a site pole (16 of 16), so every node seeds from fflex's fallback.
  - On identical rows, replays at SIMPLE_SEED 2 give XD K, fflex states, seeds and joint groups at ARI 1.00. Refinement alone costs .16-.30 (split .70-.84).
  - So the divergence across trees (rows J .70-.86) is inherited from the root. October's mechanism (depth-1 site / pole near-ties with a reproducing root) does not reproduce.
- **Joint seeding** (f2r's null-safe K, [[near-tie-is-duplicate-factor-masks-2026-11-08]]) does not fix OLE.
  - Root: K* 5/4/6/6, and the groups agree at .34-.67 against .49-.95 for single-site seeds. A raw-latent variant is worse (.15-.55).
  - Depth 1: refined joint children .76 vs single-site .77 on the same rows, and .20 vs .21 across trees. Same children kept, except r2 (3 vs 2).
- **Null-safety:** my exact-noise Gaussian null failed all 3 nulls (floor ~0; real null-fixture gains .001-.005). f2r's residual/contrast-scaled null passes olenull/qtoynull. micnull K* 3 stays.
- **TRAP:** at sub-nodes (whose rows are not 1..n), `fflex/flex_pca_coordinates.txt` numbers particles 1..n, node-local and in node order, while `rows.txt` holds project ids. Never map the coordinates' particle column through project ids.

**Why:** OLE is the panel's reproducibility dataset ([[methods-must-generalize-across-datasets-2026-11-08]]); this tells whether the RNase P r1 fix generalizes.
**How to apply:** for OLE, look at the root (XD K's flat curve, site_states seeds), not depth-1 seeding. Judge depth-1 methods on same-rows replays: independent trees confound them with the root. Extends [[ole-repro-seeding-near-ties-2026-10-28]].
