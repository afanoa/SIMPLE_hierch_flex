---
name: near-tie-rule-verdict-2026-11-09
description: "VERDICT: the 2 SE near-tie rule does not make the tree's top choices agree (XD K OLE root 4/5/5/4, RNase P root 8/7/9; site region already agrees, f1/f2 tie-break is itself a near-tie); with four XD starts + 2 SE, K is 6/6/6 at the RNase P root but not OLE (3/5/5/4); exact numpy replica of xd_select_k exists"
metadata:
  node_type: memory
  type: project
  originSessionId: 6f195b76-225b-466c-a432-7fc8360148e8
  modified: 2026-11-09T09:29:45.914Z
---

Overnight task (coordinator, 2026-11-09). Notes: `hetero_tree_20261028/handoff/ole-panel/NEAR_TIES.md`; scripts in `ole_panel/scripts/`.

- **`stagea_xd.py`** is an exact replica of flex_pca's `xd_select_k`. It reproduces every node's logged ladder to 6 decimals (24 nodes) and keeps the per-particle held-out log-likelihoods.
- **K rule** (smallest K within 2 SE, paired per-particle SE): OLE root 4/6/5/5 -> 4/5/5/4; RNase P root 8/11/11 -> 8/7/9; RNase P r1 already 5 everywhere; qtoynull 3 -> 1.
  - The boundary K sits at 2.2-3.0 SE, and the curves are non-monotone from EM optima.
- **Four starts per K** (`stagea_xd2.py`: quantile starts on the top 3 variance axes plus a warm split of the K-1 fit, best training LL):
  - RNase P root 6/6/6, which agrees;
  - OLE root 3/5/5/4, since the seeds' raw embeddings differ.
- **Site choice:** margins are z-scores with unit SE, so the tie band is 2 sqrt 2.
  - The chosen region already agrees: OLE 6/6, RNase P root 3/3, r1 28/28 replay pairs. The tie is between one region's f1/f2 masks.
  - The replication tie-break decides on reliability differences of .005-.01.
- **Most top-level disagreement is outside the choice:** the seed set (site_states pole seeds follow XD K) and root refinement (RNase P root same site: seeds .76 -> split .37).

**Why:** the user's principle "never choose on a difference below its own noise" ([[landscape-must-be-reproducible-and-splits-data-derived]]).
**How to apply:** a reproducible K needs a better-optimized XD (several starts) as well as the 2 SE rule, and only where the embedding reproduces. Do not expect a choice rule to fix the roots' seed sets or refinement. Related: [[ole-panel-joint-seeding-verdict-2026-11-09]], [[two-run-consensus-verdict-2026-11-09]].
