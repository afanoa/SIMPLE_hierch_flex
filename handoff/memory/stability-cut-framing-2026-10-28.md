---
name: stability-cut-framing-2026-10-28
description: "user 2026-10-28: agglomeration alone not worth pursuing; frame states as a partition between all-singletons and one-state (graph, 'bidirectional'); adopted framing = per-branch stability (excess-of-mass) cut over the top-down tree, replay queued"
metadata:
  type: feedback
---

2026-10-28 ~21:00, after the bottom-up replay merged almost nothing (0-1 merges per tree by equivalence, 0-5 with locality):
"dont think agglomeration is worth pursuing on its own. How would you mathematically frame this problem. Either every particle is in
its own state or every particle is in one state. Our solution lies in between. Surely there is some graph ... sort of like a
bidirection djikstra".

**Why:** at large N every significance test says "different": top-down never stops, bottom-up never merges. The criterion must be
persistence / relative effect size, not p-values.
**How to apply:** partition lattice; MST / single-linkage hierarchy (Prim ~ Dijkstra); pick per-branch levels by stability
(HDBSCAN excess of mass): keep a node if n*(s_birth - s_death) >= its descendants' selected stability, s = split effect size on a
common scale (SNR_d / locality-weighted rho). Replay by adaptive-split-merge on RNase P closed lineage, qtoy2/IgM truth, OLE pairs.
Archetypal analysis = the "extremes only" variant. Related: [[close-closed-substates-not-worth-a-split-2026-10-28]],
[[landscape-must-be-reproducible-and-splits-data-derived]].
- RESULT 2026-10-28 ~21:50 (adaptive replays, route9/stab, /repro, /agg): NO consistent knob-free consolidation. Stability cut is
  scale-dependent (units per scale rho / SNR_d / sqrt(rho): OLE ref6 31 -> 14 / 2 / 31; IgM m52M 89 -> 77 / 21 / 89); sqrt(rho)
  best on qtoy2 truth (F1 .574 -> .616) but = leaves elsewhere; SNR_d over-merges (OLE 2 units, IgM extremes lost). Reproducibility
  merge: split-half null BC 0.994 +- 0.001 too tight (OLE 31 -> 29, strong 0 merges). Bottom-up equiv+locality: 2-5 merges/tree.
  RNase P: no rule consolidates closed into one unit. OPEN: IgM extremes 14/20 (m52M) and 9/20 (m69M) at leaves vs 0/0 at depth 3-4
  cuts per adaptive's scorer, but factorial's igm_score.py says 0/20 at m69M leaves -> definitions being reconciled (maxdepth risk).
