---
name: extremes-fission-verdict-2026-11-07
description: "VERDICT 2026-11-07: self-selected tail/leaf maps place extremes 1.3-3.1x too far out and fake maps 'reproduce' (.16-.45); fission re-reads fix it; the r131 closed-corner HARD GATE passes only on own-data imprint (corner leaf 1.03/.54 own, .93/.59 honest; no leaf passes honestly); global coordinate / node latent find extremes no better than leaves; focused coordinates on data-only site masks do (IgM)"
metadata:
  node_type: memory
  type: project
  originSessionId: cc23038c-fea9-4b5f-b321-d58eb2942350
  modified: 2026-11-07T20:53:54.099Z
---

extremes-repro session (user asked "Run it": can the latent cut out extremes that reproduce?). Write-up
`hetero_tree_20261028/handoff/extremes/{VERDICT,METHOD,RUNS}.md`, tables/<node>/; runs `math_reworks/latent/extremes/`.
- Self-selection is large: at the same particles, maps built from the data that selected them place extremes 1.3-3.1x further
  along the true motion (IgM arms 1.74-2.5x even for 90-98 %-pure tails; leaves 1.75-1.8x); tails picked along a fake map no
  particle carries "reproduce" across halves at .16-.45, honestly at 0. Fission (select on X1, rebuild from X2) leaves <= 1.8 %
  (qtoyflex) and <= 3.5 % (r131 registered stacks, band-wise radial noise model) of the imprint; honest maps read the truth on one
  scale (IgM a1 beta/shift 345-360 at every tail size vs 613-793 self-selected).
- r131 HARD GATE (closed >= .95, open <= .6): the delivered corner leaf r2_s2 reads 1.03/.54 on its own data, .93/.59 honestly
  (random 5 % tails .87-.90, node mean .89-.90/.59). No leaf passes honestly; the gate as defined passes on imprint. Caveat: the
  honest leaves are a re-assignment among the tree's leaf references on X1 (agreement with delivered labels .50, chance .25), not a
  tree run on X1. Only small latent tails reach the threshold honestly (seam-free pc2 bottom 5 % .966/.574; site-2 focus top 5 %
  .953/.579; ~2.2k particles).
- Global coordinate and whole-node latent find extremes no better than the leaves (tie on r131, worse on qtoyflex: closed .197
  for the tree vs <= .13). Focused coordinates on a correct data-only site mask win where the tree is blind (IgM 10 %: a1 tails
  90-98 % pure vs 62 % for leaves; an a2 tail 87 % in its -30 deg state with no a2 leaf; null extremes +-12-15 deg vs +-1.7).
- Traps: r131's registered stacks have non-stationary noise (corners .36-.66 of the inside low/mid noise); the map-level
  even/odd cross-copy check cannot validate fission on real data, only the fake-map null can; the CTF-flip trap
  ([[rnasep-denoised-stacks-are-phase-flipped]]) was found here.

**Why:** answers the user's "cut out the latent so the extremes can be reliably reproduced" and questions how the closed-state
gate is read.
**How to apply:** any own-data map readout (tails, leaves, the r131 corner gate) needs a fission re-read before it counts; an
honest gate needs fission readouts and a threshold calibrated on them (decision for the user). Related:
[[rnasep-closed-state-must-be-retained]], [[cross-half-test-needs-band-above-assignment-2026-11-07]], [[sparse-ppca-pursuit-verdict-2026-11-07]].
- DECIDED 2026-11-08 (user: "Yes, honest re-read", decision 9): the r131 corner gate is read on a held-out fission copy.
  - The gate tree runs on X1; the corner leaf is chosen on its X1 map and read on X2.
  - Bar: the 95th percentile of matched-N random tails and random subsets, read on X2.
  - extremes-repro builds it (guide 7.10 thread 8) and applies it first to the noise-model seed table.
- 2026-11-08 afternoon, honest gate seeds 1-2 (guide 7.10 thread 8; `$H/handoff/extremes/GATE.md`):
  - Results: pp_off .893 vs bar .892 (PASS) and .870 vs .883 (FAIL); noise model .843 vs .856 (FAIL) and .864 vs .855 (PASS).
    The arms split one-all with narrow passes.
  - Every corner loses .04-.18 of closed occupancy from X1 to X2.
  - The stricter null is the latent tail along the parent's own X1 latent.
  - Root splits reproduce within an arm (77-79 %) but not across arms (50-56 %).
  - Seeds 3-4 (55084993-96) decide the noise default.
- 2026-11-08, tail-seed refinement rounds at r131 (thread 9):
  - Every round ends near a half split (~19k of 44k) at about the base closed share, whatever the seed.
  - Random seeds clear the bar 3 / 4, latent tails 1 / 3. rsub2 (both halves at the base share) still clears it.
  - IgM a1: a 97.5 %-pure 5 % tail grew to a 44 % state at 64.5 % purity.
  - The production fflex round (94 % = the tree's root split) is the most closed-enriched, at .085.
  - Rule: clearing the honest bar after a round shows non-random selection, not closed enrichment. Judge rounds by closed share
    and recall.
