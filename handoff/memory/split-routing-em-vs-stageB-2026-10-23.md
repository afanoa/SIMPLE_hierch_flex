---
name: split-routing-em-vs-stageb-2026-10-23
description: "CLOSED 2026-10-25, EM routing stays default; hetero_tree children follow the split test's cross-half EM labels (assigned on <= 20 A shells): that loses HIGH-res states (CAHRA root open: precision 0.29 vs 0.685 by stage-B state) but beats stage-B routing on LOW-res junk splits (toys); no fixed routing rule works"
metadata:
  node_type: memory
  type: project
  originSessionId: 58de2c8f-ad53-44fb-b8a2-05aadff16ec4
  modified: 2026-10-25T16:56:20.769Z
---

Measured 2026-10-23 for the coordinator, on the CAHRA m8c36 C11 trees (truth = STAR conformation) and the m8c39 toy split replays
(route.py / route_truth.py in the scratchpad).

The comparison: today's children (the EM's particle labels) vs routing whole stage-B states (each state to the child holding most of
its particles).

- CAHRA root, both seeds: stage-B state 1 is 7.3k at 0.68 open and holds 0.90 of all open. The EM moved a third of the node off
  their state's child, at the base-rate open fraction (its axis is orthogonal to open/closed).
  - The open-rich child: EM precision 0.29 / recall 0.75; state routing 0.685 / 0.895.
  - r2: state routing is also better (recall 0.68 -> 0.77).
  - The open-poor s2 r1 / r11: the EM wins on recall (0.80 vs 0.51).
- Toys: the EM wins at every truth-supported split (junk NMI 0.97 vs 0.89, 0.89 vs 0.70, 0.23 vs 0.12).
  - At p20ct2 r22 and p27ct2 r122 both stage-B states map to one child, so state routing loses the junk split entirely.

**Why:** the EM band (<= 20 A, kept below the re-test band) sees low-res differences (junk, composition) but not high-res ones
such as CAHRA's pore helices (< 19 A).

**How to apply:**
- Don't switch to state routing as a rule.
- BUILT (correction 2026-10-25): route_members in state_mergeup%test, behind the env switch SIMPLE_SPLIT_ROUTE=best (uncommitted
  prototype). After K* is fixed: one full-band cross-half pass (shells 2..khi) from the EM groups and from the mapped stage-B states;
  the higher log-likelihood decides.
- Tested 2026-10-24 (m8c42t, stree_splitrep2 jobs 54181011 ctl / 54181012 route):
  - It FAILS at the CAHRA c36 roots: it picks the EM groups by about -31k and open NMI falls 0.059 -> 0.011.
  - It wins at c30 r2: stage-B is chosen and open NMI rises 0.57 -> 0.77.
  - It never hurts the toy junk splits (p20ct2 r22 0.22 -> 0.30, p27ct2 r122 0.21 -> 0.27). The full band is dominated by low-res shells.
- 2026-10-25, the m8c77t CAHRA DESIGN root lost open this way: stage-B labels NMI 0.377 (recall 0.94) -> children NMI 0.023. The
  baseline escaped by chance. Next variant on the same replay battery: compare the starts' likelihood in the RE-TEST band only
  (kt..khi, where the split was validated).

- 2026-10-25 battery (split_replay/m81_ctl, m81_route, m82_ctl, m82_band; routerep4_list.txt). band = pick by the re-test band:
  - CAHRA: stage-B at every root. Open NMI: c36 0.06 -> 0.12, c30 r2 0.57/0.59 -> 0.77, M77d 0.023 -> 0.042 (stage-B labels 0.377).
    The full-band pass after the pick still moves open particles off.
  - Toys: that same pass after an EM pick spreads the junk into the big child (p27ct2 r122 junk P 0.750 -> 0.205, p27cd2 r221 0.800 -> 0.275).
  - OLE cross-seed ARI r1/r2: ctl 0.668/0.772, best 0.648/0.874, band 0.648/0.552 (r2 pick flips between seeds at +16 / -9).
  - band_direct (m8c84t, arm m84_bdir): CAHRA open M77d 0.377, c36 0.46, c30 r2 0.63/0.66. Toys as band (junk spread). OLE ARI
    0.516 / 0.171.
  - choose (m8c85t, m85_choose): raw label sets compared in the re-test band, stage-B only beyond 3 SE.
    - Toys = ctl (junk kept). M75b / c36 / c30 take stage-B (0.46 / 0.46 / 0.63-0.66). M77d keeps EM (+182 < 3 SE 231): 0.023.
    - OLE: stage-B wins on raw labels by 28-35 SE everywhere, so ARI falls to 0.516 / 0.590 (ctl 0.668 / 0.772).
  - LESSON: in-seed fit does not predict cross-seed agreement. OLE's stage-B states fit better per seed but reproduce worse than
    the <= 20 A EM children. All variants stay env-switch prototypes; nothing built into trees.
  - DECISION (coordinator, 2026-10-25): the routing line is CLOSED. The default stays the EM routing (the control), which reproduces
    best on real data (OLE). CAHRA's loss is a toy-specific high-res state (pore helices < 19 A). RNase P's arm swing is visible at
    low resolution. The three switches (best / band / band_direct / choose via SIMPLE_SPLIT_ROUTE) stay off-by-default prototypes.
    Revisit only if a real-data case shows a high-res state lost this way (the scorecard tracks CAHRA's open recall).
- FINAL (2026-10-25, m8c84t / 85t): band_direct (stage-B labels whole when the band picks it) and 'choose' (raw labels, pick only
  above 3 SE) recover CAHRA open (c36 roots 0.06 -> 0.46) but cost OLE cross-seed ARI (r1 / r2 0.668 / 0.772 -> 0.516 / 0.171-0.590).
  On the raw labels stage-B beats EM by 28-35 SE at every OLE split, yet reproduces worse across seeds: in-seed fit does not
  predict cross-seed agreement. DECISION (coordinator): the EM default stays; all switches stay off-by-default prototypes.
Related: [[split-test-holdout-and-drain-traps-2026-10-21]].
- REOPENED 2026-10-28 by a real-data case (RNase P closed): on route-A r111 / base fixture r1 / strongest-pole r1 the EM groups fuse
  15.5k / 7.1k / 8.2k of s01 into the closed-corner stage-B state s02; existing modes can't route (both states map to one EM group).
  route9_on_v5.diff: when K* == #stage-B states the states ARE the groups; `choose` then delivers s02 intact on all three (margins
  3-60 x the 3 SE), nulls K*=1, qtoy2 r1211 c18 purity 0.91 -> 0.98, junk split unchanged, m33 pair kept. OLE cross-seed ARI not
  measurable on 6_/7_ (ARI ~0 in every arm). Doc 10.15; beegfs hetero_tree_20261028/route9.
