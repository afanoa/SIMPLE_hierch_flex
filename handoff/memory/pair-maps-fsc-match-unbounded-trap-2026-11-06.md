---
name: pair-maps-fsc-match-unbounded-trap-2026-11-06
description: "TRAP: state_mergeup pair_maps' per-shell FSC match g = cross/fb is unbounded; past the pair's resolution (fb ~ 0) g hits 39-74 and the CLIPPED floor / guard sums gain a random positive bias. Re-measured 2026-11-06: the guard is inert at all 69 gate nodes (no verdict changes; corner 3.08 / 2.63 real)"
metadata:
  node_type: memory
  type: project
  originSessionId: c79bc433-5d25-4389-9589-67d90cec0604
  modified: 2026-11-06T06:37:48.050Z
---

2026-11-06, IgM fixture session (report tree_igm_fixture_report.md, section "TRAP found on the way"; diagnosis scripts
fixture_igm/scripts/diag_g.py, proto_pn.py GMODE=sig).

- Code: src_igm `simple_state_mergeup.f90:5208` `pair_maps`: g(sh) = 0.5 (c14 sqrt(pw1 pw4) + c23 sqrt(pw2 pw3)) / (fb sqrt(pw3 pw4)),
  no bound (g = 0 only when den <= 0). The band is 10-40 A. ex = max(smoothed(res^2) - smoothed(nse^2), 0); p = sum(ex) feeds the effect
  floor (EFF_MIN 2), and top_share(ex) feeds the locality guard.
- On four 50/50 site splits (box 200, 3.27 A/px, 100k particles) |fb| < .02 at 11-13 A, and g = 39.3 / 42.3 / 56.8 / -73.9 at one
  shell each (median .95-.98).
  - Clipped vs unclipped sphere P/N: .36 / .12 (the one run whose max g was 3.1), 6.0 / -.08, 8.9 / -.22, 32.2 / .66.
  - In-site, null: 18.1 / .30 and 23.7 / -3.96.
- Guard prototype (not ported): match only where fb > 3 / sqrt(N_shell / 2), else g = 1.

**Why:** the clipped sum turns zero-mean noise amplified by g^2 into a positive effect, so the floor or guard can pass a split on noise.
The r131 closed corner sits at P/N 2.6-3.2 ([[effect-floor-pn2-after-guard-2026-11-02]]), so the real gate margins are unknown until
they are re-measured with a guarded g.

**How to apply:** read pair_maps P/N values only next to the per-shell g (or the unclipped sum). Related: [[split-test-single-draw-variance-2026-11-02]].

RE-MEASURED 2026-11-06 (side session; note side-session-notes/local_5bd42cba.../pair_maps_guard_gates_2026_11_06.md; batch2
hetero_tree_20261028/pn_guard_20261106/, pn_guard.py mirrors the Fortran to <= 0.24 % on 166 logged values):
- The guard is INERT at every gate node: 69 nodes in the bc / df / base / guardfix arms, 231 pairs.
  - The stage-B group maps resolve the whole band: fb >= 4.97 x the bar at every shell, |g| <= 1.15, fb at 10 A 0.46-0.94.
  - So no verdict changes, floor or guard. The margins are real: corner 3.08 (unclipped 2.99) at s2 and 2.63 (2.52) at s1; qtoy2
    r3 / r5 0.22 merged; the tightest kept split is r131 s2 r2 at 2.08 (2.00).
- IgM fixture trees (stage-B box 164): the guard touches one pair, the real root's guard ratio (.480 -> .475, merged either way).
  The 39-74 g came from the xsb prototype's box-200 maps.
- Null roots: the floor alone would KEEP them (P/N 6.1 / 9.9 / 40.8, clipped = unclipped); only the locality guard (.20-.24)
  merges them.
- Porting the guard needs no new gate run.
