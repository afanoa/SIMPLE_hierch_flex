---
name: continuum-inherit-verdict-2026-10-29
description: "Continuum fixture (30 A swing, dwell + uniform): recipe17d children lose the axis after the root cut; SIMPLE_TREE_INHERIT=yes reaches the closed end (seed-reproducible) but not the far end, interior/dwell not reproducible; leaf counts != occupancy (pose bias pushes maps outward); guard blind spot found"
metadata:
  node_type: memory
  type: project
  originSessionId: 9643ce3d-fc7a-4dc1-860a-59db10867e7c
  modified: 2026-10-30T00:39:03.370Z
---

2026-10-29, continuum session (continuum_2026_10_29.md + inherit_on_r17d.diff md5 01764c0f in the shared notes; runs
hetero_tree_20261028/continuum; fixtures ctpk (dwell at t 0.65) / ctun (uniform), qtoynull recipe, 60k; ~340 core-h).
- Baseline recipe17d: root cuts at the mass median (z 25-52, guard kept); children's own discovery can't see the residual spread
  (per-particle position sd ~0.3 of the span), fflex fallback seeds an unrelated axis -> 2-3 broad leaves, ends not reached.
- Inherit (seed a pole-less node by the parent's stage-B contrast): all depth-1 children split along t; closed end reached in all
  4 runs (maps 0.00-0.08), far end 0.81-0.88 (not within 0.1). ctun seeds 1/2 identical (ARI .83); ctpk interior differs (seed-2
  root took a route-A extra seed; dwell map 0.68 vs 0.51).
- Occupancy: leaf maps sit 0.1-0.3 further out than their members (pose error grows with distance: 3.4 vs 7.6 deg); counts, hard
  routing and soft mixture weights all fail to recover the dwell (0.28-0.33 vs true 0.64). Needs conformation-insensitive poses.
- qtoy2: only the 2 pole-less nodes change; c1 F1 .54 -> .77, c2 .65 -> .56, holo .94 -> .85 (holo over-split passed via the guard
  blind spot). The fallback seeds 71 of 84 nodes in the 692k integration tree, so INHERIT would change most of RNase P.
- GUARD BLIND SPOT (in recipe17d = pushed 5709bc083): locality_guard skips silently when choose routes to EM groups and every
  stage-B state's majority lies in one EM group (no spanning pair) -> K* kept with no L/L_ref logged.

**Why:** the user asked (2026-10-29) for leaves that reach the extremes and sample the continuum by occupancy.
**How to apply:** INHERIT stays default off until the RNase P closed gates (r1211, closed_r111) and OLE pass; never read continuum leaf
counts as lifetimes. Related: [[discrete-splits-capture-the-span-states-that-matter]], [[null-split-diffuse-real-split-localized-2026-10-28]].
