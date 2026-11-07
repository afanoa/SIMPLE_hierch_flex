---
name: kway-root-binary-below-2026-10-28
description: "user 2026-10-28: no science compromise - pursue route A for closed; proposed design = initial K-way split at the root (route A extras, K data-derived) then binary below; disk secondary but route A result with less disk preferred -> recipe17"
metadata:
  type: feedback
---

2026-10-28 ~21:30: "we dont compromise science. if closed state needs route A we pursue it. What do you think of having an
initial K-split and then binary from there. the disk use is not too much of a concern. if you can get route As result without that
much disk that would be nice."

**Why:** route A's multi-way splits mostly fire at the root (OLE root K=5 both seeds, landscape W1 .022 vs null .100; 692k root
K*=4 concentrating closed into r1; qtoy2 root extras = closed + holo); binary below removes the XD-K near-tie (OLE divergence
cause) and, with no site layer below the root, keeps node projects row-local (base disk).
**How to apply:** recipe17 = recipe16 + route A only at the root + seeds capped at 2 below. First check (labels/text only) where
route A added extras in the closed lineage (memory: route A r11 K*=4 but its closed-corner state was a DRIVING-site state).
Risk: multi-arm nodes below the root (IgM r12 gained 4 arm motions at one node from extras) need more binary levels.
Related: [[site-pole-closed-state-replay-2026-10-27]], [[rnasep-closed-state-must-be-retained]], [[science-before-speed-2026-10-28]].
- 2026-10-28/29 r11 fixture chain (route A r11, 136,643 rows; recipe16 extras vs recipe17 binary, maxdepth 0 per level): extras
  path 33,646 (0.89) -> 17,153 (0.93/0.58) -> 9,241 at 1.07/0.57 = CLOSED CORNER (1,236 refs, share 0.134; r111 ref 25,528 /
  2,800 / 0.110 at 1.00/0.59). Binary path 57,788 (0.81) -> 27,033 (0.84) -> pending. That corner group = leaf of a depth-4 node
  in the full tree, so it is delivered at maxdepth 4 with d>=2. IgM: route A extras at depth 1 resolve 3-4 arm motions per node and
  recipe14 delivers 5 arm extremes mid-tree vs recipe12b 0 at 26 leaves -> kway_maxdepth=2 (recipe17b) is the working setting.
- VERDICT: binary path two splits below r11 = 13,711 at 0.95/0.62 (misses corner on open arm; 1,710 refs, F1 .152) vs extras
  path 9,241 at 1.07/0.57 (corner; 1,236 refs, F1 .138). Both leaves of a maxdepth-4 tree. kway_maxdepth=2 ADOPTED (recipe17b+).
- recipe17b GATES COMPLETE (2026-10-29): nulls pass; qtoy2 NMI 0.706 (16: 0.685), 2.6x less disk; OLE pair reproduces recipe14's
  landscape (17b seed W1 0.044 = 14's; r14_s2~r17b_s2 W1 0.003); IgM 17b NMI 0.447 (14: 0.441, 12b 0.334), extremes 2/20 (14: 3),
  truth span covered 0.54/0.70 (14 0.51/0.66, 12b 0.06/0.13), 32 vs 145 GB. Integration 692k run launched on recipe17d.
