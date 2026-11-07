---
name: effect-floor-pn2-after-guard-2026-11-02
description: "GO (577ee25fc on claude/tree-lever-floor): effect floor after the locality guard merges a child pair with P/N < 2 (children's difference power over the node's half-map noise power, the guard's own statistic); merges only the two within-conformer qtoy2 splits (NMI 0.773 -> 0.812), r131 corner margin 1.3x at both seeds, nulls untouched; Fortran = Python port to 1-10 %"
metadata:
  node_type: memory
  type: project
  originSessionId: b4673f01-1624-41b8-8544-1b32f8b5b56a
  modified: 2026-11-06T06:38:00.114Z
---

2026-11-02, session tree-lever-floor (note hetero_tree_reworks_2026_11.md section 10; state_mergeup effect_floor / pair_maps,
EFF_MIN = 2; replays math_reworks/lever/splits*, arms head vs bc):
- P/N on the same saved df nodes, Fortran vs detlim2/scripts/split_effect.py: qtoy2 root 42-82 vs 44-88, r131 root 7.8 vs 7.7,
  closed corner 3.18 vs 3.05 (s2) and 2.67 vs 2.66 (s1), OLE depth-1 5.8-6.9 vs 5.8-6.7; same verdict at every node.
- bc arm: floor fires only at qtoy2 r3 (combo 28, P/N 0.22) and r5 (combo 101, 0.22): 10 -> 8 leaves, NMI 0.773 -> 0.812,
  combo F1 0.64 -> 0.94 / 0.71 -> 0.95; r131 both seeds untouched (corner 2.6-2.7; the 0.93 closed-arm r2 split at 2.08 kept);
  nulls merged by the guard before the floor; detlim2 8 % root 9.8-9.9; OLE untouched (root pairs 31-50, depth-1 5.3-7.9; the
  8 leaves are the head arm's to the particle except 95 at r3).
**Why:** with the draw-free split test as default (c7ac3b0bd) the one marginal within-conformer split per tree needed a
map-space check; localization has no power left among delivered splits, effect size does (rare-state note 3.2).
**How to apply:** keep EFF_MIN at 2 (rare-state detections at the limit sit at P/N 1-2; the r131 corner at 2.6-3.2);
read a floor merge in the split log ("effect floor: groups g / h ... -> merged"). Related: [[split-test-single-draw-variance-2026-11-02]],
[[rare-state-detection-limit-2026-10-29]], [[took-off-site-class-is-no-pole-2026-11-02]].
2026-11-06: re-measured with a guarded FSC match; the margins are real, not g noise. The floor alone would keep every null root
(P/N 6-41), so only the guard stops the nulls ([[pair-maps-fsc-match-unbounded-trap-2026-11-06]]).
