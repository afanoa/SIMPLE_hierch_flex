---
name: tree-nu-pooled-maps-no-go-2026-11-02
description: "non-uniform tree-pooled leaf maps (nu_filter level selection over the leaf tree) lose to the tree-Wiener maps and dilute the closed corner; NO-GO, prototype on batch2 nu_pooled/"
metadata:
  node_type: memory
  type: project
  originSessionId: bc0680a2-d45d-49ac-b343-ae9e57341c60
  modified: 2026-11-02T05:37:39.840Z
---

Verdict 2026-11-02: a non-uniform final-map scheme (per voxel choose the leaf-tree level to pool at, nu_filter unary + Potts,
Wiener weight from the pool's FSC) is a NO-GO against `simple_tree_wiener`. On the 692k r17d plain leaves: held-out FSC
(even-built vs own odd, 40-10 A) NU .946/.774 vs tree-Wiener emp2 .958/.808 (plain .943/.759); the closed-corner leaf
r1311_s2 drops from closed-arm 1.00 (plain) / 1.03 (emp2) to 0.92 in every variant (raw band, 10 A band, squared loss, full-map
correction) because 99-100 % of its closed-arm voxels pool with r1311_s1 (closed 0.77); the own level wins 0 % of the envelope
for 64 of 67 leaves (noise-driven: only the 18-33k leaves keep own voxels).

**Why:** a local MSE-optimal choice between nested pools pays the sibling bias wherever it is below the 10 A-band half-map noise
in the 15 A window, which is exactly at the rare state's distinguishing feature; the tree-Wiener weights the leaf's reproducible
signal per shell globally and keeps it.

**How to apply:** do not revisit level-selection pooling; if a non-uniform final map is wanted, prototype a local tree-Wiener
(the covariance fit per window). Prototype and outputs: batch2 `/mnt/beegfs/elmlund/afan/hetero_tree_20261028/nu_pooled/`
(runs/{raw,lp10,lp10sqnc,lp10sq}, labels/unit_NN_level.mrc), repo branch `claude/tree-nu-pooled`
(`doc/implementation_notes/tree_nu_pooled_maps_2026_11.md`, `scripts/hetero_tree_nu_pooled_proto.py`). See
[[rnasep-closed-state-must-be-retained]], [[hetero-tree-branch-published-on-fork]].
