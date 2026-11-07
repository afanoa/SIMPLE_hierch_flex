---
name: pfcrt-micelle-fixture-tree-pose-scatter-2026-10-29
description: "PfCRT micelle fixture (hetero_tree_20261028/fixtures/pfcrt_micelle): the tree's own stA0 re-alignment flips ~30 % of root poses (real PfCRT 41 %) -> aligned-vs-misaligned splits; EM <= 19.5 A loses fine protein states; SIMPLE_TREE_POSEGATE / SIMPLE_TREE_SUPPORT prototypes, default off"
metadata:
  node_type: memory
  type: project
  originSessionId: c82ee840-ba2c-4b73-a246-e79b0fb2c726
  modified: 2026-10-30T00:34:36.642Z
---

2026-10-29, session tree-micelle-synth (coordinator local_97df86ed). Hand-over: shared notes tree_micelle_fixture_2026-10-29.md.
- **Fixture:** micnull and mic3, 56,808 particles, box 120 / 2.6 A, sharing per-particle micelles.
  - Protein: 6UKJ atoms; mic3 adds a half-transporter rocking 0 / 8 / 16 deg at 45 / 35 / 20 %.
  - Micelle: the real belt warped in 256 bins.
  - Real poses and CTF, band-calibrated noise.
  - Tree inputs `<ds>/src_in(_fp).simple`.
  - Scorer `micelle_score.py`, launcher `mic_gate.sbatch`.
- **Finding 1, pose scatter.** On micelle-dominated small proteins the tree's local re-alignment turns 29-32 % of the root's particles > 30 deg at stA0. Real PfCRT tree_umask runs: 0.41.
  - Splits are then aligned vs misaligned (z 15-73); masks make the guard keep them.
  - Legitimate trees: stA0 <= 0.025 (OLE, IgM, qtoy2, nulls); RNase P 692k root <= 0.094 at stage B.
- **Finding 2, the EM band.** The cross-half EM (<= 19.5 A) relabels 28-53 % from a truth start, so the fine protein states are lost (8 deg never found).
- **Finding 3, the inner belt.** The order support keeps the inner belt, so seeded micelle splits pass more strongly; the guard keeps micelle splits (0.575).
- **Prototypes** (builds recipe17d_micsup / _fp / _fpg, each gate 15/15; uncommitted worktree tree-micelle-synth-24763a):
  - SIMPLE_TREE_POSEGATE=X (X 0.165 from the survey gap): fires at the root's stA0 and switches the tree to fixposes.
  - SIMPLE_TREE_SUPPORT=order: input-free structured-support mask used as tree_mskfile.
  - Both unset = off.

**Why:** the micelle is not what breaks the tree here; pose stability is. Any band change must keep the 14.6-10.2 A re-test shells, because the RNase P closed differences live at 20-40 A.
**How to apply:** check a node's stA0 flip fraction before reading splits on small membrane proteins; reuse this fixture for those questions. Related: [[pfcrt-core-alignment-2026-10-29]], [[tree-user-mask-opt-in-not-default]], [[synthetic-toy-traps-noise-poses-2026-10-29]].
