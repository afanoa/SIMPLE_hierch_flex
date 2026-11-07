---
name: micelle-fixture-pose-scatter-verdict-2026-10-29
description: "Synthetic PfCRT+micelle fixture: the tree's own stA0 re-alignment flips 29-32 % of root particles > 30 deg (real PfCRT 41 %) and every split is aligned-vs-misaligned; POSEGATE=0.165 -> fixposes makes nulls pass; 8/16 deg protein states lost by the <= 19.5 A EM band"
metadata:
  node_type: memory
  type: project
  originSessionId: 9643ce3d-fc7a-4dc1-860a-59db10867e7c
  modified: 2026-10-30T00:34:21.204Z
---

2026-10-29, micelle session (tree_micelle_fixture_2026-10-29.md in the shared notes; W = hetero_tree_20261028/fixtures/pfcrt_micelle;
datasets micnull / mic3, 56,808 ptcls, box 120 at 2.6 A, noise band-calibrated to real PfCRT; ~196 core-h).
- Pose scatter: stA0 flips 29-32 % of the root's rows > 30 deg (real PfCRT tree_umask* root 0.41). Splits = aligned vs misaligned
  (z 15-73); masks localise them, so the guard keeps them and nulls fail; with no mask the guard merges them (null passes by accident).
  Micelle AUCs 0.48-0.52: the micelle never drives those splits.
- Flip survey (flips/flips.tsv, 321 nodes): legitimate stA0 <= 0.0245 (OLE, IgM, qtoy2, nulls); RNase P 692k root stage-B 0.094 (bound).
- SIMPLE_TREE_POSEGATE=0.165 (root-only, default off): fires on both fixture roots and switches to fixposes; nulls pass.
- EM band: the cross-half EM (<= 19.5 A) relabels 28-53 % even from truth; the 8/16 deg rocking (2.4/4.8 A) is never found.
  0 vs 16 deg passes the re-test only with SIMPLE_TREE_SUPPORT=order and a truth start (z 2.26 -> 4.70).
- The order support keeps the inner belt (95 % at 20-30 A): a seeded micelle split passes more strongly (z 3.05 -> 11.97), guard 0.575.
- Code uncommitted (worktree tree-micelle-synth-24763a); builds recipe17d_micsup / _fp / _fpg, gate 15/15. Not tested ON on RNase P/OLE.

**Why:** the user asked (2026-10-29) for a synthetic PfCRT with a micelle to develop on.
**How to apply:** treat real PfCRT tree leaves as likely pose artefacts; any combined method needs a pose-stability check and a per-node
data-derived band that keeps the 14.6-10.2 A re-test shells. Related: [[tree-user-mask-opt-in-not-default]], [[tree-data-derived-lp-verdict-2026-10-29]].
