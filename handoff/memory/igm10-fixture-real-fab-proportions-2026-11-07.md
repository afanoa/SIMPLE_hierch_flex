---
name: igm10-fixture-real-fab-proportions-2026-11-07
description: "IgM fixture with 10 % arms (fixture_igm10, real Fab proportions) built 2026-11-07 by linear reuse of the 5 % noise-free stacks (exact, x2); noise recalibrated (band power x1.12-1.43); stock roots on merge_master2 + T1 queued by the jobs; STATUS in handoff/igm10"
metadata:
  node_type: memory
  type: project
  originSessionId: fbba56ec-1a15-49c3-858f-0fa9e980370c
  modified: 2026-11-07T06:34:12.358Z
---

2026-11-07, side session of the site/mask coordinator: built `igm_cryolo/tree/fixture_igm10` (X10), `fixture_igm` (5 %) with
`ARM_KDA` 100 (each arm 10 % of the core's integrated density). The handoff is `hetero_tree_20261028/handoff/igm10/`: STATUS.md
plus `progress.txt`, which the jobs append to.
- Build in ~25 min instead of ~3 h. `simulate_particles` is linear and adds no noise at snr >= 5, so the 5 % noise-free component
  stacks are symlinked and merged with `ARM_GAIN=2`. Fresh 10 % simulations equal 2x the stacks bit for bit; the core, arm masks
  and truth are identical.
- Two things change, not one. The generator's calibration keeps the real band correlation, so the noise band power rises
  x1.43 / 1.28 / 1.35 / 1.12 (120-30 ... 12-7 A). The arms gain ~3x in power SNR and the core loses 1.1-1.4x. A fixed-noise
  variant would need fixture_igm's `measure_cn.json` and a re-merge.
- The baseline stock roots run on merge_master2 with nodepar 4. The 5 % stock roots ran on build_igm4, so a same-build 5-vs-10
  comparison needs the 5 % roots re-run (command in STATUS). Each tree job scores itself (`score_fx.py`) and queues T1 on its
  root: true-arm and stock masks only, since there are no island sites at 10 %.

**Why:** separates an SNR limit from a method limit on the IgM-regime tests ([[site-mask-rework-direction-2026-11-07]]).
**How to apply:** for a fixture variant that only rescales a component, reuse the noise-free component stacks by linearity
instead of re-simulating. Read results from handoff/igm10/score_*.txt and X10/t1_sitefocus. Related:
[[igm-fixture-xval-stageb-gapw-2026-11-06]], [[igm-root-mask-is-flex-solvent-imprint-2026-11-05]].
