---
name: igm10-fixture-real-fab-proportions-2026-11-07
description: "IgM fixture with 10 % arms (fixture_igm10, real Fab proportions) built 2026-11-07 by linear reuse of the 5 % noise-free stacks (exact, x2); noise recalibrated (band power x1.12-1.43); 10 % stock roots split both twins (fix a1 NMI .072, null .505 false positive); same-build 5 % pending; STATUS in handoff/igm10"
metadata:
  node_type: memory
  type: project
  originSessionId: fbba56ec-1a15-49c3-858f-0fa9e980370c
  modified: 2026-11-07T08:22:08.435Z
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
- The baseline stock roots run on merge_master2 with nodepar 4: 10 % jobs 55016631/2, plus the same-build 5 % re-runs
  55016840/1 (coordinator-approved). Each tree job scores itself (`score_fx.py`) into handoff/igm10/score_<TAG>.txt. The 10 % roots
  also queue T1: true-arm and stock masks only, since there are no island sites at 10 %.
- Interim at 03:17, 10 % roots on merge_master2: both twins split in two and the guard kept both splits.
  - fix: NMI a1 .072, a5 .015, others <= .007.
  - null: largest share .505, a false positive; it follows the continuous a1/a3/a5 angles (eta2 .03-.05).
  - The build_igm4 5 % roots held one leaf each.
  - Signal or build is open until the same-build 5 % scores land.
- The 10 % trees finished at 03:20 with these numbers. Their whole-molecule fflex latent saw nothing of the arms (R2 <= .006,
  LDA at chance). T1 was queued: 55020271 (fix), 55020268 (null).

**Why:** separates an SNR limit from a method limit on the IgM-regime tests ([[site-mask-rework-direction-2026-11-07]]).
**How to apply:** for a fixture variant that only rescales a component, reuse the noise-free component stacks by linearity
instead of re-simulating. Read results from handoff/igm10/score_*.txt and X10/t1_sitefocus. Related:
[[igm-fixture-xval-stageb-gapw-2026-11-06]], [[igm-root-mask-is-flex-solvent-imprint-2026-11-05]].
