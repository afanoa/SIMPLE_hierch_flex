---
name: igm-fixture-xval-stageb-gapw-2026-11-06
description: "IgM-like fixture verdict: (b) islands + cross-validated per-site stage B (swap + per-step half alignment) + split only on GAPw DISCRETE = first config with fixture gain (a1 NMI .146) AND null held; null continuum cuts are REAL (eta2 .06-.07), rigid sites split by shift error (R2 .32); not ported"
metadata:
  node_type: memory
  type: project
  originSessionId: c79bc433-5d25-4389-9589-67d90cec0604
  modified: 2026-11-06T05:43:01.911Z
---

2026-11-06, IgM fixture session (report: side-session-notes/local_5bd42cba.../tree_igm_fixture_report.md, sections BATTERY,
SCORECARD UPDATE, RECOMMENDATION; work dir batch fixture_igm/xsb/, scripts xstageb.py v2, xsb.sbatch, xsb_battery_*).

- Cross-validated stage B needs the halves' state numbers aligned every step: separate per-half seeds left every first run
  anti-aligned. That flipped the cross-half r, swapped P/N's signal and noise terms, and set the swap design into a 2-cycle that
  relabelled 50-70 % of particles per iteration.
- With alignment, swap is the best site classifier tried: a1 NMI .146 (seed 2 .138), a2 .047, a5 .039. The ±20 deg arms a3 / a4
  are not found.
- The null's continuous arms are cut for real (a1 / a5 eta2 .06-.07): every reproducibility or P/N gate passes them just as it
  passes the fixture's real arms.
- A rigid-core patch splits at z 26 by the particles' SHIFT ERROR (R2 .32), in both swap and indep, so the swap test is not inflated.
  - The 15 A site/rest P/N (< 1) and GAPw (not discrete) stop it.
  - The tree-band relative (2.5) and the cross-half z do not.
- GAPw (sibling, [[discrete-vs-continuous-site-gap-2026-11-05]]) gives fix a1 DISCRETE and no false DISCRETE on 9 null / control
  runs. Most arms are "no call" (lobes too close).

**Why:** the user's gate is "null twin unsplit, fixture's discrete arms the target". This is the first configuration that meets
both, at one arm of five.
**How to apply:**
- Gates replayed 2026-11-06, both PASS, layer inert on RNase P: [[igm-site-layer-gates-pass-add-only-2026-11-06]].
- Real tb100k has no arm islands yet.
- Related: [[pair-maps-fsc-match-unbounded-trap-2026-11-06]], [[igm-root-mask-is-flex-solvent-imprint-2026-11-05]],
  [[tree-good-pose-null-fails-halves-and-seeding-2026-10-30]].
