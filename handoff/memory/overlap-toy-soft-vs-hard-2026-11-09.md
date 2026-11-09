---
name: overlap-toy-soft-vs-hard-2026-11-09
description: "Overlap toy (2026-11-09, coordinator + user): 5 oracle-d' levels of an A/B parent cut x fix/null with a rare A-only child; stock tree vs f2r's textbook soft child refit, 3 seeds; work dir hetero_tree_20261028/overlap_toy, notes handoff/overlap-toy"
metadata:
  node_type: memory
  type: project
  originSessionId: fbba56ec-1a15-49c3-858f-0fa9e980370c
  modified: 2026-11-09T19:47:01.799Z
---

2026-11-09: the coordinator, with the user's approval ("we should be using many more cores"), asked for a two-level fixture family
with tunable parent-cut overlap, to measure WHEN the soft child refit helps.
- **Design (corrected 14:45):** RNase P conf01 base at box 64 / 3.66 A.
  - Parent: a 27 % domain hinged +-h; d'_eff .75 / .99 / 1.48 / 1.98 / 2.98 at h 1.37 / 1.84 / 2.87 / 4.07 / 7.43 deg.
  - Child: an 8 % domain, a 40.6 deg rare state on 6 % of all (A only). Its per-particle RMS information is 4.93, matched to conf18's
    5.09 at the pool noise.
  - 60k particles, SNR0 .0075 at box 64, shared draws (paired).
- **TRAP (cost one build):** SNR0 is per pixel. Binning by f averages f^2 pixels, so the box-64 equivalent of .0012 at box 160 is
  6.25 x .0012 = .0075. At .0012, refine3D lost the poses (median 98 deg). A strong reference still gave a 10 A FSC (noise
  aligned to it), so check the pose error, not the FSC. The refine job now gates trees on the median error < 15 deg.
- **Calibration:** oracle d'_eff = -2 Phi^-1(mean Phi(-D_i / 2)) from noise-free component stacks (simulate_particles is linear,
  `ov_calib.py`).
  - d'_eff saturates with swing angle (broad view-dependent D_i). Match detectability by sqrt(mean D^2), not by d'_eff: a local
    change hidden in some views has a low d'_eff at the same total information.
  - TRAP: placing domains by the smallest tip radius picked a groove (rim 18 A); place them at the outermost rim point whose tip fits
    the mask.
  - TRAP: draw each domain's texture with its own generator, or the build's child differs from the calibrated one.
- **Arms:** hard = stock tree (install5, nl_fixtree keys, maxdepth=1). Soft = f2r's textbook loop at every root child.
  Hard control = the same replay with the child's own fflex.

**Why:** RNase P r131 was a poor test of soft refitting (little lost at the parent cut); see [[textbook-flex-weighting-production-level-2026-11-09]].
**How to apply:** read results from `ov_score.py summary` and handoff/overlap-toy. Related: [[rare-state-detection-limit-2026-10-29]],
[[soft-carry-no-gain-on-fixtures-2026-11-09]].
