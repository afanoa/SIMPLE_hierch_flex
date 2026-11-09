---
name: overlap-toy-soft-vs-hard-2026-11-09
description: "Overlap toy (2026-11-09, coordinator + user): 5 oracle-d' levels of an A/B parent cut x fix/null with a rare A-only child; stock tree vs f2r's textbook soft child refit, 3 seeds; work dir hetero_tree_20261028/overlap_toy, notes handoff/overlap-toy"
metadata:
  node_type: memory
  type: project
  originSessionId: fbba56ec-1a15-49c3-858f-0fa9e980370c
  modified: 2026-11-09T19:22:51.583Z
---

2026-11-09: the coordinator, with the user's approval ("we should be using many more cores"), asked for a two-level fixture family
with tunable parent-cut overlap, to measure WHEN the soft child refit helps.
- **Design:** RNase P conf01 base at box 64 / 3.66 A.
  - Parent: a 27 % domain hinged +-h; d'_eff .74 / .99 / 1.48 / 1.98 / 2.85 at h 3.5 / 4.9 / 8.3 / 13.2 / 29.9 deg.
  - Child: an 8 % domain, a 45.4 deg rare state on 6 % of all (A only), strength d' 1.86, matched to conf18's 1.91 at the pool noise.
  - 60k particles, SNR0 .0012, shared draws (paired).
- **Calibration:** oracle d'_eff = -2 Phi^-1(mean Phi(-D_i / 2)) from noise-free component stacks (simulate_particles is linear,
  `ov_calib.py`).
  - TRAP: d'_eff saturates with swing angle (broad view-dependent D_i), so d' 3 was out of reach and 2.89 is the clean level.
  - TRAP: placing domains by the smallest tip radius picked a groove (rim 18 A); place them at the outermost rim point whose tip fits
    the mask.
  - TRAP: draw each domain's texture with its own generator, or the build's child differs from the calibrated one.
- **Arms:** hard = stock tree (install5, nl_fixtree keys, maxdepth=1). Soft = f2r's textbook loop at every root child.
  Hard control = the same replay with the child's own fflex.

**Why:** RNase P r131 was a poor test of soft refitting (little lost at the parent cut); see [[textbook-flex-weighting-production-level-2026-11-09]].
**How to apply:** read results from `ov_score.py summary` and handoff/overlap-toy. Related: [[rare-state-detection-limit-2026-10-29]],
[[soft-carry-no-gain-on-fixtures-2026-11-09]].
