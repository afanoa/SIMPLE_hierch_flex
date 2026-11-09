---
name: confident-cores-from-refinement-posterior-2026-11-09
description: "FIRST RNase P PURITY GAIN (reference labels): selecting particles by refinement responsibility beats the hard children - r13 resp>=.9 keeps 25k at .102 closed (x2.47, recall .295) vs production r1312 .083 / .209; r1 resp>=.9 30.6k at .066 (x3.38); validation (cross-seed core reproducibility; honest X2 read) running"
metadata:
  type: project
---

2026-11-09 ~03:30, f2r session, from production's own stage-B evidence (prob_state_evidence.txt; zero compute). It follows
the user's idea "refine the pdf and get responsibilities out", used to rank rather than weight (soft maps dilute contrast 11 %).
- r13 (node share .041):
  - hard child r131: .069, recall .353;
  - resp >= .9: 25,344 at .102 (x2.47), recall .295;
  - >= .95: .113 / .262;
  - >= .99: 11,518 at .139 / .184.
  - Production's next level, r1312, has .083 / .209.
- r1 (share .0194): hard r13 .041 / .380. resp >= .9: 30,567 at .066 (x3.38) / .230. >= .99: 9,285 at .090 / .095.
- The uncertain particles of the hard child are mostly open (.013-.045 closed).

**Why:** the first lever all day that raises closed purity on an external measure, beating production's own lineage.
**How to apply:** treat confident cores (with responsibilities) as the candidate deliverable for rare states. It ships only if
the cores reproduce across seeds (the hard r131 halves overlap at only .22-.39) and the honest X2 read clears same-size random
bars. Never weight maps by responsibilities.
Related: [[near-tie-is-duplicate-factor-masks-2026-11-08]], [[replicate-or-redesign-positives-only-2026-11-08]],
[[extremes-fission-verdict-2026-11-07]].
