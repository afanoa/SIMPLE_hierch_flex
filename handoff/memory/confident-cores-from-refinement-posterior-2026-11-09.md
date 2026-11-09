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
- Selection maps (own-data, imprint caveat; ~1.3 core-h):
  - r13: the >= .95 core (20,155) reads closed / open 1.023 / .604 at recall .262. Production's best leaf r1311_s2 (12,010) read
    1.004 / .577 at recall .094, so the core matches corner-level occupancy at 2.8x the recall. Halves agree to +- .004.
  - r1: contrast c rises with the threshold: hard child +.107, >= .9 +.285, >= .95 +.335, >= .99 +.423 (resolution 5.7 -> 7.8 A).
  - The confident cores are closed-enriched MAPS, consistent with the reference-row curve. Honest X2 read and cross-seed
    reproducibility are pending.
- CORE REPRODUCIBILITY (3 r13-subtree seeds plus production):
  - r13: the >= .9 cores replicate at Jaccard .88-.89 (the hard halves at .91). 16,858 particles are shared by all four runs at
    12.5 % closed (x3.03), about 1.8x purer than production's corner leaf; >= .99 gives 5,180 shared at 16.6 %.
  - r131: cores .07-.21 (the halves .25-.42); only 158 shared. The split is not reproducible, so neither are its cores.
  - VERDICT: confidence adds purity, not reproducibility.
  - Emerging design:
    - joint seeding in the latent axes;
    - accept a split only if it replicates across seeds or half-sets, which stops the descent at r131;
    - confident cores with responsibilities as the deliverable;
    - a settling pass.
  - The settling pass reaches a fixed point after 2 cold passes (3-4 % of directions moving).
- r1 cores (6 settled arm-b runs):
  - Within the hit basin, >= .9 cores are 22.6k at .069-.076 (x3.5-3.9), Jaccard .77-.87; across basins .18.
  - The intersection over all 6 runs is 12,308 at .0825 (x4.25), recall ~.12. A seed-consensus core is reproducible by
    construction and purer than production's r131.
  - Cores remain a purity lever conditional on a replicating split.
