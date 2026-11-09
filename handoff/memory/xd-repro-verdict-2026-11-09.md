---
name: xd-repro-verdict-2026-11-09
description: "VERDICT (2026-11-09): XD of a node's flex latent reproduces across seeds, but so do hard partitions of it (flex deconv ARI .95-1, XD MAP .98-.99); flattening isn't the gap (MAP sharpens closed); the split's seed dependence is downstream; Gaussian-null K rule not null-safe"
metadata:
  node_type: memory
  type: project
  originSessionId: 1218bb2b-a21a-4d25-9e48-c40ddd89ceff
  modified: 2026-11-09T05:29:15.092Z
---

**The test.** xd-repro: STATUS is `handoff/xd-repro/STATUS.md`; the work dir is `math_reworks/xd_repro`.
- The user's thesis (2026-11-08): per-particle embeddings are mostly noise, and flattening responsibilities into hard labels is
  the gap between flex and refinement. So population-level objects should reproduce where hard partitions don't.
- Method: deterministic XD on each node's own flex cache, using flex's deconvolution likelihood. K by the cross-half held-out
  gain, which must beat max(0, max of 19 Gaussian-null gains).
- Run on the 16 r1 node-levers replays, plus the fixture and null roots.

**Results.**
- **r1, same latent** (seeds, class prior, builds): XD components identical (map r 1.000, pi to .001, resp r 1.000). XD MAP labels
  (ARI .98-.99) and flex's own deconv labels (ARI .95-1.00) reproduce too. The tree's split children don't: ARI .20-.77, hit or
  miss.
- **r1, across flex's noise model off vs on**: maps r .956, |dpi| .10, resp r .75, against ARI .24-.32. This is the only place
  the thesis's prediction holds.
- **Closed at r1** (base .0194):
  - XD at K=3: best component x1.29-1.30 soft, x1.52-1.56 as MAP labels;
  - flex deconv (K 4-5): one 14-17 % label at x1.9-2.0 (recall .28-.33), on every seed and both latents;
  - the split: hits x1.84-2.00 at recall .67-.70 (arm a seed 6: x3.59 / .83); misses x1.08-1.24.
- **Fixtures:** qtoy2 K=11 NMI .502 (root split .665); qtoyflex K=2 .110 (root split .169); IgM 10 % fix K=5 carries no arm
  information (NMI <= .001).
- **Nulls:** qtoynull K 2, olenull 2, micnull 4, IgM null twin 7. The Gaussian-null ladder detects a non-Gaussian latent shape,
  not states.

**Why:** the coordinator asked whether population-level soft objects should replace hard partitions in the tree.

**How to apply:**
- Don't propose soft-vs-hard, or XD of the latent, as the reproducibility lever. The seed dependence is in the site/pole choice
  and the stage-B seeding ([[node-levers-replay-verdict-2026-11-08]]).
- flex's own deconv labels are the reproducible latent seed candidate at r1 (x2.0); untested as a seed for the round.
- Any K read off the latent needs a test calibrated on the fixtures' nulls.
- Comparing maps across runs needs [[flex-basis-and-z-signs-arbitrary-across-runs]].
