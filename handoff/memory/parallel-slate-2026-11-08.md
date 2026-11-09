---
name: parallel-slate-2026-11-08
description: "USER 2026-11-08 night: 'we have a lot of compute availble so please test many things in parallel' -> six sessions on six tests (joint seeding refined; XD reproducibility; r13 subtree reliability; 692k root reliability; IgM joint seeding; OLE panel); guide 7.10 thread 15"
metadata:
  type: project
---

Launched ~00:30 2026-11-09, with about 1,000 cores free in our caps.
- flex-refine-flex-weights: refined joint seeding (r1 x3, qtoy2, qtoyflex, nulls), responsibilities, r1 seeds 3-6, r131 honest arms.
- tree-latent: deterministic XD on full node latents, reproducibility against hard partitions. Notes $H/handoff/xd-repro/.
- extremes-repro: r13 subtree to leaves at 3 seeds, the closed-leaf reproducibility baseline. Notes $H/handoff/subtree-repro/.
- rnasep-misfit: 692k root replay at 3 seeds. Notes $H/handoff/root-repro/.
- site-chain: joint seeding at the IgM fix roots and nulls. Notes $H/handoff/site-chain/joint-igm/.
- sparse-ppca2: OLE maxdepth-1 trees at 4 seeds, node trace, joint seeding. Notes $H/handoff/ole-panel/.
**Why:** the user prefers broad parallel tests over serial ones when cores are idle
([[dont-gate-tests-on-price-cores-idle]], [[use-the-full-cpu-cap-for-tests]]).
**How to apply:** fold each verdict into guide thread 15 as it lands. Keep the shared rules (never rm, fresh folders, production
read-only, scancel by id). Watch the caps, csbdevel 642 / norm 600.
Related: [[near-tie-is-duplicate-factor-masks-2026-11-08]], [[methods-must-generalize-across-datasets-2026-11-08]],
[[replicate-or-redesign-positives-only-2026-11-08]].
- ~02:30 results:
  - r13 subtree: r13's split reproduces at 3 / 3 seeds (r131 side x1.73-1.74; partition ARI .94 between seeds, .90 with
    production). The instability is at r1, not r13.
  - Gaussian-null K floor: NOT null-safe (nulls get K 4-5). K must be decided by the map-level split test.
  - Site chain: IgM 5 % fix gains small over stock (.022 vs .001); nulls stay one leaf.
  - User "Yes, test it": a final argmax settling pass at the end of stage B, on a test build, at r1 over 6 seeds.
- r131 split (3 seeds): NOT reproducible at the particle level. Every seed makes a near-half split with a closed-leaning half
  enriched x1.16-1.57, but those halves overlap at Jaccard .22-.39 (near chance).
  - So the closed lineage reproduces r13 -> r131 and is arbitrary below. The corner leaves are not reproducible sets.
