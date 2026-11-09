---
name: near-tie-is-duplicate-factor-masks-2026-11-08
description: "FINDING 2026-11-08: the r1 seed-site near-tie (9.48 vs 9.47) is two factor masks of ONE discovered region (2f1 / 2f2, Jaccard .67, 6 A apart) scoring the same pole; they differ only in pole set (one z1 pole vs z1 + z3 poles), and the closed split needs the z3 pole; OLE's October near-tie 3 is the same duplicate-factor mechanism -> merge duplicate factor masks, keep the union of poles"
metadata:
  type: project
---

Node-lever cheap test, part 1 (f2r session, from existing r1 replay outputs).
- In every near-tie replay, the two sites within 1.0 of the top margin are factor masks 2f1 / 2f2 of the same region: mask
  Jaccard .67-.68, centroids 6 A apart, poles on the same latent tails.
- Pole sets differ: 2f1 keeps one pole (a z1 tail); 2f2 keeps two (z1 tail + z3 tail).
  - All 3-way, closed-concentrating hits came from a two-pole site.
  - The seed-1 miss picked one-pole 2f1 at 9.482 over two-pole 2f2 at 9.466.
- Production's winning site 2 (2f2, margin 16.0) used a merged pole z1hi + z3lo (172k). Its "rest" (80k) became the
  closed-enriched r13; no replay reproduced that merge.
- Same mechanism as OLE's near-tie 3 ([[ole-repro-seeding-near-ties-2026-10-28]]: 1f1 / 1f2 of one 63 A region, equal margins).

**Why:** the run-to-run instability at depth 1 is built into site discovery. Near-duplicate factor masks of one region score
the same pole almost identically, and the pick between them decides whether the closed-carrying z3 tail seeds the split.
**How to apply:** treat near-duplicate factor masks of one region as one site, with the union of their poles, before ranking.
That is the candidate deterministic fix, and it adds states rather than choosing; test it across the panel
([[methods-must-generalize-across-datasets-2026-11-08]]). Part 2, a joint mixture on z1 / z3 at r1, is running.
