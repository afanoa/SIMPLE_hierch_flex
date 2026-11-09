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
- NODE TRACE (r1 replays): divergence starts at SEEDING. Hit vs miss runs share seeds at ARI .20-.34; within a basin the seeds
  agree at .68-.95. Stage A pulls the two basins to ~.5, and the split separates them again (.3). Matching seeds give matching
  splits (.63-.74).
  - Refinement's membership switching is not the source of run-to-run instability (the user's hypothesis); it keeps each run in
    its seed's basin.
  - The fix belongs in site discovery and pole selection.
- Full trace switch rates (seed -> stage B): 26 % with two seeds (misses), 53 % with three (hits). The small z3hi seed
  (4-6k) grows to ~75k = the closed-enriched child; closed rows switch at .28 in every run.
  - Refinement's switching is reproducible given the seeds, and it is HOW r1's closed state forms.
  - Never remove refinement's re-assignment at a node unless the latent mixture forms the group itself.
  - OLE's October outputs are gone, so tracing OLE needs a rerun.
- JOINT-SEEDING PROTOTYPE (part 2): a Gaussian mixture in the latent axes the tied sites use (z1, z3), with K by cross-half
  held-out likelihood.
  - r1 miss seed: holds a closed-enriched group (x1.66, recall .20), reproducible at seed 2 (x1.67). After fixed-pose maps,
    the split test keeps 4 distinct children (one 106k at x1.69, recall .40); single-site seeding gave 2 at x1.10.
  - qtoy2 root: 6 distinct conformers against the tree's 5 children (it recovers closed conformer 18, which the tree merged;
    loses holo 101 at seed level).
  - The site endpoint posteriors are the wrong coordinate; the latent axes are the right one.
  - TRAP: not null-safe. On qtoynull, K 2 wins by .003 nats per particle and passes the split test (z 4.97). K needs a
    null-calibrated floor or an effect-size gate.
  - Next: a null-safe K, then refinement (stage A / B) seeded from the joint groups, across seeds plus qtoy2 and the nulls.
- 2026-11-08 user "Yes, run it": refined joint-seeding test (~120 core-h, f2r session, guide 7.10 thread 14).
  - Null-safe K: held-out gain above the 95th percentile of a Gaussian null with per-particle noise, plus the effect gate.
  - Stage A / B seeded from the joint groups at r1 (seeds 1-3), qtoy2, qtoyflex and the nulls, against single-site seeding.
- RELIABILITY VERDICT (r1, 12 runs, hit = closed child >= 1.5x): noise off 2 / 6, noise on 4 / 6 (Fisher p .57). Neither is
  reliable, so single-site seeding must be redesigned (user rule [[replicate-or-redesign-positives-only-2026-11-08]]).
  - Upside: one seed split 4-way with a 58k child at 6.95 % closed (x3.59, recall .47), r131-level in one step.
  - A trace over all 14 runs (55094636) places the divergence.
- 14-RUN TRACE: site discovery's candidate set varies with the seed (b hits {2f1,1f1,1f2,2f2,5}; b misses {1f1,1f2,4}), and
  the same mask can seed a hit or a miss. The pole labels decide.
  - The divergence is stochastic site discovery plus pole labelling, deeper than the duplicate-mask tie alone.
  - The latent axes are reproducible, so seeding from them (joint / XD) sidesteps it.
- Step 2, refinement from the joint seeds:
  - micnull: the locality guard folds it to one leaf (null caught end to end).
  - qtoy2: 6 of 8 conformers recovered against the tree's 4 of 8, at lower purity; the holo merges.
  - qtoyflex: WORSE (1 of 3; the holo seed swells 5.4k -> 23k).
  - r1: joint seeds IDENTICAL across seeds (K 6; two hit groups x1.66 / x1.56); refinements at seeds 1-3 running.
  - Read: joint seeding fixes seeding reproducibility and adds states, but un-priored refinement swelling still drains small
    states.
