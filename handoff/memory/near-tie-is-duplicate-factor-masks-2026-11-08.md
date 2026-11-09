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
- IgM panel FINAL: joint seeding recovers no arm at the IgM roots, and no Gaussian null (A, bin-wise, f2r) is null-safe there
  (the 5 % null gets K* 6). The split test returns one group only through its effect floor.
  - Panel: qtoy2 +, qtoyflex -, IgM -; RNase P r1 and OLE pending.
  - Joint seeding is not yet a general method.
- r1 REFINED FROM JOINT SEEDS (decisive, seeds 1-3): POSITIVE.
  - Six children with sizes within 1 % across seeds, all passing the tree's own tests.
  - 3 / 3 closed hits (x2.34-2.47, recall .23-.24); single-site seeding gave 1 / 3 at the same seeds.
  - Cross-seed ARI .68 (single-site .31 across basins); the three closed children share 35k at .052 (x2.68).
  - Refinement converges on the same 6-way partition from every seed.
  - Costs: closed spreads over 3 children (lower top-child recall than production's r13); on the toys it loses the holo.
  - Next: the panel (IgM and OLE refined from joint seeds), then a subtree method test.
- OLE panel (root): joint seeding reproduces WORSE than single-site seeding (ARI .34-.67 vs .49-.95). Cause: mostly XD's K
  choice carried into the mixture (seeds 1 / 4: identical latents, r 1.00, but XD K 4 vs 5); seed 3's latent does differ.
  A raw-latent (pre-XD) variant of f2r's rule is running.
- XD (tree-latent): reproducible within an arm, but its components are soft (responsibility .5-.67). Every null splits under the
  Gaussian-null rule (non-Gaussianity, not discreteness).
- 692k ROOT, replay seed 1 (2026-11-09): region 8's two masks tie (24.9 vs 22.5); 6 children, ARI .19 with production. Its two
  closed-rich children hold 2.24 % at recall .72 (production r1: 1.94 % at .56). Seeds 2 / 3 (same build) pending.
- PANEL (2026-11-09): IgM refinement from joint seeds is null-safe but keeps NO root split (1 leaf vs stock 2) and its stage B
  agrees across seeds only .27-.39 (nulls alike; RNase P r1 .68): NEGATIVE. OLE depth 1 has no kept pole at all (fallback
  seeds, XD K near-ties); joint groups .21 / .16 vs refined split .33 / .08; refinement from joint seeds running.
- XD-REPRO FINAL: within one latent, hard labels reproduce as well as XD's components (ARI .98-.99; flex deconv labels
  .95-1.00), so flattening isn't the gap; the split's seed dependence is downstream (site / pole choice, stage-B seeding).
  flex's own deconv labels give a reproducible r1 group at x1.9-2.0, recall .28-.33, on every seed and both latents.
- 692k ROOT seed 2 (same build as seed 1): region 9's masks tie (21.1 vs 19.5); split ARI .08 with production, .125 with seed 1.
  Closed child 148k at 3.67 % (x2.91), recall .62: closed concentration x1.54-2.91 depending on the draw. Instability starts at
  the root on wip; production's root had a dominant site (63.3 vs 43.0).
- OLE depth 1: refined joint children disagree across independent trees (.03 / .07), confounded by different rows; same-row
  replays running.
- Root tied pairs are less duplicate than r1's: seed 1's region-8 masks Jaccard .549 (centroids 12.8 A apart), seed 2's region-9
  masks .365 (29.3 A); r1's pair .67-.68 at ~6 A. OLE depth 1: same-row replays reproduce the single-site split almost exactly,
  so OLE's cross-tree disagreement comes from the root's divergence.
- OLE PANEL FINAL (244 core-h): root XD K flips (4 / 6 / 5 / 5, gaps <= .0007) change the pole seeds; depth-1 nodes replay
  exactly on the same rows, so OLE diverges at the root. Joint seeding: root .34-.67 vs single-site .49-.95; depth 1 no gain.
- PANEL VERDICT 2026-11-09: joint seeding + RNase P r1, + qtoy2; - qtoyflex, - IgM, - OLE. Not general: don't adopt it or run
  the r1 subtree with it. Common cause: flat (near-tied) discrete choices at the top (site choice, XD K).
- 692k ROOT FINAL (3 wip seeds + production, `handoff/root-repro/ROOTREPRO.md`): split ARI .08-.37; the chosen site changes with
  the seed (8f1 / 9f2 / 9f2; production 7f1) and the site tables differ. Seeds 2 and 3 chose the same site and share seeds at
  .76, yet their splits agree only at .37: at the root, stage A / B moves memberships too. Closed child x1.54 (production) to
  x2.71-2.91 (seeds 2, 3).
- ROOT REFINEMENT READ: seeds carry only x1.30-1.37; stage A / B lifts seeds 2 and 3 to x2.71-2.91 at the same recall while
  moving 50-66 % of particles. Seeds 2 vs 3: seed groups Jaccard .86, closed children .59 (closed rows .76). Freezing
  membership would keep .86 but lose the concentration. Consensus-of-replicates read pending.
- ROOT CONSENSUS (zero compute): seeds 2 and 3's closed children intersect in 110k at 4.13 % (x3.28), recall .52, purer than
  either (x2.91 / x2.71); the one-run-only particles sit at x1.47. At the root this matches production's r13 (4.11 %, .38).
  Candidate: stable site choice + two replicate refinements + consensus. Pair-to-pair reproducibility is untested.
- TIE-ROBUST RULE, STAGE A (sparse-ppca2, `handoff/ole-panel/NEAR_TIES.md`): NEGATIVE. A 2-SE parsimony rule leaves K
  disagreeing (OLE root 4 / 5 / 5 / 4, RNase P root 8 / 7 / 9: shallow, non-monotone curves). The site REGION already agrees
  in every replay, and discovery is stable; only the f1 / f2 masks tie. Most disagreement sits in the pole seed set (via XD K)
  and in refinement. Next: consensus of replicate refinements (pair-to-pair test, panel read).
