---
name: node-levers-replay-verdict-2026-11-08
description: "VERDICT (2026-11-08): 2x2 replays of production's closed lineage (r131, r13, r1): neither the flex noise model nor a class prior in the state assignment concentrates the closed state; the class prior (previous-iteration counts, no floor) empties the seed minority and removes the split; r1's closed concentration rests on the site layer's pole choice at near-tie margins (same on r17d's binary and on wip): an input/site-layer sensitivity, not a regression"
metadata:
  type: project
---

Node-levers replays (handoff/node-levers/VERDICT.md, STATUS.md), new build 84d4bb75d = wip 918215446 + class prior
(SIMPLE_STATE_CLASS_PRIOR=yes; branch nodelevers, not on wip), fast gate 15/15; nodes replayed by resuming the tree with a one-node
table at the node's production depth (parent maps regenerated: recB0 halves exactly, the recB map from the stage-B project).
- Baseline (noise off / prior off) reproduces production at r131 (latent .698 vs .699, closed child x1.25 / .226 vs x1.20 / .209)
  and r13 (x1.74 / .355 vs x1.69 / .353). At r1 it misses production's 3-way split at seed 1 (x1.10) and finds it at seed 2
  (x2.00 / .395 vs production x2.12 / .380): the site ranking sits at near-tie pole margins (9.4-9.9) where production had a
  dominant pole (16.0); production's own r17d binary on the same replay inputs misses it too (x1.08) and merge_master2 equals wip.
- Flex noise model: latent best dimension +.01-.02; children unstable (r131 x1.00 vs baseline x1.25; r13 equal; r1 x1.90 / x1.88 at
  both seeds, never above the baseline's hit or production). Not a lever.
- Class prior as specified (log of the previous iteration's fractions on the nat scale, no floor): a 2.6 % seed is charged 3.6 nats
  per particle at the first stage-A iteration and empties (4-29 particles left), the re-test fails, one leaf at all three nodes. A
  1/K floor is no prior for K = 2 and a tempered prior sits between dying and swelling (coordinator): not run.
- Reliability (6 seeds at r1): arm a 2/6 hits, arm b 4/6 (Fisher p .57); the candidate masks (3-6) and the chosen pole labels vary
  with the seed; hit-hit seed ARI .72-.98, hit-miss .16-.26, stage B follows the seeds; refinement churn 2-5 % per late iteration.
  The instability is site discovery + pole seeding, not refinement. Seed 6 arm a: 4-way split with a 58k child at 6.95 % closed.
- Joint seeding (one mixture over the tied sites' latent coordinates) hits at r1 on both seeds (x1.66-1.69) and gives qtoy2 six
  children / six conformers, but splits the qtoynull root (not null-safe as designed).
- Responsibilities: production's own prob_state_evidence.txt at every node, sharper and more closed-informative than hard labels
  ([[production-tree-has-per-node-state-evidence-files]]); not readable for wip ([[wip-refine3d-redraws-poses-each-iteration-2026-11-08]]).
- X1 honest contrast at r131 (2026-11-09): arms a/b on the fission copy MISS (x1.09 / 1.06; latent AUROC .51 vs .61 on X); the
  children's own-data closed-arm difference (.084) vanishes on the honest X2 maps (.000): self-selection on the real lineage.
  The X1 arms cost ~48 core-h each (the flex fit runs 15x slower on the noisier copy).
- Cost: r131 ~5 core-h per arm, r13 ~8, r1 ~25 (production's node steps 11 / 19 / 72 with their fuller site layer); X1 arms at
  r131 run ~15x slower in the flex fit (paired PCG probe on the noisier copy, ~50 min per probe iteration).

**Why:** the user approved the factorial to find which per-node lever concentrates closed; the answer is neither, and the r1
finding redirects attention to the site layer's near-tie fragility.

**How to apply:** do not port the class prior; do not count the noise model as a closed lever; treat the r1 site choice (and any
node whose pole margins tie within ~0.2) as a reproducibility risk and test it with seeds, not with new levers; the fix is a seed-independent site choice or a
null-safe joint seeding, not a prior inside refinement. Replay vehicle and
traps: [[never-mirror-into-reused-roots-with-production-symlinks]]; related [[stageb-oracle-references-limit-2026-11-07]],
[[oversplit-junk-ref-merge-verdict-2026-11-07]].
- r131 X1 HONEST CONTRAST (2026-11-09, ~110 core-h): splits selected on X1 (arms a / b) miss closed (x1.09 / x0.89; x0.95 / x1.06).
  Rebuilt on X2, the children's closed-arm occupancy is identical (difference .000 / .002) where X1's own maps showed .084 / .011.
  The r131 split's own-data closed contrast is imprint. Caveat: X1's latent carries less closed signal (AUROC .51-.63 vs .61-.70).
