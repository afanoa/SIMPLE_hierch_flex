---
name: shared-tree-model-rnasep-direction-2026-11-07
description: "USER DIRECTION 2026-11-07 night: pursue the shared tree model (one sparse site basis, per-site forest, fission-honest tests) until it works on RNase P; working on IgM-like regimes is already 'a great result'. The RNase P blocker is per-particle misfit (amplitude / fit quality), which is being modelled"
metadata:
  node_type: memory
  type: project
  originSessionId: cc23038c-fea9-4b5f-b321-d58eb2942350
  modified: 2026-11-08T04:28:26.677Z
---

User, after the shared latent tree failed RNase P: "We need to pursue this to work on the RNase P. If it works on IgM like
regimes that is a great result." Earlier the same evening: "our current approach lacks elegance which is why I am interested
in a shared model (of course refinements per subset are necessary)", and "port the shared tree latent if it's useful yet".

**Why:** the shared model is the user's preferred design for the method and the paper. It already beats the production tree on
the IgM-like fixtures and holds every null, but on RNase P its splits track per-particle misfit rather than the closed state
([[tree-latent-shared-basis-verdict-2026-11-07]]).
**How to apply:**
- Keep the shared model as the main methods thread.
- The RNase P route runs through the per-particle observation model: the rnasep-misfit session characterises the misfit and
  tests per-particle band scale / noise variance / interpolation-aware noise against the solvent-control floor. The tree-latent
  session runs the misfit-removed rerun (amplitude + nuisances profiled, stratified sat5).
- Port to Fortran once RNase P shows the closed state honestly (r131 corner or the 692k closed lineage), using the forest design in
  `handoff/tree-latent/PORT_DESIGN.md`.
- Present IgM-like success as a result in its own right.
- 2026-11-08: the noise model is found and works ([[rnasep-misfit-noise-spectrum-verdict-2026-11-08]]). With it, the 692k root
  split stops tracking noise but splits along amplitude (a 10 % low-amplitude side), with closed not at level 1. Deeper growth,
  honest maps of both sides, the r131 tree and the native port were running that night.
- 2026-11-08 ~05:30: at r131 neither the latent forest nor the site chain's full-image round isolates closed (7 % at d' ~.9; the
  seeds halve the node, and the halves' maps read the node mean). The production tree finds the corner only at depth 2, after a
  halving carve that enriches closed 1.5x. Open design question: allow an exploratory carve, then an honest test one level down.
- 2026-11-08 ~06:00: the two-level test fails (no site survives in the enriched halves). The production r1 corner comes from stage
  A's iterated full-image 2-state assignment, and its labels read the corner even at fixed poses. But the partition is barely
  enriched in the reference closed rows (.085-.106 vs .069), and the delivered tree's corner leaf reads .93 honestly vs 1.03 on
  own data. So the own-data r131 corner gate is largely selection imprint, and the honest gate (decision 9) must come first
  ([[extremes-fission-verdict-2026-11-07]]).
Related: [[site-coupling-factorial-verdict-2026-11-07]], [[sparse-ppca-pursuit-verdict-2026-11-07]],
[[rnasep-closed-state-must-be-retained]].
- 2026-11-08, user: "subsetting is necessary (with non random seeds derived from the latent) and refinement is necessary. I only
  know this empirically". Coordinator's answer (guide 7.10 thread 9):
  - Mixture-learning theory explains both: project to initialise, then iterate in the full space (Kumar & Kannan 2010).
  - A 7 % minority is a local optimum beside the halving of the largest spread, so it needs a seed inside its basin.
  - Subsetting raises the minority's variance share pi(1-pi).
  - Design: latent proposals (incl. minority tails) -> a full-image round per subset (hard labels on X1) -> an honest fission
    test against matched random groups, with carves as lookahead.
  - Decisive test proposed (r131: tails vs random halves vs fflex seeds; round on X1, read on X2); waiting for the go-ahead.
- 2026-11-08 lineage result ([[r17d-closed-lineage-attribution-2026-11-08]]):
  - The production tree concentrates closed by carving, discrete or not, along a closed-carrying node-latent dimension (r13 / r23
    z1, 93 % of the closed rows), after K-way levels that keep good particles. Per-subset latent refits add only .02-.04.
  - Coordinator's read: the shared model lacks those carves and the conditioning that makes closed lead, not a better basis.
  - Refinement's contribution is unproven: the corner's stage B moved particles along a direction the closed rows don't follow.
  - Per-node noise-model sparse sites on the lineage nodes are running (rnasep-misfit).
- 2026-11-08 user: "Yes, run both". Both tests are read honestly (guide 7.10 thread 9):
  - the carve test, shared model at full 692k (rnasep-misfit);
  - the tail-seed refinement test at r131 on the honest gate's X1 / X2 pair (extremes-repro).
- 2026-11-08: the shared basis carries the production tree's closed-concentrating carves (r13, r23) but has no site at the deep
  nodes (corner); see [[r17d-closed-lineage-attribution-2026-11-08]]. Carve rule caveat: at r23 the carve is on block 2, not the
  leading coordinate.
- 2026-11-08 afternoon results (guide 7.10 thread 9):
  - Carve guards:
    - qtoyflex holds (F1 .571 vs .558); qtoynull holds trivially.
    - The IgM 10 % null fails once: 1 honest pass in 73 tested nodes, z 3.3 vs bar 3.32 under a carved parent.
    - The per-node Sidak bar needs a tree-wide split under carving. Briefed: alpha (n_node / n_root) / (maxdepth + 1), then Sidak
      within the node.
  - r692 root: honest K = 3 (601k / 82k / 8.8k); the table is due ~19:00-20:00 EST.
  - Refinement: production stage A pulls every seed to about half the node, so a minority seed does not survive the round. The
    suspected cause is no class-size prior in the assignment; a per-iteration readout and a code check are running.
- 2026-11-08 ~15:00 (guide 7.10 thread 9):
  - The tree-wide bar alpha (n_node / n_root) / (maxdepth + 1) removes the IgM null's false pass (z 3.33 against 4.11). It keeps
    qtoyflex's root (3.98 against 3.75) and all 6 honest r692 splits through level 2 (tightest r23: 5.93 against 4.62).
  - Stage A has NO class-size prior: prob_neigh labels by a bare argmin over states, and state_eval by an argmax unless
    SIMPLE_STATE_PRIOR is set.
    - Every seed grows to 43-48 % by iteration 6, with poses fixed or free, on r131 and on IgM.
    - The test that separates the missing prior from the small state's noisier reference is one fixed-pose round with log class
      fractions as SIMPLE_STATE_PRIOR (briefed).
- 2026-11-08 evening, class prior at r131 (site2hi):
  - The EM class-fraction prior holds the state at 1.9-3.1k. The zero-prior control grows to 20,427, matching the earlier
    round, so the missing prior causes the absorption.
  - But members churn (39 % of the seed kept), and closed stays at about .09, honest .910 against a bar of .927.
  - Prior = necessary, not sufficient. The next round (briefed) adds the site mask.
- 2026-11-08 evening, refinement test VERDICT (ROUNDS.md; also [[extremes-fission-verdict-2026-11-07]]): the design's full round
  (EM prior + the seed's own site mask, fixed poses) works on IgM (a1 .068 -> .210) but finds nothing at r131.
  - At r131 the closed reference rows and the arm-occupancy maps disagree (blockhi: 2x closed rows, maps below random). That must
    be explained before closed share can judge rounds there.
  - What remains of thread 9 is the 692k carve table.
- 2026-11-08 user: "I think poses matter a lot but if multi state is scrambling membership too much it may be detrimental".
  Master's own multi-state benchmark agrees on scrambling ([[master-multistate-posterior-draw-2026-11-08]]). Proposed lever: hold
  membership fixed while single-state refine3D refines the poses per group (guide 7.10 thread 12).
- 2026-11-08 user: "couldnt you just freeze the memberships and just refine the poses (or only allow a very small number of state
  reassignments)? this is why I wanted to pursue somehow injecting the flex weights into the refinement".
  - The queued test has frozen membership as arm (b), and arm (c) for damped moves: an EM class prior plus an anchor (log w_flex,
    or one-hot on the current state), moves only at posterior >= .95, alternating with single-state pose refinement.
  - It retries [[flex-refinement-integration-verdicts-2026-11-06]] with the class prior and frozen-first membership.
- 2026-11-08 ~17:40, POSE TEST VERDICT (POSETEST.md, 76 core-h): at r131, MEMBERSHIP limits the closed signal, not poses.
  - Frozen membership with solo refinement: 2-5k groups lose contrast like same-size random groups (combined p .15), and the
    loss doesn't grow with closed share. Even IgM's pure a1 tail loses amplitude 172 -> 154: solo refinement of ~5k degrades maps.
  - Damped re-assignment (class prior + anchor + .95): closed shares flat. Flex weights as an anchor are too weak (median .5
    nats) and only freeze the partition.
  - Trap: a class prior plus a fixed .95 margin makes small groups shed members at even odds.
  - Implication: refinement doesn't rescue a weak subset. Purity must come from membership (per-particle evidence, the latent
    rule); the user's "poses matter" is not borne out at r131 for 2-23k groups.
- 2026-11-08 ~18:40, CARVE TREE at 692k (shared basis + noise model, honest splits + carves, maxlev 4) FAILS.
  - Main lineage: 1.26 -> 1.31 -> 1.31 -> 1.27 -> r1212 1.94 % (248k, recall .55), i.e. production's depth-1 r1 three levels late.
  - Best leaf: r1223 3.46 % at recall .055, against production's r131 6.94 % at recall .35.
  - Honest maps above random for r1223 and r1212 only.
- USER CONCLUSION (STANDING): "yeah dont think we can get out of refitting at each node".
  - The shared-basis tree is closed as a replacement for per-node refits on RNase P.
  - Future levers act inside the per-node refit framework: the node's latent fit, the membership rule, the class prior and the
    split test.
- 2026-11-08 user "Yes, run the replays" (guide 7.10 thread 13): a 2 x 2 node-replay factorial on r17d's closed lineage (r131, r13,
  r1, root) in the f2r session (flex-refine-flex-weights): noise model on/off x class prior on/off in the state assignment
  (a new env-switched code path in assign_rows). Readouts: node latent closed AUROC and the children's closed share, recall and
  concentration. Notes in $H/handoff/node-levers/.
- 2026-11-08 ~20:30, node-lever replays at r131 (corrected; guide 7.10 thread 13):
  - Baseline: x1.25 closed concentration, reproducing production's x1.20.
  - Noise model: sharper latent (.698 -> .713), but the split leaves closed (x1.00).
  - Class prior (MAP with current fractions): kills the 2.6 % seed (-3.6 nats), so one leaf.
  - Neither lever concentrates closed at r131. r13 and r1 pending.
  - FLAG: the noise model is default-on on wip. Revert it to opt-in if r13 / r1 also lose closed alignment
    ([[rnasep-misfit-noise-spectrum-verdict-2026-11-08]]).
- r13 replays: baseline x1.74 (production x1.69), noise model x1.71, class prior one leaf. Neither lever concentrates closed.
  - The noise model's r131 loss does not repeat at r13: no stable effect either way.
  - The MAP class prior removes the split at both nodes. No floored or tempered prior: a 1/K floor is no prior, and the EM-prior
    rounds already held minorities flat.
  - r1 pending; the root is held unless r1 shows a gain.
- r1 replays (2026-11-08 ~21:30):
  - Production: site 2 dominant (16.0 vs 4-5.6), closed child r13 x2.12.
  - New-build baseline: near-tie sites (9.482 vs 9.466), x1.10, no concentration.
  - Noise model: production's 3-way structure, x1.90, via a seed-site flip at the near-tie.
  - Class prior: no gain or one leaf.
  - HARD-GATE CONCERN: the wip build's baseline loses r1 -> r13. Checks briefed: arm a on production's own builds (stand-in
    maps vs code), and a second seed for a / b.
- NODE-LEVER VERDICT (2026-11-08 ~22:30, ~490 core-h): neither the noise model nor the MAP class prior concentrates closed at
  r131, r13 or r1, so the root is not run.
  - r1 build check: production's own r17d binary and merge_master2 match the wip baseline on the replay inputs, so this is not a
    regression.
  - r1 seed check: the baseline finds production's 3-way closed split at seed 2 (x2.00) but not at seed 1 (x1.10). The noise
    model found it at both (x1.90 / x1.88): a robustness hint only.
  - TRAP: r1's closed split rests on the site layer's pole choice at near-tie margins (9-10). None of six replays found
    production's dominant margin (16 vs 4-5.6). Closed recovery at r1 is a reproducibility weakness of site discovery and
    ranking.
