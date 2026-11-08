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
