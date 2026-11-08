---
name: fixed-latent-prior-on-state-eval-verdict-2026-11-08
description: "VERDICT 2026-11-08 (igm-refs Parts 7 + 7b): a fixed sparse-PPCA block prior (alpha log w) on state_eval amplifies whatever split its seed carries: +.003-.053 NMI on arms the seed carries, but +.010-.016 on a balanced 10 % null carve too (f2r failure reproduced); a polish only after GAPw, tied-refit seed, alpha 1; held-out alpha unsafe (beta* falls to the grid floor and freezes the seed)"
metadata:
  node_type: memory
  type: project
  originSessionId: 77c8a7fa-5cd0-4244-b81c-f308ec6cd639
  modified: 2026-11-08T12:48:51.016Z
---

Handoff `/mnt/beegfs/elmlund/afan/hetero_tree_20261028/handoff/igm-refs/VERDICT.md`, sections Part 7 and Part 7b; copies in
results/part7/ and results/part7b/. Drivers in `math_reworks/latent/`: `p7_chain.py` (+ `p7_run.sbatch`), `p7_chain2.py` /
`p7_run2.sbatch` (`--qpart`), `p7b_chain.py` / `p7b_run.sbatch` (tied refit seed, `--nthr`, the held-out grid widens at its
edges inside the run). Runs `math_reworks/latent/p7/` and `p7b/<fixture>/b<N>/{A,C,B}/R1-R3`.

- Setup: per arm lead block of the native lp-20 fits, the seed's XD class responsibilities w on the block latent (cross-half),
  passed as SIMPLE_STATE_PRIOR; chains A (no prior), C (alpha 1), B (alpha per round by checkerboard held-out, SIMPLE gets
  alpha*/beta*), rounds 1-3 maxits 1 from forced seed references.
- Part 7 (site chain's seed.txt, so A = lever 3): 10 % fix C > A on 5 / 5 arms (+.016-.053 at R3), labels freeze after R1;
  at 5 % seed.txt is a 1-3 % tail split and alpha 1 anchors it (C < A by .011-.081).
- Part 7b (the sklearn tied refit = v7's seed, balanced; 10 % null, 5 % fix, 5 % null):
  - 5 % fix: the seed carries a1 .181, a3 .164, a4 .135, a2 .097 (a5 .002). C - A +.003-.021 at R3 on those four, -.006 on a5.
    C holds near-true proportions (a3 36 % vs true 35 %); A pulls every split toward 50 %.
  - 10 % null: the tied refit already carves .146-.165 vs the angle median (above lever 3's R3). A holds it (slips on b4, b5);
    C deepens it by +.010-.016 over A on 5 / 5 blocks, so "the null must not keep a carve A drops" FAILS.
  - Held-out calibration: beta* at the base grid's top 0 / 45, bottom 34 / 45. Widened in-run, 9 rounds hit the 1/64 floor:
    SIMPLE alpha 51-68 and B freezes the seed (B < A on 3 / 5 arms at 5 %). alpha* > 0 on all 30 null rounds.
- Traps:
  - The site chain's `bN/seed.txt` is its lopsided INPUT seed (1-5 % minority on null and 5 % blocks); use the tied refit.
  - The held-out beta* sits at a grid edge whatever the seed: the top from tail seeds, the bottom from balanced ones.
  - The framework's alpha (20 A) does not convert to SIMPLE's (10 A) by a fixed factor: alpha 1 moved .42-.70x (tail seeds) and
    .70-1.8x (balanced seeds) the labels the framework predicts.

**Why:** the coordinator asked whether latent posteriors help refinement's state assignment where stage B fails, under the hard
gate ([[flex-frozen-poses-only-no-woodbury-2026-11-06]]); it revisits [[flex-to-refinement-coupling-verdicts-2026-11-06]].
**How to apply:** a fixed prior is a polish, never an acceptance test: it buys as much NMI on a null continuum as on a weak arm.
Use it only on a site GAPw accepted ([[per-site-discreteness-test-gapw-verdict-2026-11-07]]), seeded with the tied refit, at
alpha = 1; never weight it by the held-out alpha. GAPw should test the split the prior will hold. Contrast
[[igm-arm-flex-prior-verdict-2026-11-07]] (focused flex prior: +.003-.006 over references from the same coordinate).
Related: [[lever3-latent-refs-seeding-verdict-2026-11-07]], [[sparse-ppca-pursuit-verdict-2026-11-07]].
