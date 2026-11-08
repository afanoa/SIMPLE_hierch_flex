---
name: fixed-latent-prior-on-state-eval-verdict-2026-11-08
description: "VERDICT 2026-11-08 (Part 7 of igm-refs): sparse-PPCA block posteriors (precision A) as a fixed alpha log w prior on state_eval help only where the seed carries the arm (IgM 10 %: +.012-.053 NMI per arm over lever 3, no drift); 5 % site-chain seed.txt is a 1-3 % tail split and alpha 1 anchors it (-.011..-.081); no stop of null carving; held-out alpha > 0 on nulls too (not a null detector)"
metadata:
  node_type: memory
  type: project
  originSessionId: 77c8a7fa-5cd0-4244-b81c-f308ec6cd639
  modified: 2026-11-08T09:14:56.272Z
---

Handoff `/mnt/beegfs/elmlund/afan/hetero_tree_20261028/handoff/igm-refs/VERDICT.md`, section Part 7; copies in results/part7/.
Driver `math_reworks/latent/p7_chain.py` (+ `p7_run.sbatch`); `p7_chain2.py` / `p7_run2.sbatch` take `--qpart` (norm).
Runs `math_reworks/latent/p7/<x10fix|x10null|x5fix|x5null>/b<N>/{A,C,B}/R1-R3`.

- Setup: per arm lead block of the native lp-20 fits, the seed's XD class responsibilities w on the block latent (cross-half),
  passed as SIMPLE_STATE_PRIOR; chains A (no prior = lever 3), C (alpha 1), B (alpha per round by checkerboard held-out,
  SIMPLE gets alpha*/beta*), rounds 1-3 maxits 1 from forced seed references.
- 10 % fix (seeds .16-.39): C > A on 5 / 5 arms every round (R1 +.012-.035, R3 +.016-.053); C's labels freeze after R1
  (agreement .98-.997), so lever 3's drift (a4 .329 -> .306) is gone.
- 5 % fix: alpha 1 keeps the seed's tail split (a1 minority 3 % -> 8 %; A rebalances to 41 %): C < A on 5 / 5 arms, every round.
- Nulls: round-1 carving unchanged by the prior; alpha 1 slows later rounds only by freezing the seed.
- Held-out alpha: SIMPLE alpha .125-.25 on all 10 null blocks (gain +.004-.039 nats/ptcl), the same as on the 5 % arms. A prior
  along a real continuum predicts held-out data; f2r's alpha-0 nulls had no arm motion. B is the safe weight
  (B - A at 5 % within -.012 / +.004).
- Traps:
  - The site chain's `bN/seed.txt` is its lopsided INPUT seed: a 1-5 % minority on every null and 5 % block. Its round-1
    `labels.txt` (22-31 % at 5 %) and the sklearn tied refit (5 % fix NMI .10-.18) are balanced.
  - beta* sits at the grid top (2.0) in 39 / 45 null and 5 % rounds, so SIMPLE's alpha there is set by the grid.
  - At 20 A vs 10 A, SIMPLE moves .42-.70 of the labels the framework predicts for alpha 1.

**Why:** the coordinator asked whether latent posteriors help refinement's state assignment where stage B fails, under the hard
gate ([[flex-frozen-poses-only-no-woodbury-2026-11-06]]); it revisits [[flex-to-refinement-coupling-verdicts-2026-11-06]].
**How to apply:** use a fixed prior only as a polish on a GAPw-accepted site whose seed carries the arm
([[per-site-discreteness-test-gapw-verdict-2026-11-07]]); weight it by the held-out alpha, never as a null test. Open, proposed not
run: the 10 % null + 5 % fixtures with the sklearn tied refit seed (`--seeds ''`). Contrast
[[igm-arm-flex-prior-verdict-2026-11-07]] (focused flex prior: +.003-.006 over references from the same coordinate).
Related: [[lever3-latent-refs-seeding-verdict-2026-11-07]], [[sparse-ppca-pursuit-verdict-2026-11-07]].
