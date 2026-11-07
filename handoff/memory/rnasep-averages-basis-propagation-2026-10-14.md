---
name: rnasep-averages-basis-propagation-2026-10-14
description: user idea tested 2026-10-14: fit flex_pca to the 643 kernel-expanded sub-class averages (works: 16 comps in 77 s, beats the 16 A particle basis at matched rank in the fine band, below the lp 10 particle basis in every shell), then propagate to particles via SIMPLE_COV_COMPOSE (40 columns embedded): the averages-derived per-particle columns pool to 0.75 even/odd but their view-balanced windows carry nothing below 20 A (floor), while inheritance of the sub-class 2D coordinate keeps 2-3x the floor at 20-8 A; propagation by embedding loses the fine part
metadata:
  type: project
---

Scripts: `scripts/cls_expansion_proto/subclass_pool/{render_avgs.py, run_avgflex.sbatch, run_compose_avg.sbatch
(in compose_avg), pool_cache.py, eval_pooled.py}`; runs `rnasep_coupled/avgflex/run_{full,even,odd}` (lp 10, neigs 16,
box 80 @ 2.928 A, no CTF, member-weighted poses) and `rnasep_coupled/compose_avg` (flex2D build; compose of the lp 10
run + the averages run, 40 of 40 columns kept, cut and deconv off, 130 s embedding, cache read by `pool_cache.py`).
Notes sections 17-17.2 of `doc/implementation_notes/coupled_2d3d_math_formulation_2026_10_13.md`.

Gate (even-fitted basis on odd averages): averages 0.50 (null 0.16) vs lp 10 0.60 (0.26) vs lp 10 first-16 0.55 vs
sphere-16 0.53; at 23 A 0.31 vs 0.41 (lp 10) vs 0.15 (sphere-16); union adds ~0.08 margin at 26-20 A. Particle level
(pooled per sub-class): averages columns even/odd 0.69-0.79, tracked F 0.17-0.38 per column (permuted 0.11), class
share 0.26; balanced windows along them 0.27/0.15/0.044/0.03 vs floor 0.19/0.12/0.045/0.03 vs lp 10 columns
0.34/0.18/0.048/0.04 vs the inherited synchronised 2D coordinate 0.23/0.26/0.14/0.14.

**Why:** per particle the fine columns are noise (the compose cut found 1 signal dimension of 40); pooling noisy
per-particle coordinates recovers only their coarse part, whereas assigning each particle its sub-class's coordinate
carries the sub-class SNR intact.

**How to apply:** use averages-fitted bases as bases (cheap, decent), but propagate fine coordinates to particles
by sub-class membership, not by re-embedding; calibrate prior scales if averages columns are composed (set by the
leading variance here). Related: [[rnasep-2d-fine-3d-coarse-complementary-2026-10-14]],
[[batch2-compose-hook-and-worker-names]], [[view-imbalance-inflates-between-state-metric]].
