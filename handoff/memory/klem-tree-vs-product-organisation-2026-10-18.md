---
name: klem-tree-vs-product-organisation-2026-10-18
description: region-tuple landscape of f692e (tuple_landscape.py): 6 regions x 6-11 states per family give 57k/88k occupied tuples, only 18/6 above 1000 particles, and the populated tuples co-vary across regions (2-4-3-2-*-2 chains): regions are not independent factors; per-particle argmax tuples are too noisy; organisation = low-dimensional global latent + conditional per-region refinement (shallow tree), not a deep tree nor a flat product
metadata:
  type: project
---

User questions 2026-10-18: "is a tree the best organization? would a deeper tree be better?" and "we need an approach
that works on all regimes". Evidence: depth buys purity (closed 0.71 -> 0.94) and costs particles (30-75k -> 2-5k);
leaves across regions re-partition the same particles. Tuple test on f692e: family 1 regions region_1(11)
motion_1(7) motion_2(7) motion_3(8) motion_4(7) motion_5(6) -> 57,273 occupied tuples, 18 >= 1000 ptcls (top
1-1-1-1-1-1: 2390; a chain 2-4-3-2-{0..4}-2 ~ 10k); family 2: 87,686 tuples, 6 >= 1000. So (a) per-particle argmax
labels per region are too noisy to multiply (need confidence-gated or soft tuples), (b) region states co-vary strongly
(a few global coordinates drive them) -> a flat product overcounts; the right shape is a shallow tree: families -> a
low-dimensional global board -> conditional per-region zooms where residual heterogeneity remains, with the parent core
borrowed (compose_state_maps.py) and the data-derived merge within each level.
Unified-regime proposal (same day): region boards seeded from the 2D route (region_contrast_seed.py: per-bin principal
contrast in the projected region mask, Z2 sign sync, 3D mode, particle coordinate -> init file), tested seeded vs random
on arms2/tethered2 regions 1-3 (seed_region.sbatch, jobs 53980343-8). FSC caveat (user): composed / weighted
reconstructions need separate core (parent) and native in-region FSCs; weights for half A from half-B references.

**Seeded region boards, result (jobs 53980343-8, 2026-10-18 ~15:00):** the 2D-route seed works mechanically (100 bins
with a principal in-region contrast explaining 52-68 % of within-bin variance; Z2 sign sync 79-90 % consistent) but its
7-bin coordinate alone carries little truth (arms2 best 0.88-0.91, tethered2 0.66-0.74). As the lattice INIT of the
Fortran region board (7x1 + 1 free, focus = the region) it improves the board's argmax cells: tethered2 best site purity
0.26 / 0.23 / 0.27 seeded vs 0.37 / 0.34 / 0.45 random (regions 1-3); arms2 best 0.70 / 0.69 / 0.66 seeded vs 0.79 / 0.79
/ 0.71 random. First positive transfer of pooled 2D signal into the per-particle board. Next: 2-mode seeds (7x3
lattice), seeded boards inside the chain, and the Fab-mass arms3 set.
Two-contrast seed (7x3 lattice, jobs 53980375-80): tethered2 best site purity 0.21 / 0.19 / 0.24 seeded vs 0.28 / 0.34
/ 0.46 random; arms2 0.65 / 0.64 / 0.66 seeded vs 0.75 / 0.70 / 0.71 random: the gain is consistent across regions and
sets and grows with a second contrast. The unified route (2D-pooled seed -> per-particle region board) is validated on
synthetic data; next = the same on the Fab-mass arms3 set, then inside the chain (needs the direction-bin route run on
the real data: bin-frame stacks exist from the bins80 work).
arms3 (Fab-mass arms) seeded 7x3 region boards (jobs 53984373/5/6): best arm std ratio 0.77 / 0.64 / 0.70 seeded vs
0.78 / 0.75 / 0.79 random (regions 1-3), median 0.88 / 0.81 / 0.86 vs 0.97 / 0.95 / 0.84: the seed helps on every set
tested (tethered2, arms2, arms3); hinged arms remain the hardest regime (best cells ~0.65 = arm spread narrowed to
two thirds). Real-data (f692e family 1) seeded vs random boards pending (jobs 53984090 / 53984370 / 53984371).

**REAL DATA (f692e family 1, 335k particles, seed_rnasep.sbatch, bins80 sub-classes; 2026-10-18 evening):** region_1
board 7x3+1: random init -> 22 cells, mean cross-half reproducibility 0.585, ONE supported state (30k), docked 0/10,
closed max 0.61; seeded init (synchronised 2D contrasts, sign sync 79-87 %) -> reproducibility 0.915-0.932 (16 cells
> 0.9), SIX supported states (30k x4, 25k, 17k), docked 5/10 (mean best corr 0.52), closed max 0.71. The arm-region
(regions_motion/region_1) seeded board: 3 states, docked 2/10, closed 0.67 - but the two masks share the basename
"region_1" so both jobs wrote to seed_rnasep/region_1_seeded (muddled files); rerun with RT = parentdir_basename
(job TAG _v2, seeded + random). Seeding is validated on real data: it goes into the chain's region stage.
CORRECTION (21:0x): because of the basename collision ALL four early real-data boards used the ARM mask
(regions_motion/region_1.mrc); the "region_1" numbers above are arm-region boards. Clean arm-region comparison (TAG
_v2): seeded 3 supported states / docked 2 of 10 / reproducibility 0.93 vs random 1 state / 0 of 10 / 0.59. The
"6 states / 5 of 10" board = arm focus seeded with the contrasts computed inside regions/region_1 (a different region's
seed) - still seeded-beats-random, and the best arm board so far: a seed from a broader region may be the better init
(follow-up). regions/region_1 itself is only now being run (seeded + random, _v2).
Clean region_1 (regions/region_1.mrc, _v2, 2026-10-19 01:xx): seeded -> reproducibility 0.915, 6 supported states,
docked 5/10 (mean 0.52), closed 0.71; random -> 0.716, 3 states (30k, 30k, 6.9k), docked 3/10 (0.47), closed 0.65.
Clean arm region (regions_motion/region_1): seeded 0.93 / 3 states / 2 of 10 vs random 0.59 / 1 state / 0 of 10.
Seeding beats random on both real-data regions by every readout.
