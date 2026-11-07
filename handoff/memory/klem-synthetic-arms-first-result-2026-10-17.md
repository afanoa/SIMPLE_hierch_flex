---
name: klem-synthetic-arms-first-result-2026-10-17
description: IgM-like synthetic benchmark (4 independently hinged arms, 100 configurations, snr 0.0012): the dev chain does NOT localise arm positions (best per-state angle std ratio 0.74, lattice-coordinate/angle |corr| <= 0.48, leaf purity 0.05 vs 100 configs); the challenging regime is real
metadata:
  type: project
---

2026-10-17 18:00, synth_arms (chain9 DEV=1 LOCALDEV, 150k of 200k, MAXREG 3 -> family 1 got 1 region, family 2 three):
- configuration readout (score_states): 11 of 100 configurations are the majority of some leaf, mean leaf purity 0.05
  (chance 0.01) - but a product space of 4 arms is not a set of discrete states, so purity is the wrong metric there.
- per-arm readout (arm_angle_readout.py): angle std over all particles ~23 deg (theta) / 17 (phi); the best leaf state
  only narrows one arm to a ratio 0.74-0.78 (a localised arm would be < 0.4); boards' lattice coordinate correlates with
  one arm's theta at 0.28-0.48. => the tree picks a weak continuous trend per family, no arm-localised states.
- pending: arms_regions_check.py (do the derived motion regions coincide with the arms?); an oracle per-arm focused
  board would separate "wrong regions" from "per-image SNR floor" (the same floor seen on the RNase P arm).
Region check (arms_regions_check.py, per-voxel regression of the 100 maps on the arm angles): the chain's data-derived motion regions DO coincide with single arms (family 1: motion_1 = arm 1 at 0.99 overlap, motion_2/4 = arm 4 at 0.92/0.87, motion_5 = arm 2 at 0.80; family 2: motion_1/2 = arm 3), so the region derivation works; the arm-focused boards still do not localise the angle -> per-image SNR floor. Each synthetic arm is only ~4.4-5.3k voxels (~15k A^3, ~12 kDa) at snr 0.0012 - far smaller than an IgM Fab (50 kDa); the benchmark must be regenerated with Fab-sized arms before it says anything about IgM.

PDH-like (tethered, 11 kDa subunits at 50 % occupancy, 12 A tether jitter, 18:34): no supported region states at all; lattice coordinate vs site occupancy |corr| <= 0.08 on every board -> invisible per image at this mass/SNR; tethered2 (47 kDa subunits, 36 A core) is the fair test. arms2 first launch died: sbatch --export splits commas, so `--arm_radii 26,20,20` arrived as `26`; separator is now x (26x20x20).

**Why:** this is the regime the user cares about (IgM/PDH); the honest number today is a negative.
**How to apply:** judge the arms set by the per-arm readout, never by configuration purity; next lever = arm-specific
region focus and pooled evidence.
