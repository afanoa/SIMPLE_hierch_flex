---
name: klem-states-need-data-support-not-counts
description: user 2026-10-16: a klem state exists only if the data support it; no resolution-derived or any particle-count floor; node_states2.py --support (cross-half reproducibility >= sim and separability < sim, merge-or-no-state) is the rule, wired into chain8
metadata:
  type: feedback
---

The user rejected a resolution-derived minimum state population ("idk if resolution derived min state population is
valid. We should only have a state if the data supports it").

**Why:** a particle count is a proxy, not evidence; the board already carries the evidence (cross-half halves).

**How to apply:** states = clusters of surviving cells whose mean deviation reproduces between halves at >= --sim (0.7)
and is separable from every other cluster (< --sim); an unsupported cluster merges into its nearest neighbour only if the
merged cluster passes, else no state. Cap (maxk) lifted to 30. Implemented in toy55/node_states2.py (--support), chain
scripts f692_region_one2.sbatch / f692_regions_par2.sh / f692_family_job{,_sr}2.sbatch / f692_chain8{,_sr}.sbatch.
Test 2026-10-16 (minpop_test dirs): real region boards reproduce at 0.98-0.99, criterion merges over-splits -> 10-13
supported states per region board (cap gave 8); toy child boards are mixtures (purity 0.15-0.19) either way, so the
toy ranks it only further down the tree. Related: [[klem-f692-fortran-chain]], [[klem-priority-every-state-before-refinement]].
