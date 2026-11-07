---
name: flex-pca-downstream-scoring-protocol
description: "How to score flex_pca state maps downstream (spread + cmin, matched rank, occupancy-selected) and the three traps that make the comparison lie"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 991559ff-d34f-451c-ad5a-ac1ba95b88c7
  modified: 2026-08-10T17:15:09.004Z
---

Cross-halfset basis reproducibility ranks BASES; it cannot decide whether a change helps the
delivered states. For that, score the state maps with TWO numbers — never one:

- `spread` = mean pairwise correlation of the MOST-SPREAD triple over all C(n,3) triples, computed
  inside the soft-mask sphere. LOWER = states genuinely differ.
- `cmin` = the WORST state's fraction of positive mass inside its densest 1% of voxels.
  HIGHER = each state is still a molecule.

`spread` alone is maximised by a degenerate blob, which decorrelates from everything; `cmin` is the
guard. Scorer: `state_score.py` (rewritten 2026-08-10 from the protocol in
[[flex-pca-ab-testing-traps]]; its ABSOLUTE values do not match the older scorer's, so compare runs
within one session, never against remembered numbers).

**Three traps, all hit on 2026-08-10:**

1. **Small n is not a measurement.** At npreimages=2 (3 states) there is ONE triple. Per-shell vs
   default looked like a per-shell win (spread 0.8913 vs 0.9086); at npreimages=8 it reversed
   decisively (0.4300 vs 0.3583). Use npreimages=8+ before concluding anything.
2. **Match rank before comparing.** `spread` is a MINIMUM over C(n,3), so it falls with n on
   combinatorics alone. Two runs can deliver the same NUMBER of map files while differing in
   substantive states — per-shell delivered 8 maps but occupancies 206/38978/321/3756/256/36482/
   19817/184, i.e. ~3-4 real states against default's ~6-7. A near-empty state's map is
   prior-dominated, so letting it into the triples measures the prior. Select the top-m by
   OCCUPANCY (label-invariant), not by state index.
3. **Reproducibility is not conformationality.** Cross-halfset agreement rising with N shows content
   is SYSTEMATIC, not that it is conformational — both halves share the same imaging physics, so a
   resolution-dependent CTF/alignment residual reproduces identically. Only the downstream metric
   separates the two.

**Why:** every one of these makes a change look good that is not.
**How to apply:** run at 100k with npreimages=8, score spread+cmin, match rank by occupancy, and
report both numbers and the occupancy vector.
