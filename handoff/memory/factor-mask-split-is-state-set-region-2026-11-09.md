---
name: factor-mask-split-is-state-set-region-2026-11-09
description: "VERDICT: the 692k root's f1/f2 split depends on whether the seed's 8-state set puts f2's territory into the region core (factor-2 share .23-.31 split vs .06-.13 whole; the replication z passes every factor); f2 = inner sub-block ~40 A from the r131 closed block (which sits on f1 / the whole region); the flex basis finds f2 in all 7 seeds (J .80-.96); a deterministic region needs restarts, ~80-105 core-h per 692k node"
metadata:
  node_type: memory
  type: project
  originSessionId: 6f195b76-225b-466c-a432-7fc8360148e8
  modified: 2026-11-09T11:18:43.367Z
---

Overnight zero-compute read (coordinator, 2026-11-09).
- Notes: `hetero_tree_20261028/handoff/ole-panel/FACTOR_MASKS.md`.
- Scripts: `ole_panel/scripts/{fmask,fbasis,fbasis2,fenergy,frame}.py`, pure Python on copied masks.
- Outputs: `ole_panel/factor_masks/`.

**What decides the split.**
- `discover_sites` keeps factor f when cross-half z >= 3 AND share >= DFACT_SHARE 0.2. Half-map noise is ~0 at 692k: factors 2-4 all reproduce at r .98-1.00, so the share decides.
- The share is not near the threshold: split runs .23-.31, whole runs .06-.13.
- **The region itself differs.** Split runs' core (vs > 0.2 max over the 8 stA0 states) reaches inward over f2. Unsplit runs' region is the outer part only, holding <= 22 % of f2's core.
- **Source of the state-set difference.** The diffusion k-center has no RNG (stride subsample, farthest-point). Seeds still get different rim states, because the embedding differs slightly between seeds and stA0 reshuffles the states.

**f1 vs f2.** site_rank margin = min(repro z, distinct z) of the strongest pole. It is a near-tie in seeds 1 / 2 / 5 (gaps 2.4 / 1.6 / 1.3, all < 2.83); seed 3's gap is 8.0. f2's pole is stable: 71.5-75.5k members, z 21-24.

**What f2 is.**
- **Location:** an inner sub-block at (-19, 1, 27) A from the box centre, radius 31-35 A, 2.5-5.2k A^3. f1 sits at (-8, -19, 36); the whole region at (-5, -29, 41).
- **Not the r131 closed block.** The native nvarB blocks' energy sits on f1 / the whole region (core share .23-.40, 13-29 A away). f2's core holds <= 5 % of it, 32-56 A away. f2 is a correlate of closed, not its location.

**State-free check.** fflex's PPCA basis, factorised inside the replays' territory:
- factor 2 = f2 in all 7 seeds (share .24; J .80-.96 across seeds; .89-.95 between fit A and fit B);
- inside the unsplit runs' own region it misses f2 (J .10);
- the basis variance map does not delimit the region.

**Why:** the user's principle "never choose on a difference below its own noise" ([[landscape-must-be-reproducible-and-splits-data-derived]]). The noise that matters here is the state set, which no within-run statistic sees.

**How to apply:**
- A deterministic site layer needs a region step pooled over state-set restarts: ~13 core-h each at 692k; R 6-8 keeps f2 at P .94-.99, if f2 appears in 4 of 7 sets.
- The factor step can use particle-half replication of the flex basis, at 0 core-h.
- Keep tied sites (f1 / f2) instead of choosing.
- Related: [[near-tie-rule-verdict-2026-11-09]], [[two-run-consensus-verdict-2026-11-09]], [[rnasep-misfit-noise-spectrum-verdict-2026-11-08]].
