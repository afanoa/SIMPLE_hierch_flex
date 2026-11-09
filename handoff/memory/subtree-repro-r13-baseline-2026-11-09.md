---
name: subtree-repro-r13-baseline-2026-11-09
description: "BASELINE for the joint-seeding fix: r17d's r13 subtree to leaves at 3 seeds on merge_master2 -> 0 of 8 closed leaves reproduce (best Jaccard .13-.22, diff-map r .46-.78); r13 split reproduces (.94), r131's does not; closed-leaf unions overlap .41-.55; the all-seed intersection is 13.6k at x2.91"
metadata:
  node_type: memory
  type: project
  originSessionId: 53e38d89-ea70-4683-80d0-2600450cc355
  modified: 2026-11-09T05:13:51.182Z
---

Run: the extremes-repro session, `handoff/subtree-repro/VERDICT.md`, 127 core-h.
- **Setup.** r13 resumed from production's saved entry rows and poses, production's keys, maxdepth 4, SIMPLE_SEED 1-3, run to
  the tree's own leaves. Closed leaf = concentration >= 1.5x r13's .041; a match = Jaccard >= .5 or difference-map r >= .9.
  Both rules were set before the runs.
- **Closed leaves don't reproduce.** Seed 1 has 3 closed leaves, seed 2 has 2 and seed 3 has 3. None matches in the other
  seeds (best Jaccard .13-.22, difference-map r .46-.78), and none matches production's 2 either.
  - Each seed's best leaf reaches x2.17 / x2.98 / x2.73 (production x2.03).
  - Each seed's closed leaves recover 26-30 % of the 8,720 closed rows.
- **The divergence is r131's split.** r13's split reproduces: x1.73-1.74, partition agreement .94. r131's closed-leaning
  halves overlap at Jaccard .22-.39. All seeds also split r1312 and r1321, which production kept whole.
- **The closed side is partly stable.**
  - The unions of each run's closed leaves overlap at Jaccard .41-.55; their closed rows overlap at .64-.71.
  - The intersection over the 3 seeds is 13,630 particles at .119 (x2.91), recall .187.

**Why:** the coordinator's baseline for judging the joint-seeding fix under the user's replicate-or-redesign rule
([[replicate-or-redesign-positives-only-2026-11-08]]).
**How to apply:**
- Judge any fix on the same vehicle (`handoff/subtree-repro/scripts/st_subtree.sbatch`, readout `st_readout.py`) with the same
  rules.
- A fix must make closed leaves (or confident cores) match across seeds. Matching closed-side unions is not enough.
- Traps of the vehicle:
  - Fresh roots holding copies only. The site layer writes into ancestor folders, and children reuse production's node names,
    so a mirrored root would write into production ([[never-mirror-into-reused-roots-with-production-symlinks]]).
  - Carry no site table over.
  - Stop the master before the whole-tree finishing.
  - Raw map r is .92-.98 for any two leaves; use difference maps.

Related: [[confident-cores-from-refinement-posterior-2026-11-09]], [[parallel-slate-2026-11-08]],
[[extremes-fission-verdict-2026-11-07]].
