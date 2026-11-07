---
name: price-new-stages-lead-with-lean-scope
description: "User 2026-10-22: 'this sounds rather expensive...' after I briefed a new tree stage at full scope (all 188 site candidates). Price every new stage in core-hours and wall time BEFORE briefing it, and lead with the lean scope (dedupe, subsample-fit + score everyone, refine only what passes)"
metadata:
  node_type: memory
  type: feedback
  originSessionId: 58de2c8f-ad53-44fb-b8a2-05aadff16ec4
  modified: 2026-10-25T04:25:16.542Z
---

I relayed "adopt the factorial per-site layer" to the port as a full-scope spec: site_rank on every discovered candidate, pooling and
refinement. The user's reaction was "this sounds rather expensive...".
- The full scope: about 260 core-h, adding about 1.5-2 h to a 5.4 h tree.
- The lean scope the user picked: about 20-50 core-h, about 30 min.

**Why:**
- The user's target is a 2-3 h tree on 692k, so every added stage competes with that budget.
- Most of a naive scope repeats work: copies of a site at every depth, and root-sized EMs.

**How to apply:**
- Before briefing a new stage, count its units from the run's files (candidates, node sizes) and multiply by a measured per-unit cost.
- Present the lean variant first:
  - run each distinct item once, at the highest node where it appears;
  - fit on a capped subsample and score everyone once;
  - spend refinement only on what passes the existing test.
- Keep the full variant as an option, priced.

2026-10-23 (coordinator, after I quoted the step-0 junk replay round at 1.5 core-h and it cost ~18):
- If a price grows by more than about 2x after approval, send the new number BEFORE launching.
- Price replays and trees from sacct of the last same-kind run (ElapsedRaw x AllocCPUS), not from memory.
- Measured 2026-10-23: a junk replay (stree_junkval24, 32 cores) costs ~2.5 core-h per toy/CAHRA node and 10.6 at the
  692k RNase P root. A p27-style toy tree (8-core master + nparts=4 x nthr=6) burns ~28 core-h per wall-hour; a full
  cdisc/IgM tree runs more than 4 h, ~150-220 core-h.

2026-10-24: priced 21 split reruns at 100-200 shuffles at 4.5 core-h; they cost 18.6 (cap 5).
- I had scaled small IgM nodes' 20->200-shuffle ratio (about 3.3x) onto 50-117k cdisc nodes. Those splits had run at 8-32 threads, and a 117k node at 100 shuffles on 6 threads runs for tens of minutes.
- Price split reruns per node from a measured per-particle-per-shuffle cost.
- Run the largest node first as a probe, or watch the first runs mid-job and re-price before the rest.
- The merge-up scales with the leaf count. m33 at 692k with 81 leaves took 1 h 30 min on 32 cores (48 core-h); the 12-leaf toys took about 4.
  The capped m33 finish cost 71.6 core-h against a 45-50 estimate.
- STANDING (coordinator, after that overrun):
  - set every sbatch's -t so that cores x time is about the quoted price, and the job can't run long unwatched;
  - if a job heads past about 2x its price, stop it by id and re-price before going on.
- EXCEPTION (user 2026-10-25, on an 18 h -t for the RNase P solve2D round 2: "why is the time limit 18h thats stupid"):
  - the 2x -t rule is for tests, replays and diagnostics only;
  - production pipeline runs (solve2D, abinitio, full-data chains) keep generous limits, as in round 1:
    master -t 5-00:00:00, workers SIMPLE_DEFAULT_PARTITION_TIME=true;
  - a timeout there throws away hours of work, and the user watches them anyway.
- Price distributed SIMPLE runs from the WORKER jobs' sacct (sum of ElapsedRaw x AllocCPUS), not from cores x master wall.
  The workers hold cores only about 3 of every 16 min per iteration. RNase P solve2D round 1 was about 2,100 core-h, not the
  5,900 I had quoted from 384 cores x 15.2 h.

Related: [[port-tree-site-layer-recovers-closed-2026-10-22]], [[keep-test-loops-minimal]], [[no-4m-synthetic-tests-too-expensive]].
