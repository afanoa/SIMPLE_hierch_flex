---
name: klem-chain-speed-audit-2026-10-19
description: "klem chain wall time ~75-85 min per tree at 50k and 200k (real 692k 2.5 h); fixes applied 2026-10-19: filt_mode=none in refine3D_states (27->13 min), SERVER=1 PARMAX=2 MC=16 chain10 defaults; final rec not worth skipping"
metadata:
  type: project
---

Measured 2026-10-19 from logs/sacct (chain9, slurm mode, no queue wait):
- 50k tree 73 min, 200k tree 84 min, real 692k (f692e, 2026-10-17) 2.5 h families -> ALL DONE. The "15 h" I
  quoted for the real chain was wrong.
- refine3D_states per family (10 states, 6 its): 22 min at 19k particles, 27 min at 100k; the parts' search is
  seconds per iteration; the time is serial per-state master work (RESOLUTION/NU filter ~27 s x 10 states at
  init, every iteration, and ~9 min of final per-state post-processing).
- Region job ~18 min: focal board 563 s of which it-1 M-step 338 s = worker job launch + load (later its 3-8 s);
  level-2 board 244 s (same); recs. 6 regions at PARMAX 3 -> 28-36 min. NL2=0 saves ~8 min per region.
- Consolidation count_states.py 10-11 min, single process. Master held 48 idle cores (67 of ~330 core-h per
  200k tree).
Fixes by payoff: (1) drop/parallelise the per-state NU+resolution post-processing in the states step (~15-20
min per tree); (2) board server (SERVER=1, REQ_NEWBOARD, implemented 10-17, run_server_test NEVER RAN) to keep
workers loaded across the ~13 boards per family (~10 min per region); (3) PARMAX 6 (regions 28 -> 18 min);
(4) count_states out of the chain; (5) MC=16 masters. Target ~35 min per tree.
IgM/PDH status: tethered2 (47 kDa subunits) localised by the tree; arms up to Fab mass still weak (std ratio
0.66-0.69); seeded region boards on arms3 (jobs 53984373/5/6) = no gain over random (best 0.64-0.77 vs
0.75-0.79). Related: [[klem-speed-plan-board-server]], [[klem-synthetic-benchmarks-arms-tethered]],
[[batch2-400-core-limit-and-worker-timeouts]].

Applied 2026-10-19 (chain10 + prep defaults, scripts under $W): refine3D_states with REFINE_EXTRA default
filt_mode=none: 785 s (13 min) vs 22-27 min, volassemble it0 175 -> 86 s (nonuniform_filter 0); the final
calc_final_rec is only ~2 min, so no kernel change. f692_chain10.sbatch exports SERVER=1 PARMAX=2 (region boards
of a family queue on one persistent 10-worker server, ~112 cores per family); synth_tree_prep_full.sbatch defaults
MC=16, CHAIN=chain10. Server verified identical cells to standalone boards (corr 1.000). Not yet timed end to end.
