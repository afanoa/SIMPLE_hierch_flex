---
name: prob-neigh-candidate-buffer-corruption
description: not refine3D_states "roind out of range" crash = worker part 1 wrote a corrupted prob_neigh candidate buffer (6 copies of a particle index per record); hael 532517fed 2026-07-12; replay setup in roind_replay
metadata:
  type: project
---

2026-10-01: `not/8_refine3D_states` (default refine3D_states, iter 16, stage 6) hung after worker part 1 died with
`roind is out of range; get_rot` (roind 1936, nrots 512). Chain (addr2line on ~/SIMPLE/build/bin/simple_private_exec):
refine3d_exec -> strategy3D_prob::srch_prob (irot = assgn_map%inpl) -> refine_assignment_continuously (returns
CONT_ROUTE_INVALID leaving inpl) -> assign_ori -> get_rot.
Read-only forensics: assignment_part.dat has 2 of 61,281 records with every field = int 1936/1940 (pinds 25670/25672).
Walking all 10 dist_part_neigh_NN.dat (layout: 4xint64 header | pinds | seed_nrots, shifts(2,n), has_sh(n) | chunks of
[int n][n idx][n x 24-B prob_candidate]): all structurally sound; ONLY part 01 chunk 39 has bad candidates from entry
168,471 on (particle idx 1934; iref=175128 then all-1936). => in-memory overwrite of a thread's candidate buffer in the
worker search (fill_tab_subspace_range, omp default(shared)) before flush. Reader/writer, buffer growth, flush timing,
state_tab path all checked and clean. Code: buffers/stream hael 532517fed (2026-07-12 "refactoring of probabilistic 3D
machinery to reduce peak memory"); loops cyrilr 0a22c537d (2026-06-24); crash-visibility hael 1f2b81f40/4030f044a (Aug 7-10).
Scan scripts ~/roind_scan/{scan,walk}.sh.
Replay: ~/SIMPLE_dbg worktree (f9bbf72ce, Release + -g -fcheck=bounds -fbacktrace, own build) and
/mnt/beegfs/elmlund/afan/roind_replay/replay.sh (RUN=N -> runN/: copies iter-16 inputs = files written 18:25-18:36,
runs prob_tab_neigh part 1 = ptcls 1-80570 nthr=8, then walks both part-01 files). Project at 18:35 holds the iter-16
`sampled` set, so the worker re-derives the same particles. RNG seeding not verified.
