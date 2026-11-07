---
name: klem-speed-plan-board-server
description: speed priority (user 2026-10-16 23:55): the distributed board spends ~17 of 37 min LOADING particles (feedback board 335k: iteration 1 = 1155 s vs 111-152 s after) and every board of a family reloads the same particles (~10 boards per family) -> build a persistent per-family worker pool ("board server") that loads once and serves child/feedback/region/level-2 boards; second: profile the loader (29 particles/s/worker is 10x too slow); third: in-memory exchange
metadata:
  type: project
---

Measured 2026-10-16 (f692c fam1 klem_fb, 10 workers x 8 threads, lps 14..8): total 2229 s; it 1 = 1155 s (worker
load + first E/M), later iterations 111-152 s. Region boards "1028 s" are the same: mostly load. Nodes: csbdevel
48+ cores, 256 GB; sequential board at 32 threads ~112 s/it = the distributed 10-worker rate (exchange floor).
Plan (in order): (1) board server: REQ_NEWBOARD in simple_klem_distr (new odir, focus mask, subset as an active
mask over the loaded particles, lps/sigmas, nfree, init) with workers keeping dat; master runs a list of board
configs; chain region stage = one server job per family (feedback, regions, level 2); reconstructions stay separate.
Load once at the finest kmax (29 covers 20..8 A). Saves ~15 min x ~10 boards per family. (2) loader profile:
klem_load_particles now prints "stack reads X s, FFT + band extraction Y s" per worker (build make_ld, test
run_loadtest.sbatch: 100k subset, 10 workers, 12->8 A). (3) exchange in memory / one multi-node job. (4) dev preset
(chain9 DEV=1) already cuts regions to 3, nl2 1, level-1 lists only, refine 10->8 maxits 5.
Related: [[klem-board-convergence-and-exchange-2026-10-16]], [[klem-f692-fortran-chain]].

Loader profile 00:01 (100k, 10 workers, 12->8 A, job 53914183): whole board 120 s; iteration 1 = 34 s INCLUDING the load, then 14-19 s per iteration. So the 1155 s first iteration of the 335k feedback board was NOT I/O: it was the worker JOBS waiting in the SLURM queue (per-user cap) before loading. Every board re-queues its workers -> the board server (workers launched once per family, kept across boards) removes exactly that; implemented 00:10 (REQ_NEWBOARD/REQ_FLUSH, kmax_load, commander key boards_more = list file of extra board dirs each with board.txt + init file), build make_srv, toy test run_server_test.sbatch (A then B on the same workers vs B alone: cells corr, argmax agreement, times).
