---
name: klem-focus-step-was-the-board-floor
description: the ~150 s/iteration floor of focused klem boards was the master's single-threaded real-space focus masking (102 box-160 FFT round trips); parallel over cells it is 5 s; and the served board list (fboard.sh SERVER_LIST) design
metadata:
  type: project
---

Finding (2026-10-17): focused Fortran boards (child, feedback, region) took 100-180 s per iteration at every band and every n,
while the unfocused parent (692k/10 workers) ran 12-47 s and the 100k load test 14-19 s. The cost was `apply_focus_d` on the
master: ncl x 2 halves = 102 full-box (160^3) ifft/mask/fft round trips, single-threaded. Parallelised over cells with one
preallocated image per thread (FFTW plan creation is not thread-safe, execution is): "solve + focus 5.2 s" on the 100k focused
load test, iterations ~20 s. The per-iteration timing line (`>>> it N timing: M-step round, solve + focus, E-step round`) is in the log.

Served boards: `fboard.sh` SERVER_START=1 + SERVER_LIST=<file> starts the master in the background (board.txt keys
`server_list`, `kmax_load 29`), pid in <file>.pid; SERVER_LIST alone = served mode (append the board dir to the list, wait for
<dir>/BOARD_DONE, <file>.dead aborts); `echo STOP >> list` ends the service. Served boards must share the project (poses) and
subset of the first board: feedback + region boards yes, the child board (consensus poses) and level-2 boards (subsets) no.
Verified 2026-10-17 00:36 on the toy (run_served_test.sbatch): served board B identical to B alone (cell corr 1.000, argmax 1.000).
Toy dev chain (200k, DEV=1, SERVER=1, 10 workers) end to end 68 min on 2026-10-17: parent 3, child 4, refine3D_states 27, feedback 4, three served region boards 4 min each but ~20 min of post-board work per region (level-2 host board + rec), landscape 2.
The sequential (host) board had the same single-threaded focus: the toy level-2 board 1193 s -> 113 s once routed through the shared `apply_focus_d` (2026-10-17 01:58). refine3D_states on the toy family (87k, 10 states, 10->8 A maxits 5) = 27 min: half master-bound state reconstructions on 8 threads, half 10-part alignment rounds (queue latency ~1 s); dev preset now NTHR_REF 32 master threads, NP 20 parts, family job -c 32.
LOCALDEV (2026-10-17 13:05, default with DEV=1 or SMOKE=1): no worker jobs; chain submitted with `-c 48 --mem=120G`, parent board on the chain job's host threads, each family job on a 48-core node (host boards NTHR=48, refine3D_states/reconstruct3D through `update_project qsys_name=local` on the project copies, NP 6 local parts, regions inline via f692_regions_local.sh). User: the dev loop must be minutes, 20 workers for 68k particles is absurd; kill long test jobs.
Smoke loop (SMOKE=1 LOCALDEV, 50k subset, 3-band schedules, 1 region/family) 2026-10-17: 23 min end to end on three 48-core nodes with zero worker jobs; host boards = stack reads 117-136 s per board (5 ms/particle, one process) + seconds per iteration (50k at 14 A 3.6 s, 10 A 8.6 s; 25k focused at 8 A 46 s). The load is the floor -> load each family once (local persistent workers) next. Traps hit: refine3D_states needs maxits >= 5 (smoke preset fixed), region script needs a PS fallback when the family had no refinement (fixed).
Smoke 3 (local persistent workers, 6 x 8 threads): 34 min end to end; parent 124 s, child 216 s, feedback 220 s (served), refinement 20 min for 25k particles / 4 states (master-side, still the bottleneck); one SIMPLE segfault in refine3D_states (3 states). Traps: parts inherit OMP_NUM_THREADS (set 8, not 48); never put a `#` comment mid-line in a one-line shell block (it swallowed the region script's definitions); `set -u` + an unset $SERVER killed the dev chains at launch (use ${SERVER:-1}).
Smoke 4 (15:37, 43 min end to end, 48-core family nodes): child 228 s, refinement 20 min for 25k particles / 4 states / 6 local parts (THE floor: master-side volume assembly through BeeGFS), feedback 230 s (served), region board 450 s, region found no supported states at 25k (smoke subset too small for regions). `nproc` inside a job reports OMP_NUM_THREADS, not the allocation (use SLURM_CPUS_ON_NODE). DEV-preset families had FAMC=32 (oversubscribed 6x8 workers) -> now 48 always in LOCALDEV.
Smoke 5 (local-disk refinement dir, 50k/family, 5 states, 6 local parts): refinement 27 min - NOT I/O: worker rounds total ~340 s, master gridding 5 x 16 s, the remaining ~20 min is inside the 8 volassembles (`ML-regularized map (solve_kind=gridding_regularized)` per state and half). The refinement floor is SIMPLE's regularised volume solve, proportional to states x iterations; local disk refuted.
Smoke 6 (ml_reg=no passed through REFINE_EXTRA): refinement 19.5 min vs 27 min with the regularised solve (same 50k / 5 states / 6 local parts) - a 30 % cut, still the floor. Trap: changing a board's nfree mid-run mixes script versions (a served board must have the same lattice and nfree as the server's first board, else `board server: every board must have the same lattice and free components`).
Smoke 7 (18:54, the consistent dev config: nfree 1 everywhere, ml_reg=no, local workers + server, 100k subset): 61 min end to end, families finished at 58 min; 43 states at 0.7, 6/10 conformers (cell count). The dev loop floor at 100k is ~1 h with the refinement; ~35 min without it.
Chain: f692_chain9 (DEV=1, SERVER=1) -> f692_family_job3 -> f692_regions_par3 -> f692_region_one3.

**Why:** the board stage dominated every chain (child 22 min, feedback 37 min, region boards 25 min each); after the fix the
refinement (65-93 min at 335k, maxits 6 from 20 A) is the bottleneck.
**How to apply:** read the timing line before blaming I/O or queue latency; keep master-side per-cell loops parallel.
