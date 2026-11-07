---
name: klem-board-convergence-and-exchange-2026-10-16
description: klem board timing facts (2026-10-16): per-iteration cost 11-43 s at 20-10 A, ~100 s at 8 A (exchange floor, not particle-bound: 30k and 307k boards cost the same -> small boards on the host); single-precision exchange saves 5-10 %; early stop by argmax change (18-38 %/it) or model change (11-18 %/it at the last iteration of every band) never fires: boards are under-iterated, not over-iterated
metadata:
  type: project
---

Measured on the 692k parent (10 workers x 8 threads, lps 20 20 20 14 14 14 10 10 10): iteration times 41/11/11 s
(20 A), 24/20/19 s (14 A), 43/39/38 s (10 A); whole board 5 min. Region boards at 12..8 A: ~100 s/iteration for
30k AND 307k particles (file exchange of 49 cells x 2 halves x 59^3 cubes dominates) -> level-2 boards run on the
host threads (NPARTS unset -> nparts=1; fboard.sh now passes nparts=${NPARTS:-1} because the UI marks it required).
Single precision (pnum4/pden4) in the exchange: identical readouts, 274 vs 325 s.
Convergence metrics now printed every iteration: "argmax changed" (fraction of particles whose argmax cell moved:
0.18-0.38 at every iteration, the soft lattice keeps neighbours exchanging) and "model change" (||dV|| / ||V - cons||
over the grid cells: 0.32 -> 0.17 at 20 A, 0.18 -> 0.18 at 14 A, 0.16 -> 0.11 at 10 A). Early stop (stop_frac, default
0.03 on the model change) never fires with 3 iterations per band; if anything the fine band wants more iterations.
Refine3D_states, not the boards, is the wall-time floor (63 min / family). Trap repeated: cp -f over fboard_ps.sh
while a job executed it -> "rc: unbound variable" (use scp to .new + mv). Related: [[klem-f692-fortran-chain]],
[[rnasep-klem-boards-need-abs-ctf]], [[never-overwrite-scripts-under-running-jobs]].
