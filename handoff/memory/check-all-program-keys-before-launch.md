---
name: check-all-program-keys-before-launch
description: "Before submitting any SIMPLE chain, check EVERY key of EVERY program call against that build's UI (add_input list + globally required nthr/nparts/mskdiam/pgrp); 2026-10-25 the RNase P consensus 3D failed 3x one argument at a time"
metadata:
  type: feedback
---

2026-10-25, RNase P consensus 3D: the chain stopped three times in a row at launch: a wrong step dir, then `scale` without `nthr`,
then `refine3D_auto` rejecting `sigma_est` (gone in release 4). The coordinator: check every key of the refine3D_auto line first,
"so we don't relaunch one error at a time". The user was waiting. Earlier the pilot also failed twice the same way (motion_correct
needs nparts, stackops needs nthr).

**Why:** each relaunch costs the user 10-20 min of waiting and a round of messages.

**How to apply:** for each program in a script, before sbatch:
- `awk '/subroutine new_<prg>\(/,/end subroutine new_<prg>/' src/main/ui/simple/*.f90 | grep add_input` in the build that runs;
- the globally required keys `nthr`, `nparts` (when the program lists them), `mskdiam`, `pgrp` must be passed;
- any key not in the program's list must be a known global parameter; never assume an old key (e.g. `sigma_est`) still exists.
Related: [[keep-test-loops-minimal]], [[newest-step-dir-sort-trap]].
- 2026-10-28 TRAP: recipe13/14 gates launched WITHOUT the split env (SIMPLE_SPLIT_CONTRAST=maps SIMPLE_SPLIT_RHOMIN=0.03; the rho
  floor only applies when the env is set) -> bare split test, qtoynull root split K*=2 (z 14.4). Env switches are keys too: every
  launcher writes its env into the run folder and refuses to start without the production split env.
