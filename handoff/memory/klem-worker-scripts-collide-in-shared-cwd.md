---
name: klem-worker-scripts-collide-in-shared-cwd
description: TRAP: the klem board master writes klem_worker_NN scripts into its cwd; concurrent masters with the same cwd (parallel region jobs) overwrite each other's scripts -> "Fortran runtime error: End of file" / hung boards; fboard.sh now cd's into the board dir
metadata:
  type: project
---

2026-10-17 night: 6 of 24 region boards of the corrected 692k runs (f692c/f692d) died with `Fortran runtime error: End of file`
on the master, and the family jobs then waited 12 h for DONE markers. Cause: the SIMPLE queue layer submits `<cwd>/klem_worker_NN`;
the parallel region jobs inherit the family job's cwd (the fam dir), so up to 4 masters overwrote each other's worker scripts
(two SLURM workers serving the same part of one board = garbled reply files = EOF; the other board never gets its worker).
Fix: `fboard.sh` runs the master inside `( cd $OUT; ... )` so scripts and logs live in the board dir. Failed regions were
resubmitted by hand (f692_region_one2 with the family env, NPARTS from a sibling board's log).

**Why:** any two boards launched from the same directory at the same time will hit this.
**How to apply:** never launch two board masters from one cwd; on an EOF in a board log check `distr/` for duplicated
part files and the worker script names; the region wait loop should also stop early when a region job has exited without DONE.
