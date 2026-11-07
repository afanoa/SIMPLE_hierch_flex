---
name: simple-chdir-failure-does-not-abort
description: "TRAP: a plain file named N_prg does NOT stop a mkdir=yes SIMPLE step; simple_chdir only warns and the program runs in the project root (ribosome_10509 extract, 2026-10-20)"
metadata:
  type: feedback
---

2026-10-20, ribosome_10509: to hold a chain (pick -> extract) for a box-size decision I placed a FILE named `10_extract`
so extract's mkdir=yes would fail. simple_mkdir failed, but simple_chdir only printed "failed to change path" and did NOT
abort. Extract ran in the project root: it overwrote the root project file and wrote ~3,600 particle stacks, 20 worker
scripts and stderrout/ there, and launched 20 slurm workers before I cancelled it 1.5 min later.

**Why:** the user wants clean, readable project dirs ([[persist-project-lineage-no-new-projects]]); a guard that relies
on SIMPLE failing is not a guard.

**How to apply:** never chain a step that waits on a user decision. End the job script before it, or have the script
test for a sentinel between steps (`[ -e $W/scripts/HOLD ] && exit 0`). To stop a running chain, cancel it by job ID
([[scancel-only-by-job-id-never-by-name-pattern]]) between steps. Never count on a SIMPLE program erroring out.
