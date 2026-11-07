---
name: newest-step-dir-sort-trap
description: "TRAP: `ls -d $W/[0-9]*_prg | sort -n | tail -1` on FULL paths sorts lexically (leading '/' is not numeric), so 10_selection loses to 5_selection; sort basenames numerically or name the step dir explicitly"
metadata:
  type: feedback
---

2026-10-25 (RNase P full data): the consensus-3D job looked up "the newest selection step" with
`ls -d $W/[0-9]*_selection | sort -n | tail -1`. SIMPLE had made `10_selection` (9 was a manual experiment dir), but with full
paths `sort -n` compares from the leading '/', finds no number, and falls back to lexical order, so it returned `5_selection`
(the older round-2 selection). The particle-count check caught it (3,972,259 vs expected 3,927,919) and nothing ran.

**Why:** step numbers pass 9 quickly in a long chain, and a wrong project silently refines the wrong particle set.

**How to apply:** `cd $W && ls -d [0-9]*_prg | sort -n | tail -1` (basenames), or `sort -t/ -k<last> -n`, or use the exact
step dir the program printed (`>>> EXECUTION DIRECTORY:`). Keep a count/identity check after every selection step.
Related: [[persist-project-lineage-no-new-projects]].
