---
name: never-overwrite-scripts-under-running-jobs
description: TRAP (2026-10-15) - scp'ing a patched .sbatch/.sh over one a running job is still executing made bash fail with "syntax error near unexpected token done" mid-job (bash reads scripts incrementally); patch under a new name or wait
metadata:
  type: feedback
---

Overwriting `run_expand_toy.sbatch` on batch2 while the annealing job was still executing it (the spectral control
pass) broke the job at the loop's `done`: bash reads a script file incrementally, so the running shell saw the new
file's bytes at its old offset. The children had finished; only the merge step was lost and had to be re-run.

**Why:** batch2 jobs `bash $W/run_*.sbatch` the shared scripts directly; any live edit races with every job that is
mid-way through that file (also the nested calls: run_expand_toy -> run_klem -> klem.py; python reads the whole file at
start, bash does not).

**How to apply:** before `scp` of a shell script into routes/toy55 (or any job dir), check `squeue` for jobs that call
it (directly or via a parent script) and either wait, or copy the patched version under a NEW name and submit that.
Python files are safe to overwrite for running processes but not for jobs that will start a later step with them
(then the new version runs, which is usually what you want). Related: [[klem-kernel-latent-em-2026-10-15]],
[[chain-steps-inside-one-job]].

2026-10-25: I clobbered stree_port/stree_noderep.sbatch (an Oct 22 node REPORT script) with a new replay script of the same name. My
`ls X; scp ... X` printed the existing path but the scp ran anyway. I restored it from the session scratchpad copy (Oct 22 14:04, md5
78320f76) and moved the new tool to stree_recb0rep.sbatch. Use a guard that STOPS: `ssh host "[ -e X ] && echo EXISTS"` first, then
copy only on its absence, or `scp` to a fresh, unique name.

