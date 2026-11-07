---
name: plain-ssh-sbatch-lacks-gcc-runtime
description: sbatch from a plain (non-login) ssh command exports no GCC 15.2 runtime, so SIMPLE dies at launch with GLIBCXX_3.4.32 not found but tee makes Slurm say COMPLETED; every script must export the gcc lib path
metadata:
  type: feedback
---

2026-10-02: the `not` default rerun 53549694 was submitted via `ssh batch2 '... sbatch script'` (no `bash -lc`), so the
job inherited no module environment; `simple_exec` failed instantly (`/lib64/libstdc++.so.6: GLIBCXX_3.4.32 not found`,
`GFORTRAN_10 not found`) and, because the script pipes through `| tee`, sacct reported COMPLETED 0:0. It sat dead ~9.5 h
before I looked. Same trap hit the ~/SIMPLE_dbg build's UI validation step (nohup without the module env).
**How to apply:** every job script I write exports
`PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH` and `LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64:...`
itself (never rely on the submitting shell), and after every submit I confirm the log shows real SIMPLE output
(EXECUTION DIRECTORY / first stage line) within a minute, not just RUNNING in squeue. Related: [[flex-refine3dstates-gridding-and-knobs]]
(tee hides failures), [[box-build-needs-gcc16-ld-library-path]].
