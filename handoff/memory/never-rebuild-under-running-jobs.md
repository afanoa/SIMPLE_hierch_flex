---
name: never-rebuild-under-running-jobs
description: "`make install` into a build tree whose binaries are running on the cluster kills those jobs with SIGBUS (exit 7) minutes later; use a separate build/prefix per experiment or wait for the jobs"
metadata:
  type: feedback
---

2026-09-22 on batch2: four refine3D pose-test arms (10180, 10 slurm workers each) lost most workers with
`Bus error (core dumped)` (sacct ExitCode 7:0, MaxRSS 6.5 GB, so not OOM) about 10 minutes after I ran
`make -j32 install` in `~/SIMPLE_flexw/build` to test a new consensus module. `install` overwrites
`bin/simple_private_exec` in place; a running process that page-faults on the modified executable gets
SIGBUS. The masters then waited forever for the missing `algndoc_NN.simple` files. Cost: the whole pose
test had to be resubmitted after the last build.

**Why:** the batch2 experiments run straight from the build directory (`SIMPLE_PATH=$HOME/SIMPLE_flexw/build`),
so the binary under test is the binary being rebuilt.
**How to apply:** before any `make install`, `squeue -u atifao` must show no job using that tree; otherwise
build into a second tree (`~/SIMPLE_<name>` worktree + its own build) or wait. Cheap flex resumes can share
a tree only when nothing else is running from it. Related: [[flex-tests-run-on-cluster-only]],
[[build-on-batch2-login-node]].
