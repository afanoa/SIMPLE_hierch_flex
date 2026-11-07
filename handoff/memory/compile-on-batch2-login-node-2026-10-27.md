---
name: compile-on-batch2-login-node-2026-10-27
description: "user 2026-10-27: compile directly on the batch2 login node (ssh batch2.ncifcrf.gov = fsitgl-xfer03p, full hostname); SLURM build jobs only when batch2 is down; head01p never builds"
metadata:
  type: feedback
---

User, 2026-10-27, relayed verbatim by the coordinating session: "you can compile on batch2. you do not need to submit a slurm job
unless that node is down".

**Why:** build jobs wait in a full queue; the login node builds at once. The 2026-10-24 SLURM-build rule was only for when batch2 was down.

**How to apply:**
- Build with `ssh batch2.ncifcrf.gov` (the short alias `batch2` fails host-key verification from the Mac).
- Export PATH and LD_LIBRARY_PATH for /mnt/nasapps/production/gcc/15.2.0 and cmake/4.0.3 inside the build command.
- Fall back to a SLURM build job only when batch2 is down. batch.ncifcrf.gov (head01p) never builds.
- No SIMPLE runs, scorers or tests on any login node (sbatch only).
- Never rebuild an install tree that running jobs use ([[never-rebuild-under-running-jobs]]).

Supersedes the SLURM-only part of [[build-on-batch2-login-node]].
