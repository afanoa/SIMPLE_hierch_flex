---
name: build-on-batch2-login-node
description: "batch2: compile ON batch2.ncifcrf.gov login node (user 2026-10-27; slurm build only if batch2 down), but NEVER run, test or score anything on the login/head node - every run, test, scorer goes through sbatch/srun (user, 2026-09-22: 'DO NOT submit on head node. use slurm' ... 'you can build on slurm just dont test there')"
metadata:
  type: feedback
---

2026-09-22: the user objected twice to work on the batch2 head node. Clarified: builds are fine (login node or a
slurm job), tests are not: nothing that runs SIMPLE (or a python scorer over maps) may run on the login node.

**How to apply:** `make install` on the login node (or `srun -p csbdevel -c 16 --mem=32G -t 60 make ...`) is OK;
everything else - simple_exec/simple_private_exec calls, scorers, quick "does it start" checks - via `srun -p csbdevel`
or `sbatch`. Compute nodes are Skylake-or-newer; a compute-node ARCHOPT build may SIGILL on the login node
([[batch2-archopt-sigill-login-node]]), another reason nothing runs there.
Related: [[flex-tests-run-on-cluster-only]], [[never-rebuild-under-running-jobs]].

UPDATE 2026-10-24 (user, after the batch2 login node fsitgl-xfer03p went down about 00:00): ssh to **batch.ncifcrf.gov** (head node
fsitgl-head01p; the same /mnt/beegfs and SLURM queue, and running jobs are unaffected). The user: "idk how to compile code. maybe you
can use slurm to submit a compile job because you cant compile on head node". So on the head node, do NOT compile: run builds as a
SLURM job (e.g. `sbatch -p csbdevel -c 16 --mem=32G -t 1:00:00` wrapping `make install`), and keep the head node to sbatch / squeue
/ sacct and file reads.

UPDATE 2026-10-27 (user): "you can compile on batch2. you do not need to submit a slurm job unless that node is down".
batch2 is back (batch2.ncifcrf.gov = fsitgl-xfer03p; the short alias `batch2` fails host-key checks from the Mac, use the full
name). Compile there directly (export the gcc 15.2.0 / cmake 4.0.3 PATH+LD_LIBRARY_PATH); SLURM build jobs only when batch2 is
down. batch.ncifcrf.gov (head01p) still never compiles; no SIMPLE runs/scorers/tests on either login node.
TRAP (2026-10-27): the batch2 login node is a Xeon Gold 6150 (Skylake); builds use -march=native, so a login-node build
targets Skylake while x6342 compute nodes are Ice Lake. It runs, but its codegen differs from x6342-built baselines (e.g. the
clean d3513ff70 install). For A/Bs (timing, near-tie decisions) build every arm on the same host type; a SLURM build on x6342
(amp_mstep_v1/build_variant.sbatch, ~4 min at 32 CPUs) matches the existing baselines.
