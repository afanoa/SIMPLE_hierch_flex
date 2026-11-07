---
name: particle-sieving-job-model
description: particle_sieving is the offline sieve program (not "sieve"); its nparts = concurrent chunk SLURM jobs, and script=yes strips mkdir so you must re-add mkdir=no
metadata:
  type: reference
---

The offline sieving program is `prg=particle_sieving` (simple_exec). Searching for
"sieve" misses it — `sieve_cavgs` is the *streaming* program (simple_stream, needs
`dir_target` from a running pick_extract) and is not what you want for an existing
project of extracted particles.

Job model (verified 2026-09-18 on batch2):
- The master job is a **lightweight orchestrator** — it only needs `nthr` cores.
- `nparts` on the CLI means "number of chunk abinitio2D jobs run concurrently".
  The commander remaps `nparts -> nchunks` and forces per-chunk `nparts=1`.
- Each chunk is submitted as its **own SLURM job** via
  `qenv%exec_simple_prg_in_queue_async`, getting `--cpus-per-task=nthr` and
  `--mem=job_memory_per_task` from the project compenv. So `nparts=10 nthr=8`
  is 10 separate 8-core jobs, NOT 80 cores on one node.
- Preconditions: `os_mic`, `os_stk` and `os_ptcl2D` all non-empty and
  `nstks == nmics` (otherwise use prune_project).
- Chunking boundary: `nmics` (default 50) micrographs OR `nptcls_coarse`
  (default 5000) particles, whichever hits first.

Two traps:
1. `script_exec` **deletes `mkdir` from the cline** before writing the script, so
   the commander's own `mkdir=yes` default applies at runtime and the sieved
   project lands in `<n>_particle_sieving/`. To get a deterministic path for a
   chained job, `sed` `mkdir=no` back into the generated script's command line.
2. The sieved project is safe to chain abinitio2D off: chunks downscale with
   `box_crop` only, so the merged project still points at the **original
   full-box stacks**. Sieving rebuilds the project, though, so re-run
   `update_project` before the chained job to guarantee compenv survives.

See [[simple-script-yes-slurm-submission]], [[simple-script-yes-compenv-trap]],
[[batch2-csbdevel-node-classes]].
