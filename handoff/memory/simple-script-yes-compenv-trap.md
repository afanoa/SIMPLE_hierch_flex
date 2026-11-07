---
name: simple-script-yes-compenv-trap
description: "script=yes writes #SBATCH directives ONLY from the project's compenv segment; params and SIMPLE_QSYS* env vars are ignored when a projfile is present, so fix the queue with update_project first"
metadata:
  type: reference
---

`script=yes` (simple_exec_helpers.f90 `script_exec`) calls `qsys_env%new`, which does:

```
if( params%projfile /= '' ) then     ! <-- true for every real run
    self%qdescr = compenv_o%ori2chash()   ! queue settings come from the PROJECT
else
    call init_qdescr_from_runtime(...)    ! only here are params/env consulted
```

So with a projfile present **`qsys_name=slurm` on the command line and `SIMPLE_QSYS` /
`SIMPLE_QSYS_PARTITION` in the environment are all ignored**. If the project's compenv says
`qsys_name=local`, `generate_script_3` takes the local branch and emits a bare shell script with
no `#SBATCH` lines at all, which looks like script=yes is broken.

**Recipe on a SLURM cluster:**
1. `simple_exec prg=update_project projfile=proj.simple qsys_name=slurm qsys_partition=<part> job_memory_per_task=<MB> nthr=<n> nparts=<n>`
   (verify with `prg=print_project_info`; `job_memory_per_task` is MB, default 16000)
2. `simple_exec prg=flex_pca projfile=proj.simple script=yes <args>` -> writes `<prg>_script`
3. `sbatch <prg>_script`

**Three gotchas after generation:**
- `update_project` sets `user_email=nil`, which emits `#SBATCH --mail-user=nil`; drop that line.
- The header writes `./stderrout/outfile.%j`, so `mkdir -p stderrout` first or the job cannot start.
- `script_exec` **deletes the `mkdir` key** before writing the script, so the run defaults to
  `mkdir=yes` and works in `1_<prg>/`, not the directory you launched from. Passing `mkdir=no` at
  generation does NOT survive into the generated script.
- Env vars the run needs (e.g. `SIMPLE_COV_*`, `SIMPLE_PATH`) must be injected into the generated
  script or exported in the submitting shell; sbatch exports the submitting environment by default.

Master thread boost = nparts x nthr, but the master's own allocation is only `cpus-per-task`, so
10x16 oversubscribes a 16-core master 10-fold on master-only stages (deconvolution, UMAP);
`SIMPLE_COV_MASTER_NTHR` caps it. Related: [[flex-pca-pcg-branch-state]].
