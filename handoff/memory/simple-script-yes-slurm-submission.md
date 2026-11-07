---
name: simple-script-yes-slurm-submission
description: "How to run SIMPLE distributed programs on a slurm cluster: script=yes makes simple_exec WRITE a <prg>_script sbatch file and exit (it does NOT run), you then sbatch it. Two headers in that generated script are always wrong and must be patched: --time=0-0:0:0 (would be killed instantly) and --mem, which ignores job_memory_per_task entirely."
metadata:
  type: reference
---

**The mechanism (verified 2026-09-13):** `script_exec` runs at the TOP of `simple_exec` for EVERY program
(production/simple_exec.f90, implemented in src/utils/simple_exec_helpers.f90). If `script=yes` is on the command line it:
requires `projfile`, forces `mkdir=no`, strips `script`/`mkdir`, builds a qsys_env with nparts=1, calls
`qenv%gen_script(...)` to write `<prg>_script` (output tee'd to `<PRG>_OUTPUT`), then `call exit` -- the program never runs.
You then `sbatch <prg>_script`. This applies to shared-memory programs too, NOT just distributed ones.

Do NOT conclude from grepping `params%script` that script=yes is inert -- this path reads `cline%get_carg('script')`
directly, so it does not appear in a `%script` grep. I got that wrong once.

**Why you must sbatch rather than run the master directly:** a master left on the login node gets killed. sbatch the
generated script; the master then submits its own part-jobs from the compute node (nested sbatch works on batch2).

**Two headers in the generated script are always wrong:**
1. `#SBATCH --time=0-0:0:0`. job_time is `time_per_image * partsz`, but script_exec builds the qenv with no image
   count so partsz=0. Patch it (e.g. `--time=1-00:00:00`). The PART jobs the master spawns later compute their own
   sane walltime from the real partsz, so only the master script needs this.
2. `#SBATCH --mem=` ignores `job_memory_per_task` on the command line. `init_qdescr_from_runtime`
   (src/utils/qsys/simple_qsys_env.f90) hardcodes `JOB_MEMORY_PER_TASK_DEFAULT` = 16000, and there is NO env override
   (only SIMPLE_PATH / SIMPLE_QSYS / SIMPLE_QSYS_PARTITION / SIMPLE_EMAIL are read). Passing it on the cline only
   reaches the inner command, not the header. Editing the generated script is the only lever for the MASTER; the
   part-jobs are generated at runtime and always get 16000 MB -- fine for ctf_estimate, watch it at extract/3D.

**How to apply:**
```
cd <project root>                      # numbered dirs are FLAT under the root, plus one stderrout/
mkdir -p stderrout                     # --output=./stderrout/... ; job fails without it
simple_exec prg=<prg> projfile=<abs path to .simple> nparts=10 nthr=8 <args> script=yes
sed -i 's|^#SBATCH --time=.*|#SBATCH --time=1-00:00:00|; s|^#SBATCH --mem=.*|#SBATCH --mem=32000|' <prg>_script
sbatch <prg>_script
```
`qsys_name` is the key (NOT `qsys`); on batch2 `SIMPLE_QSYS`/`SIMPLE_QSYS_PARTITION` are already set in the shell
profile, so slurm+csbdevel land in the project compenv automatically. See [[empiar-10330-cs-corrected-optics]],
[[announce-commands-before-running]], [[simple-run-parallelism-preference]].
