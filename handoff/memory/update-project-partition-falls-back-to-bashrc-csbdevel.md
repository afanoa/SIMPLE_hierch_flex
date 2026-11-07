---
name: update-project-partition-falls-back-to-bashrc-csbdevel
description: "TRAP: update_project without qsys_partition writes $SIMPLE_QSYS_PARTITION (csbdevel, from ~/.bashrc on batch) into compenv, overriding a project's norm; pass qsys_partition on every update_project"
metadata:
  node_type: memory
  type: reference
  originSessionId: 371593ba-f05b-45f2-801d-d68e405eef0b
  modified: 2026-11-05T22:11:15.168Z
---

`sp_project%update_compenv` (src/main/project/simple_sp_project_core.f90) sets `qsys_partition` from the command line, or else from the environment variable `SIMPLE_QSYS_PARTITION`. It never keeps the value the project already holds. `~/.bashrc` on batch exports `SIMPLE_QSYS_PARTITION=csbdevel`, and sbatch passes the submitting environment through.

So an `update_project` that only changes `job_memory_per_task` silently moves every later distributed worker to csbdevel.
- 2026-11-05: the IgM disc 2D (dj_ws2 job 54961009) was copied from a norm project, but its 10 x 32-cpu workers went to csbdevel. The coordinating session caught it.
- A copy of a project that went through such an update carries csbdevel into the next step unless that step's update_project sets the partition again.

**How to apply:** give every `update_project` the whole queue set: `qsys_name=slurm qsys_partition=<part> job_memory_per_task=<MB> nthr=<n> nparts=<n>`. Check a generated worker script's `#SBATCH --partition` (or `print_project_info`) once the run fans out. IgM work stays on norm.

Related: [[simple-script-yes-compenv-trap]], [[norm-partition-when-csbdevel-capped]].
