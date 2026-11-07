---
name: scancel-only-by-job-id-never-by-name-pattern
description: 2026-10-16 I cancelled "rnasep_f692_fam*/reg*" by squeue name pattern to kill the f692b replicate and killed the TOY chains' family/region jobs too (the family-job scripts reuse the job name across runs); cancel only explicit job ids traced to the run's own log lines ("family 1 -> job N")
metadata:
  type: feedback
---

**Why:** f692_family_job.sbatch and f692_region_one.sbatch are submitted with -J rnasep_f692_fam$j / rnasep_f692_reg*
by every chain (real, replicate, toys), so a name pattern spans runs; the toys lost ~1 h of region work and had to be
resubmitted with their environment recovered from /proc/<pid>/environ via srun --overlap.

**How to apply:** build the cancel list from the chain's own log ("family j -> job N", "region x -> job N") or from
squeue filtered on the WorkDir/run dir, never on the job name; show the list before scancel.
