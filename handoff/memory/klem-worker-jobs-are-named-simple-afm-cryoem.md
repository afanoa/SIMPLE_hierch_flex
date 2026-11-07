---
name: klem-worker-jobs-are-named-simple-afm-cryoem
description: TRAP: the distributed klem board workers (simple_private_exec prg=klem_board_part) and reconstruct3D workers are submitted by SIMPLE's qsys with the PROJECT compenv's job name "simple_afm_cryoem" (+ --mail-user=my.name@uni.edu --mail-type=FAIL from the same compenv); they are klem jobs, not an afm project; count them as klem cores and never cancel/ignore by that name
metadata:
  type: project
---

2026-10-16 22:30: 480 of the user's 640 running cores were "simple_afm_cryoem" jobs = the klem_worker_NN scripts of
f692c/f692d (WorkDir rnasep_klem/f692c|f692d, Command .../klem_worker_06, prg=klem_board_part) and the
1_reconstruct3D worker jobs. I misreported them as another project. The name and the mail lines come from the
692k project's compenv (rnasep_coupled/full692k_2d/rnasep.simple, copied into every stage project); fix for the next
run = update_project on the source project (job_name, user_email) inside a short sbatch (not on the login node).
Related: [[scancel-only-by-job-id-never-by-name-pattern]], [[no-slurm-email]], [[klem-f692-fortran-chain]].
