---
name: no-slurm-email
description: never put --mail-type/--mail-user in batch2 job scripts; use --mail-type=NONE (user, 2026-10-02)
metadata:
  type: feedback
---

2026-10-02 the user: "stop emailing me via slurm". Every sbatch script I write gets `#SBATCH --mail-type=NONE` and no
`--mail-user`; when taking over or resubmitting a job that has mail on, clear it with
`scontrol update JobId=<id> MailType=NONE`. Supersedes the "SLURM mail must go to afan.atif@nih.gov" line in
[[trpm4-flipqr-sieve-campaign]]. SIMPLE's own worker scripts carry `--mail-type=FAIL --mail-user=nil` from the
project compenv (user_email=nil), which does not reach the user.
**Why:** job-end/failure mails from my many runs were noise. **How to apply:** check new scripts for mail lines before sbatch.
