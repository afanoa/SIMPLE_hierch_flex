---
name: autonomous-build-and-submit-2026-10-27
description: "STANDING (user 2026-10-27): work autonomously across this chat and spawned sessions; no permission needed to build, deploy or submit jobs; push/delete rules unchanged"
metadata:
  type: feedback
---

User, 2026-10-27, after asking for harder testing on the 1242-core budget: "Please work autonomously over this chat and
others. You do not need to ask for permissions to build, deploy or submit jobs".

**Why:** the user wants the test campaign (amplitude A/B, split/merge rules, start-point sensitivity, porting/optimizations)
to keep the cluster busy without waiting on approvals; idle cores are the complaint, not overuse.
**How to apply:** build (SLURM job, never the head node), install builds on the cluster and sbatch tests without asking;
keep sessions busy and report results/core-h after. NOT covered: pushing to origin (never), pushing to the fork or
committing (still on the user's ask), deleting cluster data (never), rebuilding under running jobs (never). Keep the
per-partition running+pending check (csbdevel <= 642, norm <= 600). Related: [[batch2-400-core-limit-and-worker-timeouts]],
[[never-push-to-origin-without-permission]], [[never-delete-only-read-or-move]].
