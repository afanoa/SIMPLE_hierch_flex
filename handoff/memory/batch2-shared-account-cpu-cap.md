---
name: batch2-shared-account-cpu-cap
description: "my batch2 jobs run under the user's own account and its ~642-CPU QoS cap, so a big experiment (8 masters x 10 workers x 8 threads) queues the user's own jobs behind mine; keep my footprint to roughly a third of the cap when they are running"
metadata:
  type: feedback
---

2026-09-22 ~00:40: while my paired pose test (8 refine3D masters, 80 worker jobs at 8 threads) plus three
flex resumes were running, the user's own `flex2D_bt16..19` jobs sat PENDING. Same account, same QoS
(csbdevel caps one user at ~642 CPUs, see [[batch2-csbdevel-node-classes]]), so my jobs were blocking
theirs. I cancelled the ribosome pose arms and their workers.

**How to apply:** before submitting a multi-arm experiment, `squeue -u atifao` and count CPUs already in
use by the user's own jobs (names not starting with `flexw_`/`simple_` from my scripts); keep mine to
about a third of the cap and stagger arms rather than fanning out. Never cancel a job I did not submit.
Related: [[flex-tests-run-on-cluster-only]], [[never-rebuild-under-running-jobs]].
