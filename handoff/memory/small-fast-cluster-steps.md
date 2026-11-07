---
name: small-fast-cluster-steps
description: the user rejects long chained ssh/batch2 pipelines and slow multi-step git operations; do each step as one short announced command and report between steps
metadata:
  type: feedback
---

Keep cluster and git work in small, fast, announced steps: one ssh call per step (push, checkout+cmake,
make, sbatch gate, sbatch arms, compare), each reported as it finishes. Never a single 30-minute
pipeline that pushes, builds, gates, runs six arms and compares in one Bash call.

**Why:** on 2026-10-01 a rebase took 39 minutes and the follow-up one-shot batch2 verification
command was rejected: "not sure why you are taking so long... stop wasting my time and sort this out".

**How to apply:** check `git fetch` + `merge-base --is-ancestor` first (origin/master moves between
turns), rebase, gate, push -f the scratch refs, then build in the background and submit jobs as soon
as the build finishes. Related: [[keep-test-loops-minimal]], [[announce-commands-before-running]],
[[head-node-runs-nothing-but-sbatch]].
