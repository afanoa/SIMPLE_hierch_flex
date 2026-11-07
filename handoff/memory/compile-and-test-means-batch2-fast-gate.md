---
name: compile-and-test-means-batch2-fast-gate
description: "compile it and run tests" = build on batch2 (not the Mac) + the ~30 s fast gate via sbatch; no multi-hour library suites unless asked
metadata:
  type: feedback
---

2026-09-28 (flex-GPU removal on master): I started a local Mac build and planned the fast gate plus
lib_heterogeneity (2 h cap). The user redirected twice: "just compile it on batch2..." and "dont need
a 2 hour test...".

**Why:** batch2 is where the user builds and tests; the Mac toolchain is also currently broken (the
Sep-22 CLT update points MacOSX.sdk at SDK 27, whose libSystem.tbd the conda ld64-956 cannot read;
SDKROOT/CMAKE_OSX_SYSROOT=MacOSX26.5.sdk works around it). Long nightly-tier suites are not wanted
for a routine change.

**How to apply:** for "compile and test", ship the diff to a fresh batch2 worktree (fetch origin there,
`git worktree add --detach`, `git apply --index` over ssh), build on the login node, then one short
sbatch on csbdevel running `scripts/run_fast_gate.sh build` (-t 20 min). Offer library/highlevel
suites only as an optional follow-up. Related: [[build-on-batch2-login-node]],
[[head-node-runs-nothing-but-sbatch]], [[announce-commands-before-running]].
