---
name: never-push-before-build-and-gate
description: never push to origin before the commit has been built on batch2 and passed the fast gate, even for comment-only changes; the Mac syntax check is not a test
metadata:
  type: feedback
---

Order is always: commit -> push scratch ref to batch2 -> build -> fast gate (13/13) -> only then
push to origin. This holds for comment-only or "obviously trivial" changes too.

**Why:** 2026-10-01 I pushed the flex comment-port commit (c961fb432) to origin/master with only
the Mac syntax gate, and started the batch2 build afterwards. The user: "not sure why you would
push before testing... you need to learn better development practices."

**How to apply:** treat the batch2 fast gate as the push gate, per
[[compile-and-test-means-batch2-fast-gate]]; the syntax-only check is a pre-commit convenience,
never a release check. Related: [[small-fast-cluster-steps]], [[push-means-master-after-rebase]].
