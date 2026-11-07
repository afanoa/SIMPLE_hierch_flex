---
name: push-means-master-after-rebase
description: when the user has had a branch rebased onto origin/master and says "push it", they mean fast-forward origin/master, not publish a feature branch
metadata:
  type: feedback
---

"Push it" after a rebase-onto-master means `git push origin <branch>:master` (fast-forward), not
`git push -u origin <branch>`. Before pushing, show the user exactly which commits will land
(`git log --oneline origin/master..HEAD`) and the shortstat, then push on their go.

**Why:** 2026-10-01 I pushed flex-pca-refactor as a remote branch; the user replied "delete that
branch. i wanted you to push onto master", then asked to verify only one commit would appear.

**How to apply:** the SIMPLE workflow commits to master directly (no PRs); squash + rebase +
fast-forward is the pattern. Related: [[no-coauthor-trailer-on-commits]], [[small-fast-cluster-steps]].
