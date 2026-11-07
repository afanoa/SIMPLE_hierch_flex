---
name: rebase-check-upstream-edits-landed
description: after rebasing a restructuring branch over upstream, verify line-by-line that every upstream edit to moved/split files landed; patch --fuzz silently dropped ~40% of a comment rewrite
metadata:
  type: feedback
---

After a rebase that resolves conflicts on files the branch moved or split, run a line-level
check of every upstream commit in the range against the new tree: each line upstream added must
exist somewhere, each line it removed must be gone (scratchpad `comment_gap` / `portcomments.py`
pattern: anchor each upstream change group on its exact removed lines wherever they now live).

**Why:** 2026-10-01 the flex refactor rebase re-applied upstream's "comment universe refactoring"
(adaaea4c7) with `patch --fuzz=3` and hand-handled rejects; 143/335 added and 291/1104 removed
comment lines were missed, and it was pushed to master before the user asked "did that transfer
correctly?". Fixed in a follow-up commit.

**How to apply:** before declaring a rebase done, list `git log base..upstream -- <moved paths>`
and run the line check per commit; report the counts to the user with the rebase result.
Related: [[small-fast-cluster-steps]], [[flex-pca-refactor-branch-state]].
