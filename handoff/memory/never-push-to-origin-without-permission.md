---
name: never-push-to-origin-without-permission
description: "NEVER push anything to origin (github hael/SIMPLE) without the user's explicit instruction for that push; pushing to the user's private batch2 build remote is fine. 2026-10-01: I pushed flex-cls-split to origin unasked and the user objected"
metadata:
  node_type: memory
  type: feedback
  originSessionId: 03e8b343-4213-4208-890c-56346e7de4b2
  modified: 2026-10-24T22:24:37.497Z
---

2026-10-01: after merging flex-cls-split onto master I pushed the branch to origin because the user
had said "work on master". The user: "who said you could push, delete that branch?" and
"you DO NOT push without my permission". The delete itself was then blocked by the permission
classifier, so the user had to remove it.

**Why:** origin is the shared repository; a push is outward-facing and publishes work the user
has not reviewed. "Work on master" means base the work on master, not publish it.

**How to apply:** never `git push origin ...` (branches, tags, deletes) unless the user asks for
that specific push in that turn. Pushing to the user's own batch2 build remote
(`git push batch2 <branch>:refs/heads/<branch>-in`) is part of the build loop and is fine.
When work is ready for origin, say so and give the one-line push command instead of running it.

**STANDING (user, 2026-10-24): "we push only to the fork".** Publishing pushes go ONLY to the `fork` remote
(git@github.com:afanoa/SIMPLE_hierch_flex.git), never to origin. The state-tree work lives on fork branch `hetero-tree`.
- 2026-10-24, on the user's ask: pushed d3513ff70 there as a fast-forward, no force. Commits at or below it must not be rewritten.
- Still only on an explicit ask, and only after the batch2 build + fast gate ([[never-push-before-build-and-gate]]).
- Show which commits land (`git log --oneline fork/<branch>..<commit>`) and their shortstat, no trailers ([[no-coauthor-trailer-on-commits]]).
Related: [[no-coauthor-trailer-on-commits]], [[act-on-the-literal-ask]].
