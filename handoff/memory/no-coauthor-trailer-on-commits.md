---
name: no-coauthor-trailer-on-commits
description: never add a Co-Authored-By trailer to git commits in this repo — the user does not want Claude listed as a co-author
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 34cf3bd2-fff1-4dd4-bda2-2ec6d52c01a4
  modified: 2026-08-03T17:49:35.754Z
---

Do not append `Co-Authored-By: Claude ... <noreply@anthropic.com>` to commit messages in
this repository. Commit messages should contain only what the user asks for.

**Why:** stated 2026-08-03. The user saw the trailer on a pushed commit, asked for it to be
removed, and confirmed the preference for future commits. GitHub renders that trailer as a
second author avatar on the commit, which they do not want on their work. Note the trailer is
the *only* thing that produces this — the `Author` and `Committer` fields were already the user
alone.

**How to apply:** this overrides the default harness guidance that says to end commit messages
with the trailer. Write the message the user asked for and stop there. Removing it after a push
costs a force-push and a rewritten SHA on a shared master, so get it right the first time.
See [[flex-pca-rename-and-cleanup]].
