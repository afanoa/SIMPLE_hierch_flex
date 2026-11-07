---
name: start-sessions-without-app-worktree
description: "STANDING: start side sessions with use_worktree=false; app-made worktrees fail on this repo ('shadow-unavailable', 4 failures 2026-11-07); a code-writing session makes its own git worktree with `git worktree add`"
metadata:
  node_type: memory
  type: feedback
  originSessionId: cc23038c-fea9-4b5f-b321-d58eb2942350
  modified: 2026-11-07T23:11:49.789Z
---

2026-11-07: every side session I started with `use_worktree: true` after ~15:20 failed before receiving its task ("Background
full checkout failed: shadow-unavailable" — the app's own git setup; 3 sessions at 15:21, the fission-split session at 23:06,
which failed again when the user retried it, because a retry reuses the original settings). Every session started with
`use_worktree: false` ran. I had already switched at 16:08 and then repeated the mistake. User: "I have had 4 sessions fail,
yet you keep doing this", "you are an orchestrator. Please orchestrate correctly".
Context: the repo is large (3.2 GB .git, 87 worktrees, big untracked files in the main checkout) and the disk is ~92 % full.

**Why:** a failed start loses the task and the user's time; the user expects the orchestrator to learn from the first failure.
**How to apply:**
- Always `use_worktree: false` for side sessions in this project.
- A session that writes code: first step in its brief = `git -C /Users/atifao/SIMPLE worktree add /Users/atifao/SIMPLE_<task> <commit>`
  and work only there; never edit, stage, stash or switch branches in `/Users/atifao/SIMPLE` (the user's dirty checkout).
- After any session start failure, read the failure reason and change the setting before retrying; never retry unchanged.
