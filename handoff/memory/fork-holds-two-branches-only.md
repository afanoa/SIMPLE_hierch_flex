---
name: fork-holds-two-branches-only
description: "STANDING (user 2026-11-07): the fork holds only hetero-tree (known to work) and ONE wip branch; no archive, per-session, experiment or handoff branches on GitHub"
metadata:
  node_type: memory
  type: feedback
  originSessionId: cc23038c-fea9-4b5f-b321-d58eb2942350
  modified: 2026-11-07T06:13:41.847Z
---

On 2026-11-07 I had pushed 66 `mac/*` archive refs, five `wip/*` branches, a 7,300-file handoff branch and eight doc commits onto
hetero-tree. The user: "not sure what happened to the hetero_tree branch. there is too much nonsense on my fork now", then "Most of
these branches are useless. Really just wanted a hetero tree we knew works. And then a wip. Most other things can go".

Cleanup the same night: hetero-tree reset to the gated `e1d5fdd4c`; one `hetero-tree-wip` (ports + testing guide + `handoff/`
memory and notes); everything else bundled to batch2 (`mac_archive_20261107/fork_archive_20261107.bundle`) and deleted from the fork.

**Why:** the fork is the user's public, working view of the code: one branch that works, one in progress.
**How to apply:**
- New work goes onto `hetero-tree-wip` as commits; promote to `hetero-tree` only when it is gated and the user says so.
- Archives and backups go to batch2 (git bundle or tarball), never as refs on the fork.
- Keep doc updates to one commit per change set, not a commit per edit.
- `master` and `klem-consolidated` on the fork predate this rule; leave them alone unless the user asks.

Related: [[never-push-to-origin-without-permission]], [[mac-retired-handoff-2026-11-07]], [[hetero-tree-branch-published-on-fork]].
