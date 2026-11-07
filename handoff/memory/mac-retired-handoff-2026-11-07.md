---
name: mac-retired-handoff-2026-11-07
description: "2026-11-07 the dev Mac was retired: start from the testing guide on fork/hetero-tree-wip (with handoff/ memory and notes); every other ref is in a git bundle on batch2, every untracked file in a batch2 tarball; the fork is PUBLIC (manuscripts never pushed)"
metadata:
  node_type: memory
  type: project
  originSessionId: cc23038c-fea9-4b5f-b321-d58eb2942350
  modified: 2026-11-07T06:14:10.749Z
---

User 2026-11-07: "I am not going to have access to this computer anymore (I think) so I would like to be able to pick up from
wherever we leave off". After the cleanup the same night ([[fork-holds-two-branches-only]]):
- fork/hetero-tree = `e1d5fdd4c`: the gated master merge (`ee26d9bbf`) + verdict notes, the version known to work.
- fork/hetero-tree-wip = hetero-tree + the ports (`d142ea8af` depth_wiener, gated; `043b6bcfc`..`89604e8ac` site_split, gating) +
  the guide `doc/implementation_notes/hetero_tree_testing_guide_2026_11_07.md` (builds, gate fixtures and pass rules, harnesses,
  verdicts, open threads, decisions, how to restore this memory) + `handoff/` (`memory/`, `side_session_notes/`, `batch2_scripts/`,
  `mac_untracked_list.txt`, `ARCHIVE.md`).
- batch2 `/mnt/beegfs/elmlund/afan/mac_archive_20261107/`:
  - `fork_archive_20261107.bundle` (114 MB, incremental on origin/master + hetero-tree): 36 `mac/<branch>` (Mac branches with
    commits not on the fork), 30 `mac/dirty/<worktree>` stash snapshots (`mac/dirty/SIMPLE-main` = flex2D's 20 uncommitted files),
    `wip/state-eval`, `wip/prob-state-marg`, `wip/tree-ports-charts`, `wip/prior-parentpose`, `handoff-2026-11-07`, and hetero-tree's
    earlier guide tip `1d474d6e0`; restore with `git fetch <bundle> refs/remotes/fork/<ref>:refs/heads/<name>`.
  - `mac_untracked_20261107.tar.gz` (1.93 GB, 16,623 files + macOS `._` companions; extract with `--exclude='._*'`): every untracked
    file of every Mac worktree, this memory folder, side-session notes.
- Not archived: resume, GRFP application memories, session transcripts. The manuscript drafts (`doc/manuscripts/hetero_tree_2026_11`,
  worktree magical-babbage-701c7b; `doc/manuscripts/flex2D_cls_split_2026_10_03`) and the GRFP science notes are only in the
  tarball: the user said not to push the drafts (2026-11-07).

**Why:** later sessions run elsewhere; the Mac paths in older memories (`/Users/atifao/...`, `.claude/worktrees/...`) are gone.
**How to apply:** start from the guide on hetero-tree-wip; pull a ref from the bundle or a file from the tarball when needed.
Related: [[hetero-tree-branch-published-on-fork]], [[never-push-to-origin-without-permission]], [[site-mask-rework-direction-2026-11-07]].
