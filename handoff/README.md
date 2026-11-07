# Handoff material, 2026-11-07

Made when the development Mac was retired. Start with the guide:
[`doc/implementation_notes/hetero_tree_testing_guide_2026_11_07.md`](../doc/implementation_notes/hetero_tree_testing_guide_2026_11_07.md).

| Path | What |
| --- | --- |
| `memory/` | Claude Code's project memory from the Mac (`MEMORY.md` = index, `hub-*.md` = topic indexes). Restore: guide section 10. The two GRFP application memories are left out. |
| `side_session_notes/<session id>/` | status notes and handovers written by side sessions (2026-10-30 to 2026-11-07): hetero_tree threads, the IgM fixture, OLE, port proposals, diffs |
| `batch2_scripts/` | small scripts and records from batch2 `/mnt/beegfs/elmlund/afan`: the gate harness (`hetero_tree_20261028/math_reworks/{launch_gates.sh,sst_gate.sbatch,readout.sbatch}`), scorers (`fixtures/scripts`), build scripts, gate run manifests and readouts, every thread's scripts (`math_reworks/{latent,shrinkB,anchor,mest,ppinfo}`), the sessions' handoffs (`hetero_tree_20261028/handoff/`), `routes/toy55` generators |
| `mac_untracked_list.txt` | every file in the batch2 tarball `/mnt/beegfs/elmlund/afan/mac_archive_20261107/mac_untracked_20261107.tar.gz` (1.93 GB): all of the Mac's untracked files, including result tables, embeddings and the manuscript drafts |
| `ARCHIVE.md` | the refs kept in the git bundle `/mnt/beegfs/elmlund/afan/mac_archive_20261107/fork_archive_20261107.bundle`: the Mac's other branches, uncommitted-work snapshots and parked experiments |

Not on GitHub: the manuscript drafts and result data (tarball), the Mac's other branches and worktree snapshots (bundle). Not
archived anywhere: the resume, the GRFP application memories, Claude Code's session transcripts.

One file from the tarball: `tar xzf mac_untracked_20261107.tar.gz --exclude='._*' "SIMPLE/<path>"` (macOS tar added a `._`
metadata companion per file).
