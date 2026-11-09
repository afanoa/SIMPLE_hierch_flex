---
name: never-mirror-into-reused-roots-with-production-symlinks
description: "TRAP (2026-11-08): a replay harness rerun in a reused root whose node_* entries were symlinks into production wrote 1,599 self-links into production's subdirectories and an rm -rf through node_r13/recB0 deleted production's step directory; ln -sfn and rm -rf descend through a symlinked parent path component"
metadata:
  type: feedback
---

2026-11-08 19:13 EST, node-levers replays of f692_r17d_s1: harness v1 made `$R/node_r13 -> production/node_r13` (a symlink per
node at the replay root). Harness v2, started in the SAME roots, ran `mkdir -p node_r13` (no-op on the symlink) and
`ln -sfn production/node_r13/<entry> node_r13/<entry>`: the final component resolved through the symlinked parent to production's
real subdirectory, so ln created `node_r13/stageA/stageA -> itself` inside production (1,599 such links over all nodes); then
`rm -rf $R/node_r13/recB0` resolved through the symlink and deleted production's recB0 and recB contents (b.simple recovered from a
copy made minutes earlier; logs, node_rows and pcg/RESOLUTION files lost). The user's standing rule is never to delete cluster data.

**Why:** `ln -sfn` and `rm -rf` follow a symlink that sits in a parent path component; `-n` only protects the final component. A
reused replay root carries the previous layout.

**How to apply:** a harness that mirrors production (1) requires a fresh root (refuse an existing one), (2) refuses a destination
that is or sits under a symlink (`pwd -P` check against the root), (3) never runs `rm -rf` on a path that can contain a symlinked
component (remove a symlink with plain `rm` on the link itself), (4) is tested once on a scratch copy before touching anything
near production, and (5) asks the coordinator to make finished production trees read-only (chmod -R a-w) before replays start
(this session's tool permissions refused the chmod). Related: [[never-delete-only-read-or-move]],
[[never-overwrite-scripts-under-running-jobs]].
