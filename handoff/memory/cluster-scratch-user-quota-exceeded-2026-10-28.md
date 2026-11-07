---
name: cluster-scratch-user-quota-exceeded-2026-10-28
description: TRAP: /scratch/cluster_scratch has a hidden per-user quota; at 15:09 2026-10-28 every atifao write failed ("Disk quota exceeded") with the FS 79 % full and `quota` showing none; all sessions' jobs died in seconds
metadata:
  type: project
---

2026-10-28 ~15:09: every write under /scratch/cluster_scratch/atifao failed with "Disk quota exceeded" (touch at the top level too),
while `df` showed the filesystem at 79 % (55T free) and `quota -s` said "none": the filer enforces a per-user/tree limit that no
login-node tool reports. Worker jobs of every session (route A, recipe9, recipe11 gates/subtrees, recipe12 and d3513ff70 builds,
equiv_port) FAILED in 2-19 s; tree controllers kept RUNNING but could not write their logs (logs freeze at the failure time).

**Why:** after beegfs filled ([[beegfs-transient-write-failure-2026-10-27]], [[test-runs-moved-to-delete-folder-2026-10-28]]) all runs moved
to cluster_scratch; full 692k trees with crops still write 0.5-2 T each, and several sessions share the one user quota.

**How to apply:** before a launch, `touch` a file in the run folder (a failing touch = quota); grep step logs for "quota" when jobs fail
in seconds; size trees to the shared quota; freeing space needs the user (no rm by me; mv within the quota frees nothing).
See [[never-delete-only-read-or-move]], [[site-pole-closed-state-replay-2026-10-27]].

Size when it hit (du, 15:3x): ~5.1 TB = site_states 2.9 TB (route A sst7 1.93 TB, build-7 base 631 GB, recipe9 692k 315 GB and
growing: recipe9's prune_partials deletes nothing), ole_repro 1.4 TB, recipe4 648 GB. So the limit is about 5 TB. Parts fail with
exit 0:53 or 1:0; distributed masters whose parts failed stay RUNNING forever: cancel them by job ID and rerun once space is back.
- 16:3x user: "can you just go back to beegfs". New root /mnt/beegfs/elmlund/afan/hetero_tree_20261028/{fixtures,recipe12,equiv,
  from_scratch}; beegfs Default pool 318 TiB, 11 TiB free (96.4 %), shared -> small per-session budgets, cleanup on. Full route A /
  recipe9 trees not resumed (closed subtrees replace them). Sessions copy keep sets to from_scratch/ before the user wipes scratch.
  My scancel of a hung controller was blocked by the permission classifier: cancelling running trees is the user's call.
