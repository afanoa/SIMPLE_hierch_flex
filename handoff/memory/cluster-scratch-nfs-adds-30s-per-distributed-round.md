---
name: cluster-scratch-nfs-adds-30s-per-distributed-round
description: "TRAP: trees on /scratch/cluster_scratch (NFSv3, default attribute cache) pay ~24 s per distributed round to notice parts finished; 7-10 min of a 20-min small node; beegfs 4-10 s"
metadata:
  node_type: memory
  type: project
  originSessionId: 0f67cd42-16a3-4392-a859-a59fdac654fc
  modified: 2026-10-28T11:54:11.088Z
---

**Measured 2026-10-28** on the full692 build-7 base tree (site_states/full692/base_build_sst7_it3crop_emb_s1):
- Almost every distributed phase logs "4 PART(S) COMPLETED IN ~30 s", while its parts ran 4-10 s.
- The kept beegfs tree (ole_12707/6_hetero_tree) logs the same phases at 4-10 s.
- A small node (7-21k rows) runs about 22 rounds:

| Item | Per small node |
|---|---|
| Detection lag | 357-491 s |
| Tree-controller hand-offs (polling hetero_step.exit), median 7-29 s each | 100-260 s |
| Part compute | only 151-303 s |
| Queue wait | 30-45 s |

**Why:** /scratch/cluster_scratch is NFSv3 with no actimeo or lookupcache option. The directory attribute cache (acdirmin 30 s) keeps the part's sentinel file "missing", so qsys_ctrl's 1 s file_exists poll sees it up to 30 s late.

**Validated fix (recipe9, 2026-10-28):** `simple_sync_dir(dir)` before the done-file test, in qsys_ctrl `update_queue` and hetero_tree `step_status` / `dir_status` (poll_nfs_revalidate.diff).
- Probe lag: 10.6 s -> 0.07 s.
- qtoy2 maxdepth-1 tree: 56.8 -> 38.7 min, about 8.7 min per node, leaves identical.

**How to apply:**
- Before blaming compute or queue for slow small nodes on scratch, check the logged COMPLETED IN times.
- Builds without poll_nfs_revalidate.diff (recipe7 and earlier) still pay the lag. For those: run tree directories on beegfs, or ask the admins for lookupcache=positive.
- Details are in the side-session note small_node_wall_profile_2026-10-28.md.
- Related: [[m33-tree-wall-profile-2026-10-25]], [[ole-timing-triple-crops-verdict-2026-10-27]].

**beegfs space (2026-10-28):** Default pool (our dirs) back to ~94 % full, ~460 GiB free on each of 39 targets (~18 TB) after the user deleted the Oct 27 runs; check with `beegfs-ctl --listtargets --storagepoolid=1 --spaceinfo`.
