---
name: hetero-tree-4p5m-12h-serial-profile-2026-10-23
description: "TARGET: full tree on ~4.5M clean particles in ~12 h wall. The 692k root is 60 % master-only time: flex_pca's all-particle contrast embed in its master (0.42 ms/particle), one-job site_rank and state_split; plan v1 = 12.4-14.7 h with the approved changes only"
metadata:
  node_type: memory
  type: project
  originSessionId: a4258345-b249-48ff-a0d3-d33e99ae33ac
  modified: 2026-10-24T02:28:44.275Z
---

2026-10-23. The user's target: the full data (~4-5M clean particles) as one hetero_tree in ~12 h of wall. The coordinator and the
port asked me for a step profile. Tables and scripts: batch2 `rnasep_klem/fast_tree/prof/`
(tw1_timeline.txt, tw1_profile.txt, timeline.py, profile_tree.py). Design note section "Step profile of the 692k tree and the
4.5M plan".

- Twin tw1 (m8c20, 3+3): the root (692k) takes 2.32 h, 60 % of it master-only work that more parts cannot shorten.
  - flex / fflex: ~920-960 s each. site_rank: one 32-thread job, 787 s. state_split: one job, 640 s.
  - refine3D masters: 85-115 s per iteration.
- flex_pca's master embeds EVERY particle after its part rounds ("CONTRAST EMBED SECONDS"): linear, 0.42 ms per particle, so
  ~32 min per flex step at 4.5M. Distributing it is the user's code (flex_pca), so it goes to the user first.
- Size scaling, from node sizes 47k-692k:
  - site_rank sublinear (787 / 532 / 437 s at 692 / 300 / 183k);
  - state_split about linear above 200k;
  - the root-only sigma2 bootstrap is capped at 25k particles, so it stays flat. Don't extrapolate every master term linearly.
- Plan v1, joint with the port: with only the approved changes (node_proj node-sized projects, memory sizing) the projected
  wall is 12.4-14.7 h. Then stA0 one pass -0.6 h; a split subsample -1.8 h and a site_rank cap -0.7-1.7 h (both method changes,
  so truth checks first); the flex embedding -1.2 h.

- USER DECISIONS 2026-10-23:
  - NO result-changing shortcuts (no split subsample, no site_rank cap).
  - The port distributes flex_pca's contrast embedding (identical output; the user reviews the diff).
  - The port runs a merge-up replay at 32 threads, which settles the end phase: if merge-up doesn't parallelise, it is ~6.8 h at
    4.5M.

- Plan v2 final (2026-10-23; node_proj on by default; one pass for stA0 / A / B; the embedding distributed):
  - 11.6-14.6 h in total;
  - with merge-up near-linear in threads: 10.1-13.1 h (merge-up measured 1,962 s at 32 threads at 692k, only 1.9x for 4x);
  - with the split's image load fixed as well: 9.1-12.1 h. The 692k split replay took 900 s: half of it loading 85 GB at
    190 MB/s, half compute.
- TRAP: a root split ran 5 h+ because its job inherited OMP_NUM_THREADS=8 with nthr=32. Check the environment of tree step
  jobs.

**How to apply:** for wall-time work, profile the master's serial phases (BENCH plus the sacct round gaps) before adding parts.
Refresh the plan with the m33 run (full692_m33_s1) numbers. Related: [[hetero-tree-schedule-cuts-fail-unconverged-drift-2026-10-22]], [[flex-delivered-state-maps-one-process-2026-10-24]].
