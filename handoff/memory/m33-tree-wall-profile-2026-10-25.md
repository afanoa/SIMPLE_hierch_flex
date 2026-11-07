---
name: m33-tree-wall-profile-2026-10-25
description: "692k RNase P hetero_tree wall profile (sacct by WorkDir): serialization-bound, not cluster-bound. nodepar 2 + the design's sequential site ranking dominate; user target 5-7 h per 700k tree; ranked levers"
metadata:
  node_type: memory
  type: project
  originSessionId: 58de2c8f-ad53-44fb-b8a2-05aadff16ec4
  modified: 2026-10-25T23:21:43.319Z
---

2026-10-25, for the coordinator, after the user's priority: "5-7 hrs is the goal for 700k particles because otherwise I can't really test/develop".

Method: one `sacct -u atifao -S .. -X -n -P -o JobID,JobName,Submit,Start,End,ElapsedRaw,AllocCPUS,State,WorkDir`, mapping each job to node and step by WorkDir. The script is the scratchpad's prof/tree_profile.py (pure python, runs locally).

Findings:
- Queue waits were about 0 and controller gaps 7-11 s.
- full692_m33_s1: 28.4 h wall, 987 core-h. 5.5 h of it was the old m8c33 root split, which takes 0.1 h today. Most of the rest was nodepar-2 slot waits: 1.5-3.5 h between levels on the critical path, about 14 h in all. Its merge-up took 1.5 h (pre accumulate fix; about 19 min since).
- m8c80t roots: the baseline takes 1.57 h. The design takes 3.23 h, of which its 7 site_ranks run sequentially at 0.18-0.40 h each, 1.93 h in all.
- Per step, at 692k / 200-500k nodes: flex+cheap at the root (infile) 0.23; fflex 0.17-0.20 / 0.16; stA0 0.1; stageB 0.12-0.18 / 0.1; split 0.08-0.11 / 0.07.
- Small nodes (10-50k) take 0.25-0.5 h for about 10 sequential SLURM steps.

Wall model (level sizes 1/2/4/8/16/27/31/30/10, node walls 1.57/0.9/0.6/0.5/0.45/0.4/0.35/0.3/0.25 h): the baseline takes about
26 h at nodepar 2, 9.2 h at 8, 6.4 h at 16 and 5.3 h at 32, plus the merge-up. The design with only concurrency takes about 10 h, about 7 h
with concurrent site ranking, and about 5.5-6 h with the lean subsample site layer plus the flex levers.

Ranked levers:
1. nodepar 16-32 with workers on csbdevel,norm (a key, no code change).
2. subtree_max jobs for nodes of 50k or less (built, never timed).
3. Concurrent site ranking in the design.
4. Lean subsample site ranking.
5. flex: state-map delivery with nparts (the user's code, 599 -> 189 s) and fits on 100k or fewer rows.
6. merge-up pair tests in parallel.
7. The root mask from fflex maps instead of flex+cheap (0.23 h).
8. A calibrated motion_floor instead of the nulls.
9. Motion bins only for the strongest candidate.
Not levers: the 50 shuffles, stA0's iterations, controller latency.
Status 2026-10-25 14:00: (3) and (9) are built in m8c83t (results identical). (5a) was already in the numbers: a027c6d03 is active, so
the tree's fflex ran nparts=12. (6) was skipped (one shared-memory merge-up job). Run A (design83, nodepar 24) was cancelled under the
freeze ([[tree-design-feature-freeze-2026-10-25]]). The timing test is base80 resumed at nodepar 8 (54278813). It times nodepar only:
the baseline uses neither the site batch nor motion.

- base80 RESULT (finished 2026-10-25 19:12): 10.02 h, 768 core-h, 144 nodes. A hybrid: nodepar 2 (m8c80t) for 0-4.73 h, 8 (m8c83t) after.
  - Depth ends: 1.57 / 2.87 / 4.07 / 5.23 h, then depths 4-7 together at about 8.7 h (nodepar-bound). Node work ends at 8.78 h.
  - End phase 1.24 h: final_maps 18 min, serial in the controller before the merge-up is submitted; merge-up 54 min (top level 19).
  - The clean timing run is base90 (54300456, m8c90t, nodepar 16).

**Why:** the user needs a fast test loop; every design feature has to earn its wall time.
**How to apply:** time trees from sacct by WorkDir, not from logs. Check concurrency (nodepar) before blaming compute.
Related: [[hetero-tree-4p5m-scaling-2026-10-23]], [[package-tree-tests-2026-10-25]].
