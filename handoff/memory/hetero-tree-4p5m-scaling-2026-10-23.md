---
name: hetero-tree-4p5m-scaling-2026-10-23
description: "hetero_tree at 4-5 M particles (user target ~12 h): NOT on m8c33 - master-sized step projects make the controller serial (~50-60 s per step at 4.5 M rows, x nodes x 12 steps), site_rank and state_split are single jobs over their 62.5 GB ask; fix = node-sized node projects + memory sizing"
metadata:
  node_type: memory
  type: project
  originSessionId: 58de2c8f-ad53-44fb-b8a2-05aadff16ec4
  modified: 2026-10-24T16:09:19.547Z
---

2026-10-23 rough projection for the coordinator (the user: "if 4-5 M particles can be done in 12 ish hours I can accept that"; the grow-on-subset-then-assign plan was PARKED). Firm numbers wait on full692_m33_s1 (job 54157123).

Measured basis: the cancelled m8c20 twin root at 692k took 2 h 19 min and 90 worker core-h.
- site_rank is one 32-thread job: 13.1 min / 18.6 GB at 692k, 8.9 min / 8.9 GB at 211k. That is about 7 min + 0.53 ms per particle and about 20 GB per M.
- state_split is one job: 10.7 min / 17.9 GB at 692k, 4.1 min / 8.8 GB at 211k.
- Both request 62.5 GB (hetero_step.sh --mem=64000), so both OOM above about 3 M rows.
- Distributed steps get one part per 20k particles, capped at nparts (NPARTS_MIN 4, PTCLS_PER_PART 20000).
- Controller memory was 1.05 GB at 932k rows.
- The controller's step transitions (write the step project, cut to node rows, expand back) take about 10 s at 932k rows (488 MB).

**Why it matters:** every node's step projects are master-sized, so each transition costs about the same at any depth. At 4.5 M rows (2.4 GB) that is about 1 min per transition, serial in the controller. With 60-120 nodes x 12 steps, that adds up to 10-20 h. Disk would be about 3-4 TB per tree (sigma2_state.bin is master-sized too).

**How to apply:** before any multi-M tree:
- Make node projects node-sized: extend node_rows from per-step cuts to per-node cuts, with a row map back to the master for the leaves' labels.
- Size site/split memory from the row count.
- Run nparts=150, nthr=8, nodepar=8, with workers on -p csbdevel,norm.
Projected result: about 8-10 h and 3,500-5,000 core-h. Expect 1-3 more levels and about 80-150 nodes.

Measured later on 2026-10-23:
- Merge-up at 32 threads (m8c34, 692k): 1,962 s, 19.2 GB, against 3,790 s at 8 threads. That's only 1.9x, so about 3.5 h at 4.5 M; the end phase is about 4 h.
- node_proj (m8c34, key, default no): the controller per transition went from 8.3 to 1.7 s (r3, 932k master), and disk fell about 36 % per node.
- BUT node_proj is NOT result-identical. Flex's paired fits split by mod(pind, 4) on the PROJECT ROW, which node-sized projects renumber (it also breaks flex's e/o balance). The same-seed floor is exact (IgM no1/no2 ARI 1.000); yes vs no gave ARI 0.62 on IgM (same structure) and 0.03 on RNase P r3.
- The fix proposed: a stable half key (the input row via the node's rows.txt), which needs the user's OK since it touches flex.
- Seed noise (m8c34, RNase P r3 151k, node_proj=no, seed 1 vs 2): ARI 0.195, and different tree shapes. Seed 1 splits 123k / 28k
  into 4 leaves; seed 2 does not split (126k leaf + 25k junk leaf); the mean maps already differ. So yes vs no (0.03) can't be
  separated from seed noise on RNase P; only a same-seed identity pair (IgM) can validate the half key.
- Half key validated (m8c37k, IgM same seed, no vs yes): ARI 0.9995 (was 0.62); root and r1 make identical decisions (seeds,
  split children); float noise of ~1e-6 is amplified by r2's fflex (|dz| up to 0.34) -> 28 of 200k final labels differ; maps 1.0000.
  Not bitwise identical (summation order, like an nparts change). User approved; committed as m8c40 (3e8395c95, gated,
  snapshot ~/SIMPLE_stree/snap/m8c40). node_proj default still no.
- RNase P same-seed check (m8c40, r3, maxdepth=1, no vs yes): ARI 0.867, NMI 0.731, leaf maps 0.9998, same tree shape. fflex
  equal to float level; 227 seeds and 2,740 split children differ, at the real/junk leaf boundary. FAILS the ARI >= 0.99 bar, so the
  default was NOT flipped (user's call). The "parts cut by row ranges" cause was WRONG: refine3D, rec3D and prob_align balance
  parts on state > 0 (split_nobjs_active), and flex cuts over its own pinds. The 2026-10-24 diagnosis: flex_pca is not
  reproducible run to run in ANY layout, and r3 yes40's fflex is byte-identical to the m8c34 node_proj=no run's. So the
  0.867 is flex noise, and no same-seed pair meets 0.99 until flex is deterministic
  ([[flex-pca-ab-needs-simple-seed]]). The subtree identity check can't be bitwise for the same reason.
- 2026-10-24, from the speed-up session: flex_pca's delivered state maps run as ONE in-process reconstruct3D
  (deliver_state_maps passes no nparts). That is 34-47 % of each m33 fflex step below the root; nparts=12 cut r1's 599 s
  to 189 s (maps corr 1.000000000), about -2 h over L1-L5 at 4.5M. It's a two-line flex_pca change, so it goes to the
  user. The root's flex/fflex/mean maps use only updatecnt > 0 rows: rnasep.simple's "flat 25k" is its 25,000 updated
  rows, not a cap. A 4.5M input with all (or no) rows updated would rebuild 4.5M rows in one process, about 80 min per
  flex step, unless that change lands.
- 2026-10-24 levers (692k m33 root, 32 threads, HEAD m8c55t):
  - split 507 s (load 117 s I/O-bound, shuffle null 68 s, about 320 s of few-map phases that barely scale);
  - merge-up 31:17;
  - site_rank 809 s, pole EMs scaling about 3x per 4x threads, so no lever there.
  - Lever (c), state_mergeup accumulate_jobs (fewer maps than threads: chunked parallel accumulation, sum_reduce in slot
    order): split 237 s, merge-up 19:23, numerically close (ARI 0.9997 / 1.0000). Uncommitted:
    .claude/patches/mergeup_parallel_accumulate_on_16e52fbdd.patch (user decision).
  - beegfs, measured: one stack 625 MB/s at 1 stream, about 1.5-2.6x at 2-4 streams, 8 is worse; 8 different stacks 2.35 GB/s.
    RNase P has 10 stacks of 96k images, so every batch sits in one stack.
  - TRAP: loader-local striping (8 stretches per batch) was output-identical but loaded 2.6x SLOWER (115 -> 305 s) through
    stack_read_session (single prefetch producer). The across-file gain needs a shared-reader change.
- TRAP: hetero_tree resume scripts must use projfile=<the input project's basename, copied into 1_hetero_tree> (rnasep.simple), not root.simple; stree_f692twin/stree_f692m33 had this bug. A resume with infile keeps l_subtree.
- Safe resume: stop the controller by job id only while its in-flight step has a hetero_step.exit.job start record and no transition is running. A step submitted but not yet started would be resubmitted.
- Resume of a LOCAL-queue tree (OLE 6_hetero_tree, 2026-10-24: m8c33 local -> m8c55t slurm, job 54206504, ole_12707/scripts/ole_tree_resume.sbatch).
  - The cancel kills local steps. A killed step keeps its start record with no exit, so the resumed controller waits on it forever: mv that step dir aside. start() then reruns it from prev_step.
  - Also set aside a node project written at the cancel.
  - Re-stamp qsys_name=slurm in every project the controller copies from: the tree project, node_*/n.simple, node_*/stageA/r.simple, node_*/stageB/b.simple. Each step takes its queue from its project's compenv.
  - Result: finished in 4 h 49 min wall, all 37 nodes, no aborts. 31 leaves over 80,598 particles, merge-up 14 units.
- Capping a running tree at depth D and finishing it (m33, 2026-10-24, job 54214929, stree_port/stree_f692m33cap7.sbatch):
  - cancel by id;
  - in the job, drop every node of depth > D from the node table, and set the depth-D nodes that had children back to step 12 (they re-decide into leaves);
  - drop the leaves of removed nodes, and of those depth-D nodes, from the leaf table;
  - reset FAILED nodes to step 4 and set their stA0 aside;
  - keep *.before_cap copies, then resume with maxdepth=D.
  - TRAP: an apostrophe inside a single-quoted ssh command ends the quote. The rest then runs LOCALLY and the heredoc lands half-written. Write such files locally and scp them.

Related: [[tree-contrast-and-identity-junk-2026-10-23]], [[size-job-memory-to-measured-rss]], [[batch2-400-core-limit-and-worker-timeouts]].
