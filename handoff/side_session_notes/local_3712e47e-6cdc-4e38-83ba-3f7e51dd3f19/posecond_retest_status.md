# Pose-conditioned split re-test prototype (session amazing-tesla, worktree amazing-tesla-c2b47f, branch claude/amazing-tesla-c2b47f = claude/tree-dist a9b6507a7 + commit e8da598cd)

Work dir on batch2: /mnt/beegfs/elmlund/afan/hetero_tree_20261028/posecond/ (scripts/, cases/<case>/{full,stable,turned}/ half maps + stats.tsv,
counts.tsv, flags.txt; results/ tables; cases_bad_selection/ = the first attempt whose prg=selection aborted; logs/). Repo: scripts/hetero_tree_proto_posecond/
(+ results/ TSVs), note doc/implementation_notes/tree_pose_conditioned_retest_2026_11.md. Memory tree-pose-conditioned-retest-no-go-2026-11-02.
The branch claude/tree-dist can be fast-forwarded to e8da598cd (one commit on top of a9b6507a7). Nothing pushed.

## 2026-11-02 ~16:10 (batch2 time): DONE, verdict NO-GO
- 14 nodes / 30 pairs: mic_het root + depth 1 at both seeds, micnull root (stage-B states), qtoy2 root + r3 (quality), r131 root + r1,
  real PfCRT roots x2, 692k r17d root (r4 vs siblings). Per child: full / stable (angle <= 30) / turned subsets reconstructed at stage-B
  poses (reconstruct3D nstates=K box_crop 88 ml_reg=no), draw-free T at the node's logged n + guard L/L_ref (Python ports; the full
  condition reproduces the tree: T .966 vs logged .955, L/L_ref 1.31 vs 1.36).
- T_stable / T_full = 0.99-1.01 on EVERY pair, truth irrelevant: mic_het real root .970/.966, micnull pose null .876/.880, PfCRT
  .975/.976 and .987/.989, depth-1 real .910/.908, quality .167/.165, r131 .503/.503, 692k r4 pairs .955-.963 vs .957-.964. Turned
  subsets carry it too (.89-.99). Guard unchanged per condition (micnull .24/.24/.28). Sphere mask: same ratios. Within-child
  stable-vs-turned T .41-.79 < between (.88-.99) except the weak mixture splits.
- Why: stage B assigns particles by map fit; any subset of a child reproduces its map; the turn is a correlate of the worse-fitting
  population. Keep dflip report-only; no pose gate; the guard stays the only thing separating micnull from mic_het; PfCRT's
  partitions are reproducible population differences made local by the user mask (open question about the populations, not poses).
- r4: its difference from r1/r2/r3 persists among stable particles (T .963/.963/.955, L/L_ref .985/.997/.998) = population (quality)
  difference, not alignment; the stage-B corr gap persists in the stable subsets (.3904 vs .3989/.3958/.3940) -> quality flag is the lever.
- TRAP for anyone relabelling node projects: prg=selection aborts on node-sized projects (STK/PTCL2D/PTCL3D consistency check) and
  the reconstruction silently proceeds with the old states; posecond_setstates.py writes float32 slot 23 (I_STATE) of segments 3/6 in
  place and the job verifies by a dump. Cost 26.5 core-h (jobs 54870765, 54871205-17, 54873616-24 on norm).
