---
name: package-tree-tests-2026-10-25
description: "2026-10-25 whole-tree tests for the coordinator: the m8c75t design trees (one motion per node + split-half gate + inherit) vs same-build baselines, and the m8c71t multi-factor package trees kept as the comparison; where each lives"
metadata:
  node_type: memory
  type: project
  originSessionId: 58de2c8f-ad53-44fb-b8a2-05aadff16ec4
  modified: 2026-10-25T17:15:00.967Z
---

2026-10-25, for the coordinator "2D to 3D heterogeneous state translation". Dirs are under /mnt/beegfs/elmlund/afan/.

(1) The intended design: the trees started on m8c75t, and all 12 now run on m8c77t.
- m8c77t = m8c75t + pole-rule fix (A) + the symmetry-mate fold. Fix (A): a kept pole keeps exactly the seed rows it was tested on; a rejected pole's rows go to their best kept pole at posterior >= 0.5, else the rest. The fold drops a site matching an earlier one at Dice >= 0.99 after a group rotation.
- Keys: node_motion=yes motion_multi=no motion_gate=yes motion_inherit=yes strongest_pole=yes two_poles=yes, SIMPLE_TREE_REENTER=yes, 4e on.
- Baselines: the same build with every one of those keys off. The flex component sort changes what site_dirs=1,3 seeds from, so a baseline is only valid on the same build.

| tree | design | baseline |
|---|---|---|
| OLE s1 / s2 (resumed) | ole_12707/10_tree_design_s1, 11_tree_design_s2 | 12_tree_base75_s1, 13_tree_base75_s2 |
| m33 full (692,392 rows, junk_screen=no, maxdepth 7): RELAUNCHED fresh on m8c80t (the split-screen fix), 54253316 / 54253317 | stree_port/full692_design80_s1 | full692_base80_s1 |
- 2026-10-25 13:56: design80 CANCELLED by id (its root split on site 6, which the scorecard calls a per-particle quality axis);
  full692_design83_s1 (m8c83t, nodepar 24) cancelled at its root mean step. base80 RESUMED in place on m8c83t at nodepar 8
  (stree_resume_pkg 54278813, PF=rnasep.simple, INFILE=rows01.txt) = the concurrency timing test.
| cdisc M (restarted) | stree_port/fulltoy_m77Mdes_cd1 | fulltoy_m77M_cd1 |
| CAHRA M s1 | stree_port/arm_M77d_c11_s1 (restarted) | arm_M75b_c11_s1 (resumed) |
| IgM M (restarted) | stree_port/fulltoy_m77Mdes_igm | fulltoy_m77M_igm |

- The m8c75t cdisc, IgM and CAHRA-design dirs (fulltoy_m75*, arm_M75d_c11_s1) were abandoned when those trees restarted.
- The mirrored swing at the m33 root (replay stree_port/gate/m33_root_m79, re-split without the screen in m33_root_m79_split80_nojunk)
  gives K* 3: 372,159 middle / 251,181 B side / 69,052 A side, re-test z 6.80. The deviation maps put B at bin 1 (+0.89) and
  A at bin 5 (+0.85), so B is the other extreme. With the screen on, A was junked and K* was 2.
- The m8c77t m33 pair (full692_design_s1 / full692_base75_s1) was cancelled by id on 2026-10-25. The split screen had junked the root's site-6 extreme (124,999); the split-only re-run with junk_screen=no gives K* 2 (597,146 / 95,246), stree_port/gate/m33_root_split80_nojunk.
- The Gamma-pairing patch is out of these 12. On the OLE s2 r2 fflex step (m8c76t vs m8c74t) it gave Spearman(Gamma, var z) 0.927 -> 0.988, lower -2logL/N, components 9/10 regularized again. Its 5-state partition matched the old one at ARI only 0.21, with 15 A deviation maps 0.59-0.81 under the best pairing.
- m8c78t = m8c77t + the Gamma patch, as a snapshot only; the local worktree and ~/SIMPLE_m8c stay at m8c77t. Its OLE design pair is ole_12707/14_tree_design_g_s1 and 15_tree_design_g_s2 (54234679/80, same keys). Score with ole_12707/scripts/tree_pair_verify.py against 10_/11_ and 6_/7_ (unit ARI 0.18).

(2) The multi-factor package (product cells), the comparison:
- On m8c71t: the seeding rule is a cell's core intersection when it holds minpop, else the whole cell.
- Resumed in place: stree_port/full692_pkg_s1, ole_12707/8_tree_pkg_s1 and 9_tree_pkg_s2, fulltoy_m69Mpkg_cd1.
- full692_pkg_s1 (controller 54227856) was CANCELLED by id on 2026-10-25 and its dir kept: product cells are superseded, and the split
  screen had junked 107,655 at its root.
- Restarted from scratch: arm_M71p_c11_s1 and fulltoy_m71Mpkg_igm, because their m8c69t hand-backs had cores under minpop.
- Their baselines: fulltoy_m69M_cd1, arm_M69b_c11_s1, fulltoy_m69M_igm (m8c69t).

Floors with matched fakes dropped (2026-10-25, m8c75t): IgM 0.098, cdisc 0.093, m33 0.061, OLE 0.054. Every non-root real site axis passes except IgM r2 site 9. The m33_design replay's second candidate per node is a SKIPRANK duplicate of site 1; see [[gate-replay-skiprank-duplicates-site-1]].
- CLOSED: making a reproducible floor fake a candidate. CAHRA's check: the OLE fakes (0.30-0.32) match the other seed's missed axes at |r| <= 0.18 and are partial copies of known axes. The floor rule stays as built.
cdisc pair on m8c77t (stree_port/gate/toyscore/cd1_series_*):
- Design unit NMI 0.369 vs baseline 0.323. j2 junk recall 0.62 vs 0.14. Closed recall 0.19 vs 0.17.
- REGRESSION, older than the design: the closed state was lost between m8c58 (site-layer closed pole AUC 0.94, closed precision 0.78,
  13 leaves) and m8c69t (AUC about 0.51, 5-8 leaves). RULED OUT 2026-10-25 (m8c81t, gate/toyscore/cd1_4e_*): split_ntests=1 and
  node_proj=no each give a tree identical to m77M_cd1 (every cd1 re-test z is far from both z_crit; node_proj is I/O only). Seed 2
  (m81Mnt1_cd2) finds the closed pole again (r11_s4, AUC 0.937). The cause is another m8c58 -> m8c69t change.
CAHRA pair on m8c77t (stree_port/gate/toyscore/cahra_series.txt):
- The design lost the open state: recall 0.527 vs the baseline's 0.857, NMI 0.12 vs 0.29.
- Both roots seeded from site 1 with near-identical seeds. The design's split cross-half EM routing mixed the stage-B states (15,735 /
  14,265) where the baseline's followed them (23,106 / 6,894); see [[split-routing-em-vs-stageB-2026-10-23]].
Gate test replays are in stree_port/gate:
- *_floor75: the floors re-measured with matched fakes dropped;
- m33_design: the m33 floor plus the r2 read;
- *_fb_b: the fallback series on m8c73t.

Tools:
- Scripts in stree_port: stree_gate_replay2.sbatch (STOPNODE, SKIPRANK, FROMFFLEX), stree_resume_pkg.sbatch (resume in place), stree_toyscore.sbatch.
- gate_matched.py re-scores gate maps.

**Why:** the user's goal is "fewer and better states". Multi-factor product cells over-fragmented: m33 r2 gave K* 11, IgM r1 K 19.
**How to apply:**
- Compare each design tree only with its m8c75t baseline. Judge the OLE pair by cross-seed ARI against the 6_/7_ unit ARI of 0.18.
- Never resume a tree on a build whose hand-back would differ from the one it already ran. Check the seeded counts first.
