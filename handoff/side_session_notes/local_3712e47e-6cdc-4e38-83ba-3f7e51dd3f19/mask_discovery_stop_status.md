# hetero_tree mask / discovery / stop (session tree-mask-discovery-stop, worktree tree-mask-discovery-stop-7503a4, branch claude/tree-mask-discovery-stop-7503a4 off claude/tree-dist 9005a8c77)

Work dir on batch2: /mnt/beegfs/elmlund/afan/hetero_tree_20261028/mask_disc (proto/<node>/report.txt, builds/{base,m23,m23d}, gates/<arm>_<ds>[_sN], scripts/).
Note on the branch: doc/implementation_notes/hetero_tree_mask_discovery_stop_2026_11.md. Prototype: scripts/cls_expansion_proto/tree_mask_discovery_stop/vmask_proto.py
(port of half_variance / falloff_mask / discover_sites, Dice 0.98-1.00 vs the trees' mask.mrc, identical site lists to the Fortran).

## 2026-11-02 ~02:30 (batch2 time)
- Item 1 (share-rule mask core + acf falloff): NO-GO at the prototype on 19 nodes. The cross-half variance is diffuse everywhere (top 5 % of
  support voxels hold 0.11-0.14 of it on nulls, 0.16-0.35 on real nodes); a level set holding 0.5-0.8 of it is LARGER than the 0.2-of-peak core
  on every node, so the mask grows (r131 root 93 -> 98-100 %, qtoy2 root 51 -> 56-66 %, depth-1 nodes 25-60 -> 65-99 %). The r131 93 % mask is
  a 13 % core in 14 scattered pieces joined by the 3 px dilation + 6 px falloff, not the core rule. No Fortran change, no gate; MFRAC/MDIL/MFALL stay.
- Item 2 (SVD site factors, DFACT_SHARE 0.2 + Fisher z >= 3 on the odd half's projection, n_eff = half the voxel count x the 15 A band share):
  implemented (discover_sites; DFACT_COS/DFACT_ITS/top_dirs gone). Prototype: identical sites where the 2-means agreed, one fragment site fewer
  at the r131 root (10f1/10f2 -> 10), arms 5 -> 4 (same arm recall 3/4), nulls unchanged (20 / 7 / 0), qtoy2 root unchanged (6).
- Item 3 (diffuse-variance stop): stage-B variance only; level set holding 0.7 of it covers 33-41 % of the molecule on the three nulls, 31 % on
  the refused qtoy2 r1, <= 22.5 % on every real split (ole r2 the closest). Coded at the split step (MSHARE 0.7, MDIFF 0.3): diffuse -> leaf, no
  state_split. Saves only the split job (0.15-0.72 core-h per null root, ~1 % of a node); cannot shorten deep 692k lineages (they are localized).
- m23 gates (items 2+3, old ranking): nulls one leaf each (stop fires at the root: 40.7 / 37.8 / 34.6 %); r131 seed 2 corner KEPT (r1 s02 0.99/0.52,
  T 0.195 vs base 0.154) although the root site order flipped (site 1 margin 12.2 -> 5.0, site 3 7.0 -> 7.5; site 1 same pole, mask Dice 0.95);
  qtoy2 seed 1 WORSE (NMI 0.667 vs 0.712, closed combo 18 F1 0.59 vs 0.87) through site_states extra states from the changed factor masks.
- Running: m23d arm (items 2+3 + cherry-picked variant D ranking dd7fa0f45 = 6b4697c38): r131 s2 + s1, qtoy2, three nulls; m23 ole; a qtoy2 seed-2
  pair (base build_dump vs m23) to tell the qtoy2 change from the seeding lottery. Readout -> gates/final_readout.txt.
- Draw sd: the split test's subsample draw is seeded by the constant RND_SEED in state_mergeup (SIMPLE_SPLIT_EMNULL_SEED seeds a log-only switch);
  I did not touch that file, so verdicts carry the single-draw T and threshold only.
- Commits pending the readout: item 2 (discover_sites), item 3 (stop; may be reduced to the log line if the gates only show no harm), note + scripts.

## 2026-11-02 ~03:10 (batch2 time): DONE
- Commits on claude/tree-mask-discovery-stop-7503a4 (on dd7fa0f45 = cherry-pick of the rework's 6b4697c38 onto 9005a8c77), nothing pushed:
  413c0ce39 item 2 (discover_sites: SVD factors kept at share >= 0.2 and cross-half Fisher z >= 3, DEFAULT); fcc699788 item 3 as a log line
  only (variance_spread at every split step, no leaf rule); 9b26a8076 the note (sections 4.2-5 filled) + the prototype scripts.
- Item 2 verdict GO (user may drop 413c0ce39): nulls 3/3 one leaf; OLE 10 leaves as base, landscape within base's null (W1 .002), 16 sites
  not 22; qtoy2 under variant D NMI .796 vs .764 with the SAME root partition (gain = r1 seeding, combos 1/2 F1 .73/.59 -> .85/.80, net
  of an r5 draw over-splitting combo 101 .94 -> .71); seed-2 pair under the old ranking .733 vs .715 (seed 1 .667 vs .712).
- r131 corner under D (root site 1 first at both seeds, root partitions within 130 particles of the reference, depth-1 seeds identical =
  fflex states): seed 2 kept (r1 T .147 own draw; 2/5 draws pass, base 2/5; T_df .126 vs .117). Seed 1: the corner sits in r2's stage-B
  map (.97/.53) and the single-draw test refuses it (0/5 draws; the rework's own sub4 r2 passed 1/5, its tree's draw was the lucky one);
  T_df .122 vs .121, z 6-10 both. So the corner is a split-test draw in EVERY arm at these ~20k nodes, not a site-rule effect; the
  draw-free statistic (their 9ee5d5c72, not default) passes all of them. Replays: mask_disc/draws (build_df + SIMPLE_SPLIT_SEED 1-4).
- Item 3: stop fired on 3/3 null roots (33-46 %), never on a real node (7 trees), saves ~1.1 core-h total -> log line only (feature freeze,
  science before speed). Final build builds/final: fast gate 15/15; gates/final_qtoynull md0 one leaf, line printed (40.7 %).
- Cost: m23 arm 151 core-h, m23d 84, qtoy2 seed-2 pair 48, draws 3, prototypes 4.5, builds 6 (about 300).
- Open questions for the user: note section 5 (keep item 2 as default or hold until the split test is settled; add the stop or not).
