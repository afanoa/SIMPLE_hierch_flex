# 2D sub-class mask test (session mystifying-mcclintock, branch claude/tree-2d-subclass-mask off claude/tree-dist 9005a8c77)

Status 2026-11-01 evening: jobs running on batch2, work dir /mnt/beegfs/elmlund/afan/hetero_tree_20261028/subcls_mask/.
- 54793729 arms_tree: hetero_tree root (maxdepth 0, site_layer=no, KEEPMAPS) on the light arms set (160k non-junk rows of
  routes/toy55/simarms at the lat_arms refined poses), dist_sdisc/build12. -> arms/tree/a/1_hetero_tree/node_r
- 54793730 / 54793731 subcls_prep FIX=r131 / arms: node stack (r131) or oritab, direction bins (n/1300 as cheap_children),
  coupled-snapshot project, make_cavgs, cls_expansion pose mean 8/bin -> <fix>/prep/children.txt, bins.txt
- then subcls_rec (tree build, chunked reconstruct3D of every sub-class + every bin at the node poses) and subcls_eval
  (Python port of half_variance / falloff_mask / discover_sites; bin-centred covariance for the single-view children).
r131 comparator = route9/guardfix/runs/r131_md1_s1/1_hetero_tree/node_r (maps kept): its cheap-children mask covers 93 % of
the molecule support (tree.log), sites 3f1 / 1 / 3f2.
Scripts: repo scripts/cls_expansion_proto/tree_subclass_mask/ (copies in subcls_mask/scripts).

## Verdict (2026-11-01 ~21:30)
NEGATIVE on both fixtures: the bin-centred sub-class variance is noise peaking at the molecule centre (peak/null 3-7; cheap children 19,
peaking inside the closed arm / arm 4); matched-size recall of the arm regions 0.00 (r131 top 2-20 %) / 0.00 (arms, all sizes); discovery
48 fragments at r131 (best closed-arm coverage 0.22), 0-2 off-arm sites on arms; a 100-particle floor changes nothing; 7-20 core-h vs 1.
The cheap children's variance IS localized; only the 0.2-of-peak core rule fails (93-100 % masks) -> that is audit 6.1 (math-reworks session).
Do not wire sub-class children into cheap_children. Note: doc/implementation_notes/tree_2d_subclass_mask_test_2026_11.md (branch
claude/tree-2d-subclass-mask). Port validation: Dice 0.98/1.00, variance corr 1.0000, sites identical.
