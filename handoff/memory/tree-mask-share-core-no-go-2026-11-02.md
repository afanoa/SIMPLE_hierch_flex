---
name: tree-mask-share-core-no-go-2026-11-02
description: "VERDICT 2026-11-02 (audit 6.1 relative rule): a node-mask core = the level set holding a fixed share (0.5-0.8) of the cross-half variance ENLARGES every mask (r131 root 93 -> 98-100 %, qtoy2 root 51 -> 56-66 %, depth-1 nodes 25-60 -> 65-99 %) because the variance is diffuse everywhere (top 5 % of voxels hold 0.11-0.35 of it); the 93 % r131 mask is a 13 % core in 14 pieces joined by the 3 px dilation + 6 px falloff, not the 0.2-of-peak rule; MFRAC/MDIL/MFALL stay"
metadata:
  node_type: memory
  type: project
  originSessionId: 433d3409-f4c3-4781-bfb8-d178578802df
  modified: 2026-11-02T04:47:32.656Z
---

2026-11-02, session tree-mask-discovery-stop (branch claude/tree-mask-discovery-stop off claude/tree-dist 9005a8c77), note
doc/implementation_notes/hetero_tree_mask_discovery_stop_2026_11.md, prototype scripts/cls_expansion_proto/tree_mask_discovery_stop/
(vmask_proto.py, Dice 0.98-1.00 against the trees' mask.mrc), data hetero_tree_20261028/mask_disc/proto/<node>/report.txt on batch2.
- Share of the positive cross-half variance in the top 5 % of support voxels: nulls 0.11-0.14, real nodes 0.16-0.35 (r131 root 0.25,
  arms 0.16, qtoy2 root 0.28); a level set holding 0.7 of it covers 25-41 % of a real node's molecule and 44-52 % of a null's; after
  the dilation and falloff the mask is 65-100 %. The 0.2-of-peak core (6-16 % of the molecule on real nodes) is smaller than any share
  core at S >= 0.5 on every one of 19 nodes tested.
- The arm coverage readouts cannot discriminate (a 95-100 % mask covers everything); the recB0 state maps after the focused stA0 are the
  localised ones (arms: top-5 % share 0.48, 0.2-of-peak mask 23 % holding all four arms) and the children already use them (oneflex).
- The variance map's autocorrelation half-lag is 3.5-5.6 px on nulls / cheap children and 6.6-11 px (18 px at 0.95 A/px) on real
  depth-1 nodes: a falloff from it would widen real nodes' masks; not adopted.
**How to apply:** do not propose share / quantile / relative cores for the node mask again; the lever for a smaller root mask is the
map source (recB0-like) or the dilation+falloff, not the core threshold. Related: [[mask-core-null-calibration-no-go-2026-11-01]],
[[tree-2d-subclass-variance-mask-fails-2026-11-01]], [[tree-svd-site-factors-and-variance-stop-2026-11-02]].
