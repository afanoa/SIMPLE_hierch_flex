---
name: quality-matched-null-core-score-2026-10-28
description: "SIMPLE_SPLIT_QNULL (quality pseudo-splits, qnull_v4.diff): 'yes' (contrast/corr ranks) refuses the RNase P closed pair; 'core' keeps it 3/3 + qtoy2, refuses qtoy null 3/3 but OLE null only 2/3 -> cannot replace the 0.03 floor; no default changed"
metadata:
  type: project
---

2026-10-28 overnight 2x2 (contrast unset/maps x qnull off/on, 3 seeds, standalone state_split, floor off;
/scratch/cluster_scratch/atifao/qnull_port/battery4, qnull_table.py). z_q = (T_real - T_pseudo)/(sqrt2 sd null) >= 1.67.
- 'yes' (max of contrast / consensus-corr / core pseudo-splits): closed pair refused 3/3 - contrast and corr percentiles
  carry the closed state (pseudo re-test 0.20-0.23 > real 0.16-0.18).
- 'core' (corr to other half's consensus outside the focus mask, junk band), maps: closed pair keep 3/3 (z_q ~2.2),
  qtoy2 3/3, qtoy null refuse 3/3, OLE null refuse 2/3 (z_q 2.57/0.78/1.46). unset+core keeps OLE null 3/3.
- With the floor: same decisions as floor alone at the tree count; adds margin on qtoy null (rho 0.0289 vs 0.03).
- Node-local Bonferroni count (2 N_node/minpop - 1) keeps deep RNase P splits (r2111112 2/3 maps) the tree count (460 at
  692k, z_crit 3.71) refuses; nulls unchanged. User decides; no default changed.
Related: [[split-test-accepts-quality-split-on-null-2026-10-27]], [[contrast-switch-kills-real-rnasep-splits-2026-10-27]],
[[rnasep-closed-state-must-be-retained]].
- 692k route-A root (build 7, 6 states; children 5/6 read as a quality axis by stA0 corr AUC 0.28/0.25): with
  SIMPLE_SPLIT_GUARD=pairs (qnull_v5_on_v4.diff; per nearest-pair floor + core null) NO state is dropped, 3 seeds: nearest
  pairs z_q 9-13.5, rho 0.24-0.54. The old floor is all-or-nothing on the weakest STARTING pair at K>2 (a defect).
