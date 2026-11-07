---
name: test-runs-moved-to-delete-folder-2026-10-28
description: "2026-10-28: the Oct 27 session run folders (~14 TB) were renamed into /mnt/beegfs/elmlund/afan/_to_delete_20261028 for the user to delete; scripts backed up on scratch; new runs live on /scratch/cluster_scratch/atifao"
metadata:
  node_type: memory
  type: project
  originSessionId: 9643ce3d-fc7a-4dc1-860a-59db10867e7c
  modified: 2026-10-28T04:35:45.964Z
---

On the user's request ("move these runs to one folder and I will delete"), renamed (os.rename, no copy) into
`/mnt/beegfs/elmlund/afan/_to_delete_20261028/`: hetero_tree_regression_20261026/{amp_mstep_v1, start_point_v1, recipe1},
ole_12707/split_merge_replay, ole_12707/18_tree_recipe1_s1_ioerr_2333, rnasep_klem/site_states, swing_switch_20261027,
tree_port_optim_20261027, tree_data_lp_20261027, pool2_test (about 14 TB). KEPT in place: the regression suite in
hetero_tree_regression_20261026 (fixtures, suites, scripts, runs/ incl. runs/real_root_split_contrast_m8c52_s1/closed_reference_rows.txt,
null_controls/ole_null_v1 = the OLE null fixture), ole_12707/6_hetero_tree (Oct 24 reference, 1.8 TB). Small scripts
(*.sh *.py *.sbatch *.tsv *.md *.json *.diff < 5 MB) of the moved folders: /scratch/cluster_scratch/atifao/tools_backup_20261028/<full path>.

**Why:** beegfs Default pool 99.9 % full; our two days of uncropped test trees wrote ~14.7 TB.
**How to apply:** old paths under those names are gone; write new runs to /scratch/cluster_scratch/atifao (recipe4/recipe5,
site_states/full692 there). Builds: recipe5 = fixm source + ITA=ITB=3 + validated crops + embedding fix + map cleanup
([[beegfs-transient-write-failure-2026-10-27]], [[integration-plan-2026-10-27]]).
