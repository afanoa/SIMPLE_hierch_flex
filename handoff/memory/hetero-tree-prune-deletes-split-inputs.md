---
name: hetero-tree-prune-deletes-split-inputs
description: "TRAP: a finished hetero_tree deletes split/s.simple, site masks and every stage-B map unless SIMPLE_TREE_KEEPMAPS=yes; replays then lose the locality guard"
metadata:
  node_type: memory
  type: project
  originSessionId: f9f21a10-8a4f-494f-920d-b5692d106e48
  modified: 2026-10-30T00:41:08.162Z
---

hetero_tree prunes on S_DONE (prune_node) and at the end (prune_maps): split/s.simple, sites/*.mrc and all node .mrc not used by a leaf go.
stageB/b.simple survives (s.simple is its copy), so a state_split replay can still run, but the locality guard reads ../stageB half maps and
logs "half maps not found; skipped" -> the replay's K* is pre-guard. site_rank replays need the site masks.
**Why:** hit 2026-10-29 replaying L10_12k's root split (detection-limit study); had to rerun the node with maps kept.
**How to apply:** export SIMPLE_TREE_KEEPMAPS=yes on any run whose decisions may be replayed (5-14 GB per 60k run, mind the ~50 GB folder norm);
detlim/scripts/launch_battery.sh does. See [[rare-state-detection-limit-2026-10-29]].
