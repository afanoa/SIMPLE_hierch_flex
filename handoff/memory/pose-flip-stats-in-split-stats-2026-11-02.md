---
name: pose-flip-stats-in-split-stats-2026-11-02
description: state_split writes report-only pose-change lines to split_stats.txt ('# pose node flip30 ...', '# pose split|states pair g|h ... dflip d z z ... flag'); port of the coordinator's pose-artefact study; on the gates the flag (dflip >= 0.04, z >= 5) fires only on micnull's refused root pair (0.071, z 20.6)
metadata:
  type: project
---
2026-11-02, session tree-lever-floor (claude/tree-lever-floor, note hetero_tree_reworks_2026_11.md section 12; state_mergeup
pose_lines, exec_state_split reads ../n.simple): per particle the geodesic angle input pose -> stage-B pose over the point group
(euler2m, sym%apply on the right), flip = > 30 deg; pairs = delivered children after guard + floor (stage-B states where none);
z from 200 label shuffles inside the split test's own strata (stage-B view bins x defocus tertiles, not the study's input-pose bins).
Gate nodes: flip30 qtoy2 .0009, r131 .001-.002, nulls .002 / .00006, OLE root .021, micnull .273; pair dflip <= .0014 on needed
splits (z <= 5), OLE root .008-.021 (z to 14.7), micnull stage-B pair .071 z 20.6 = the only flag. Nothing gates on it.
**How to apply:** read the lines from split/split_stats.txt (comment lines after the numeric line); a 692k tree's r4 (0.106) is the
case the study expects to flag. Related: [[effect-floor-pn2-after-guard-2026-11-02]].
