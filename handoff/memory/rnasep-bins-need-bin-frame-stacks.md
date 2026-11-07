---
name: rnasep-bins-need-bin-frame-stacks
description: RNase P direction bins: the disc bin averages of 2026-10-15 under every sign convention were the SCRAMBLED IMPORT (see rnasep-692k-project-rows-are-a-permutation-of-stack-order), not a convention problem; the transform rotate(-e3) shift(-(x,y)) verified on raw members (read by keys) is applied to the images (bin-frame stacks, e3=x=y=0); cls_expansion overwrites ptcl2D class with sub-class ids; import does not carry state
metadata:
  type: feedback
---

2026-10-15, direction bins on RNase P (notes 21.1-21.3). The consolidated stacks are pre-aligned to 2D-class frames;
to group particles by projection direction one must bring each image into its bin's projection frame. make_cavgs with
ptcl2D (e3, x, y) copied from ptcl3D, negated, half-negated or zeroed all gave featureless discs (frame_score.py
correlations 0.56-0.59 at arbitrary rotations); on raw members the transform rotate(-e3) then shift(-(x, y)) (scipy
conventions, box 160; shifts halved at box 80) gives sharp averages matching the consensus reprojection
(frame_analysis.py). Second trap: cls_expansion writes its sub-class ids into ptcl2D class, so any rerun on the same
project splits the sub-classes again (22k sub-classes) and any orientation table for the dumps must be exported
BEFORE cls_expansion. Third: import_particles does not carry state; set state 0 for never-refined rows with
selection on ptcl2D and ptcl3D before make_cavgs / cls_expansion.

**Why (corrected 04:40):** every make_cavgs/cls_expansion test on an imported bins project paired poses with the
wrong images (import re-keys rows sequentially; the 692k project's keys are a permutation of stack order), so the
"no convention works" finding is void; the raw-member tests read images by keys and stand. Whether the ptcl3D
in-plane parameters copied into ptcl2D would work directly on a row-ordered import is untested; the bin-frame
stacks (v2, row order) are the route in use.

**How to apply:** build bin-frame stacks (run_bins_transform.sbatch, box 80), import them with e3 = x = y = 0 and
class = bin, run the 2D side there, reconstruct with the original project by particle index. Related:
[[rnasep-stacks-are-class-prealigned]], [[rnasep-bins-overnight-2026-10-15]],
[[rnasep-692k-project-rows-are-a-permutation-of-stack-order]].
