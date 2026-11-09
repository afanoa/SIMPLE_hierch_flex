---
name: fixture-box-must-exceed-tree-box-crop
description: "TRAP 2026-11-09: synthetic fixtures for the stock hetero_tree need native box >= the tree's stage-B box_crop (88); box-64 stacks fail set_sigma2 (spectra short of the plane grid's shells). Also scale SNR0 by (160/box)^2 when binning"
metadata:
  type: project
---

The overlap-toy fixtures failed twice on 2026-11-09 (igm10-fixture session, `hetero_tree_20261028/overlap_toy/`).
- Noise: RNase P's SNR0 .0012 is per box-160 pixel. A binned stack needs .0012 x (160 / box)^2, so .0075 at box 64 and .00333
  at box 96. At .0012 on box 64, the poses were lost (98 deg median error).
- Box: the tree's stage B crops to box_crop 88 (smpd_crop 2.66 A). A native box-64 stack gives sigma2 spectra only to shell
  32, short of the 44 needed, so stage B fails in set_sigma2. Use a native box >= 88; box 96 / 2.44 A keeps qtoy's 234 A field.

**Why:** I briefed "box 36-64" for a cheap toy, and both mistakes cost a rebuild.
**How to apply:** a fixture for the stock tree uses a native box >= 88, scales the noise to the box, and gates its trees on the
refined pose error (< 15 deg) before launching. Related: [[use-the-full-cpu-cap-for-tests]].
