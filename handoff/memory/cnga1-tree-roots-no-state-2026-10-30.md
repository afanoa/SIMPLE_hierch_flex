---
name: cnga1-tree-roots-no-state-2026-10-30
description: "CNGA1 root trees (C1/C4 x no mask/support_c3, 245k): all split, none a state. C1 = micelle halo / view imbalance; C4 = a 25-35 % worse-particle group (dB +160-190); low-res pore-tilted residual not gating-shaped. Roots cost 85-125 core-h each (408 total, over budget)"
metadata:
  type: project
---

2026-10-30, hetero_tree_20261028/cnga1_tree (21fac34db + legacy switches, HALVES/HALFSEED off, flex without pgrp).
- c1_nomask K*4, z 13.1: belt 3.3x enriched + corr. c1_support K*4, z 4.97: tied to view tilt (AUC 0.255).
- c4_nomask K*2, z 24.8; c4_support K*2, z 13.9: the minority group is low-corr, pose-unstable and blurrier (dB +160 / +193).
  - B-matching removes the 15-8 A difference.
  - A 40-15 A residual stays, tilted to the pore (2.1-2.2x), but only 24-27 % of it lies within r<20 A (open-closed gating: 71 %).
- Verdict: no state, consistent with [[cnga1-is-the-apo-closed-dataset]] and [[cnga1-continuum-is-the-micelle-halo]].

**Why:** the user asked for CNGA1 root trees with C1/C4 arms and a protein support. Maps go to the user only for a split that survives nuisance checks.
**How to apply:**
- Price CNGA1 root trees at ~85-125 core-h each, not the 45-75 scaled from RNase P, and track spend in flight.
- On C4 data read the split groups against corr, pose change and dB before any state claim; the B-match test plus the cylindrical profile vs the gating reference separate blur from structure.
Inputs need [[cnga1-legacy-project-eo-and-indstk-traps]].
