---
name: refine3d-auto-ref-pose-init-only-for-foreign-vol1
description: "TRAP: ref_pose_init=cc re-searches every pose from scratch; never pass it after solve3D on the same particles"
metadata:
  node_type: memory
  type: feedback
  originSessionId: 851a7b77-9a45-48b1-b2b6-396c1a7b5829
  modified: 2026-11-01T20:12:45.619Z
---

`refine3D_auto ref_pose_init=cc` runs a fixed-reference CC pose search at 15 Å that discards the project's poses. It is only for a vol1 of independent provenance. After solve3D on the same particles, run refine3D_auto on the solve3D project with no ref_pose_init and no vol1 (it picks up the project volume and keeps the poses).

**Why:** on TMEM16A (2026-11-01) it collapsed cFAR from 0.44 to 0.027; the user called it "stupid". I had copied the flag from earlier NOT/RNase P launchers without checking what it does.
**How to apply:** check every refine3D_auto key against its UI description before launch ([[check-all-program-keys-before-launch]]); never carry a recipe's flags to a new chain blindly.
