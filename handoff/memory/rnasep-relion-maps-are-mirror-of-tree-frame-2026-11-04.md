---
name: rnasep-relion-maps-are-mirror-of-tree-frame-2026-11-04
description: "RNase P RELION per-state maps (raw SIMPLE poses) are the MIRROR of the tree maps' frame, turned 165 deg; one docked transform (rel_xform.json) registers all of them; RELION hand = deposited hand"
metadata:
  node_type: memory
  type: project
  originSessionId: dc17548d-1b94-4fe8-ad02-93f812b3f2a7
  modified: 2026-11-04T22:33:38.451Z
---

2026-11-04. RELION refines from the raw-frame SIMPLE poses (Refine3D/leafd_<L>, leaf_direct.sbatch; the user's preferred
workflow: "all the volumes in the slides to come from this workflow") land in a different frame from the f692_lf_s1 tree maps:
mirror + 164.6 deg rotation + <1 voxel shift. Docked once (union consensus Refine3D/job018 onto the tree root map, both hands,
20k rotations + local refine, cc 0.974 at 8 A); rnasep_slides_20261103/rel_xform.json, applied by rel_apply.py ->
maps/relx/<leaf>.mrc (+ relx.json levels); every state then matches its tree map at 0.95-0.98 (8 A).
The tree frame itself is the mirror of the deposited models (models_frame.npz hand -1), so the RELION maps have the
DEPOSITED (true) hand; renders made in the deck frame are mirrored left-right to show the true hand.
The 51 non-union states run as leaf_direct_r.sbatch (array 54948410, csbdevel); outputs Refine3D/leafd_<L> too.

**How to apply:** never compare RELION and tree maps without rel_xform.json; rerun rel_apply.py (idempotent) as new
leafd maps land. Related: [[rnasep-relion-polish-pilot-2026-11-03]], [[deliver-maps-unsharpened-for-judging]].
