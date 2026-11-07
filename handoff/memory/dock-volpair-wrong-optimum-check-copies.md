---
name: dock-volpair-wrong-optimum-check-copies
description: "TRAP (2026-11-06): dock_volpair can land in a wrong optimum (s6 full/even onto P_asis cc .72, 60-90 deg off; odd .878 right). Dock every copy and cross-check, or anchor on the best and dock the rest onto it (.998)"
metadata:
  node_type: memory
  type: project
  originSessionId: 371593ba-f05b-45f2-801d-d68e405eef0b
  modified: 2026-11-06T13:29:10.668Z
---

Found 2026-11-06 (IgM DJ session, `ph_tests/s6score*`).
- What happened: `dock_volpair lp=20 mskdiam=280 trs=15` docked DJ r300k state 6 onto P_asis (EMD-13921 as-is in the consensus frame).
  - Full map and even half: cc .722 and .722, in one orientation.
  - Odd half: cc .878 (right: it equals s6 vs EMD in EMD's own frame, .877).
  - The first two sat 60-90 deg from the odd half (corr .09). Same reconstruction, deterministic search, two optima.
- Fix that worked: anchor on the copy with the highest cc, then dock the other copies onto the anchor (near-identical maps: cc .998, halves vs full .97 at 40-16 A).
  - The anchored placement reproduced the EMD-frame register test exactly (core CC k=0..5 .740 .578 .542 .606 .597 .501 vs regemd .740 .571 .539 .602 .597 .490).
- Also failed: a python own-axis search (max 60-deg self-corr, 30-8 A) on a map without mid-resolution signal (s6 self-corr .10, axis 6 deg / 9 px off).
  - The EMD control placed that way matched P_asis only .716.

**Why:** a wrong-optimum docking silently mis-places every downstream reference (masked scorer references, J-slot readouts); the first s6 scorer run read J at 114 deg instead of ~197.
**How to apply:** dock all copies of a map (full, even, odd) and require mutual corr >= .95 in the target frame before using any of them; or chain a validated docking (map -> EMD frame as in regemd, then the EMD -> target transform). Related: [[igm-11077-consensus-poses-tree-ready-2026-11-04]], [[refine3d-lp-set-scores-merged-map-trap]].
