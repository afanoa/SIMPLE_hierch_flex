---
name: relion-small-state-local-refine-blurs
description: "RELION per-state local refine of 1.5-9k-particle tree states is WORSE than reconstructing from the SIMPLE consensus poses (OLE: +0.1 to +2.5 A, every state); cause = SSNR<1 current-resolution criterion, not the poses; best variant --tau2_fudge 4 still loses"
metadata:
  node_type: memory
  type: reference
  originSessionId: a25e6f6d-c94c-4ba6-92f7-c1fb7553af83
  modified: 2026-11-05T06:25:11.389Z
---

2026-11-04/05, OLE 12707, 31 tree leaves, RELION 4.0.1, C2, box 360 at 0.954 A. The RNase P leaf_direct recipe (own sum low-passed at
12 A, ini_high 12, healpix 4, sigma_ang 3) stalled 1.7-3k-particle states at 6.5-10 A, although the same particles from the SIMPLE poses
give 4.1-4.9 A (relion_reconstruct + postprocess, same mask). RELION moved the poses only ~0.6-1 deg (C2-aware: a C2 mate is rot+180,
compare A vs B.C2 with C2 on the RIGHT). Its iteration-1 gold-standard FSC still crossed 0.143 at ~4.5 A, but the per-shell SSNR stayed
0.2-0.7 from 10 to 5 A, so auto-refine's current resolution (SSNR = 1) sat at ~10 A. It then aligned against an over-regularised
reference, Pmax fell 0.69 -> 0.02 and the maps blurred. A full-resolution consensus reference (ini_high 3.5) stalled the same way.
--tau2_fudge 4 (ini_high 6, healpix 5, sigma_ang 1.5) unstuck it, but every state still ended worse than the consensus poses (r11111
9.1k: 3.43 -> 3.54; r22x states: +1.7 to +2.5 A). RNase P showed the same effect, milder (10-20k particles).

**How to apply:** for per-state maps of tree leaves, reconstruct from the consensus poses after one CTF refinement of the consensus
(relion_reconstruct halves + postprocess); do not local-refine states under ~10k particles unless asked, and if asked, use
--tau2_fudge 4 and report the consensus-pose resolution beside it. Related: [[rnasep-relion-polish-pilot-2026-11-03]],
[[relion-reconstruct-maps-need-fsc-weighting]], [[ole-relion-state-maps-2026-11-05]].
