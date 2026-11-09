---
name: pose-weighted-flex-dropped-pose-spread-2026-11-09
description: pose-weighted flex latent dropped: the wip draw-to-mode pose spread at r1 (direction p50 2.9 deg, shift .2 A) moves a 15 A projection by a tenth of a 1-sd latent step (ratio p50 .13-.17; 15-16 % of particles >= .5 vs the 20 % gate) and carries no latent or closed information (|r| <= .04)
metadata:
  type: project
---

node_levers/nl_posespread.py, handoff/node-levers/POSE_SPREAD.md (2026-11-09). Six arm-b r1 replays, drawn vs settled poses:
direction angle p50 2.88 deg (42 % unchanged, p90 4.9, p99 7.7), in-plane .45 deg (p99 36), shift .2 A (p99 6.9 A). Fourier-slice
projections of the r1 mean at both poses, 15 A low-pass, node-mask weighted: 1 - r p50 .0038 vs .019-.024 for a 1-sd step along
z1 (mode volumes regressed from the eight flex state maps); ratio p50 .13-.17, share >= .5 is 15-16 %, >= 1 is 5-6 %; the effect
follows the shift (r .92), not the direction (.27). Leakage of every spread measure into z1..z4 and the closed rows |r| <= .04.

**Why:** the coordinator gated a pose-weighted flex (frozen w_nk per candidate pose) on the pose spread being wide at flex's band;
it is not.

**How to apply:** do not prototype pose-weighted flex at r131 from this spread; flex on the settled (or drawn) single pose is
adequate at 15 A. Related: [[wip-settling-pass-prob-beta-cold]], [[wip-refine3d-redraws-poses-each-iteration-2026-11-08]].
