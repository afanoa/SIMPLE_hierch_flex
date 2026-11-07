---
name: fixed-pose-node-mask-em-is-not-refine3d-2026-10-22
description: "TRAP (fast-tree switch 1, r1): a K-class site-EM-style classifier at fixed poses over the broad NODE mask does not reproduce refine3D's states (stA0 / stage A+B) - soft+priors collapses, hard CEM jitters; discovery then picks the wrong site and r1 becomes one leaf; works only inside the local SITE mask"
metadata:
  node_type: memory
  type: project
  originSessionId: a4258345-b249-48ff-a0d3-d33e99ae33ac
  modified: 2026-10-22T21:17:50.152Z
---

2026-10-22, fast-tree session (worktree interesting-diffie-091064, batch2 ~/SIMPLE_fast, runs rnasep_klem/fast_tree/{fp1,fp2,emvar}):
`prg=state_em` = site_rank's 2-class EM generalised to K classes (nn maps, consensus + FSC-weighted masked deviation, other-half refs,
per-shell scale/noise E-step), at the node's input poses inside the node falloff mask, as a replacement for stA0 and stage A/B.
- From the 8 fflex labels: soft with class priors collapses (one class 17% -> 85%, never converges, 10% still changing at 50 its);
  flat priors or hard CEM keep 8 classes but 1-4% of labels flip every iteration. NMI vs refine3D stA0 from the same labels 0.03-0.07.
- Discovery on those state maps picks 1f2 (repro 0.68 soft / 0.94 hard; site_01 Dice 0.72-0.79 to the 1f1 closed-arm site), then
  site_rank accepts no pole and r1 becomes ONE leaf (current chain: P = 2, closed 0.90/0.55 vs open 0.53/0.93).
- From the build-47 seeds: NMI vs stage B 0.03 soft / 0.10 hard; even at stage B's own refined poses 0.19 / 0.24 -> a different
  classifier, not only a pose effect. Inside the SITE mask the soft K = 2 EM reproduces the site pole (r = 0.82 with the endpoint
  posteriors); hard there matches stage B's sizes (40k vs 38k closed) but NMI 0.27.
- The current chain reproduces itself: v47a vs v47b site Dice 0.997, seeds NMI 0.98, stage-B NMI 0.95.

**Why:** the user's speed plan assumed fixed-pose classification gives the same states with blurrier maps (true for the 1-D site
ranking in a local mask, see [[minority-soft-em-fixed-poses-no-rescue-2026-10-21]]); in the broad node mask it does not.
**How to apply:** do not replace refine3D at interior nodes with a fixed-pose node-mask EM; cut refine3D cost instead (fewer stA0
iterations, skip stage A, real CPUs for the site step) and validate against v47a/v47b. Related: [[poles-are-minority-vs-rest-k-includes-rest-2026-10-22]].
