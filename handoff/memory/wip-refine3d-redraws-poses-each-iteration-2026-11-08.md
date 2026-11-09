---
name: wip-refine3d-redraws-poses-each-iteration-2026-11-08
description: wip tip (918215446) refine3D samples poses, no fixed point (72 % of directions re-drawn per iteration); r17d build settles; a cross-build evidence pass is not that refinement's posterior
metadata:
  type: project
---

Measured 2026-11-08 at r1 (252k particles, box 88, lp 10, nspace 5000): one extra stage-B iteration from a run's final project
moves, with the r17d (production) binary on its own run, 5 % of directions by > 1 deg and 2.4 % of states (a fixed point); with
the wip binary on its OWN run 72 % of directions by one grid step (p50 2.92, p90 5.8 deg) and 4.6 % of states; wip on r17d 62 % /
4.6 %; r17d on wip 62 % / 13 %. The wip's topk / keyed-draw / polar-shift commits (since 2026-10-29) sample poses from the
posterior instead of taking the argmax. Consequence: the r17d binary's prob_state_evidence.txt computed at wip poses is flat
(median gap .3 nats vs 2.0 on its own poses) even though both builds' state maps are interchangeable (swap test).

**Why:** a "responsibilities" readout of a wip run with the production binary looked uniform and would have been reported as
"the states are indistinguishable"; the swap and the one-iteration movement tests showed it was the pose convention of the build.

**How to apply:** per-particle state posteriors must come from the same binary that refined the poses; the wip tip needs
write_state_evidence ported (code change) before any responsibility readout of its runs. Refinement churn (2-5 % of states per
late iteration) is not the source of cross-run label disagreement (see [[node-levers-replay-verdict-2026-11-08]]).
- POSE-WEIGHTED FLEX LATENT, DROPPED (user idea 2026-11-09; f2r measured drawn vs settled poses at r1): spread p50 2.9 deg
  direction, .2 A shift. At flex's 15 A band the pose effect is a tenth of a 1-sd latent step (ratio >= .5 for 12-16 %, bar
  20 %), driven by sub-pixel shift. It correlates with z1-z4 and closed within |r| .04.
  - REOPENED by the user ("flex band should be data derived ... 2.9 degrees is quite a bit"): the proxy fixed 15 A, and the
    ladder's 15 A pick was made at single poses. Direct test at r131 (8 posterior draws, information-averaged latent vs
    settled pose, ladder band per arm) briefed to rnasep-misfit.
  - Stage 1 at r131: draws median 2.92 deg from the TEMPERED posterior (beta .3, N_eff 5.3 vs 2.3 at beta 1). Early stage 2:
    the equal-weight 8-draw latent is far weaker than single pose. TRAP: pose weights must be the beta-1 posterior,
    not the search temperature's draws. A draw-1 control decides between redrawing at beta 1 and finishing.
