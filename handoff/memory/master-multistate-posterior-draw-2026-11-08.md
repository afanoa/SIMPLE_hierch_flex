---
name: master-multistate-posterior-draw-2026-11-08
description: "origin/master cac010e40 (Hans, 2026-11-08): multi-state state = one draw from the tempered joint (state, direction) posterior, final assignment at prob_beta 1, NO population prior; its 3-state benchmark never isolates the 10 % state and half the labels flip per iteration at 20-23 A; don't merge into the tree's stage A/B ungated"
metadata:
  type: reference
---

The user asked me to look at master's multi-state work, suspecting that "multi state is scrambling membership too much".
Master has 9 commits since our merge base 6634307d9; only cac010e40 concerns states.
- Design:
  - assign_rows draws state and direction jointly, p ~ exp(-T (L - L_min)), with labels hard. The neighbourhood table pools
    every state's candidates.
  - The fill budget and the annealing stall are per state, and multi-state runs end with a final assignment at prob_beta = 1.
  - nstates = 1 is bit for bit unchanged, and no key was added.
  - There is no class-size term.
- Report (completed/multistate_posterior_state_draw_report.md), on 20k simulated ribosome assemblies at 30 / 60 / 10 %:
  - The 10 % state never gets a state of its own; 60 % splits over two. De novo accuracy is .42-.59, and .65-.69 from random
    labels at true poses.
  - At 22.6 / 20.2 A about half the labels change every iteration, with the posterior spread over 2.2-2.7 of 3 states.
  - Labels organise only at 20 A and finer.
- Planned there, not built: a learned "pose student" (single-state) that proposes directions for refine3D.

**Why:** it is the same failure as ours (no class prior, membership scrambled at low resolution, minorities absorbed). It backs
the user's hypothesis and the fixed-membership lever (guide 7.10 thread 12).
**How to apply:**
- Never merge cac010e40 into the tree's stage A/B without the gates and the r131 honest read. A tree-wide prob_beta = 1 is the
  one arm that lost the corner ([[tree-master-merge-gates-2026-11-06]]).
- When proposing state rounds, use a class prior and keep membership fixed while poses refine.
Related: [[shared-tree-model-rnasep-direction-2026-11-07]], [[extremes-fission-verdict-2026-11-07]].
