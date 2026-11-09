---
name: discrete-splits-capture-the-span-states-that-matter
description: "STANDING (user 2026-10-23): the state tree's discrete splits are the way to go; the goal is the conformations that MATTER across the entire span (extremes included), not a map at every point of a continuum; don't push continuous-landscape output as the main deliverable"
metadata:
  node_type: memory
  type: feedback
  originSessionId: 25d47be4-7d1d-48d8-abdf-97b1c3c946ab
  modified: 2026-10-24T23:25:52.529Z
---

User, 2026-10-23, after I suggested sampling the landscape continuously along flex's latent directions instead of "forcing
discrete leaves": "well the discrete split is the way to go, we just have to be smart about it so that we can capture the entire
span. The goal is not to produce every conformation at every point, but the ones that matter. And it does work on the RNAseP."

**Why:** the deliverable is a set of distinct, interpretable states that cover the full conformational range (e.g. RNase P's closed
end through the most extended open conformer), which the tree already achieves on RNase P.

**How to apply:**
- Improve how the tree proposes and keeps discrete states: latent-first clustering, the split test, the per-site layer.
- Judge it by span: are the extremes reached? And are the states that matter recovered (purity / recall against truth on toys and
  CAHRA)?
- Judge reproducibility by whether the same key states reappear across seeds (map matching), not only by label ARI.
- Keep landscape_sample / continuous profiles as a diagnostic, not the main output.
- 2026-10-24, user, after the 692k tree split nearly every node to depth 8: "with more splits you lose resolution, thats why I am
  interested in producing extremes of a certain motion/composition". So deep splitting is COUNTERPRODUCTIVE: each extra split
  divides a state's particles among its data-derived groups (K* chosen by the split test, uneven sizes; in practice mostly 2,
  because a node is usually seeded with two groups) and costs resolution. Aim for the extremes of each motion / composition axis with enough
  particles to resolve them. Prefer caps, stricter split tests and merging over finer slicing; on that request the 692k tree was
  capped at maxdepth 7 and finished.
- 2026-10-24 about 19:20, the user restated the GOAL; this is the current target:
  - The tree plus a particle floor already gives good high-resolution states, organised meaningfully. Keep it.
  - The open questions are QUANTITATIVE: how much a region swings; whether the swing is concerted (a two-state switch) or a
    continuous motion; and whether particles are spread non-uniformly along it, which reveals the dynamics.
  - At EACH SPLIT, tell discrete change (split) from continuous swing (don't slice it; describe it with poles = the latent
    extremes, plus the occupancy profile along it). The user called this "an extension of the pole idea".
  - Post-hoc binning after the tree is "ad hoc". So is reconstructing every point of a swing. Neither is wanted.
  - The user invited graph theory here, e.g. spectral graphs / diffusion maps, where a bottleneck = discrete.
  - The user added: "the goal is to still get states which I can refine afterwards". So the output stays refinable particle sets:
    - a discrete split gives states, as now;
    - a continuous swing gives two POLE STATES (the most extreme particles at each end of the swing coordinate, each at least
      the leaf particle floor), with the middle kept as one state;
    - the swing's amplitude comes from the refined pole maps, its dynamics from the occupancy profile.
  - Slicing source in the code: flex's continuum macro-cluster is cut into k-means regions (simple_flex_pca_gmm.f90
    continuum_region_weights), and the tree seeds its splits from those regions.
- Related: [[landscape-must-be-reproducible-and-splits-data-derived]], [[klem-priority-every-state-before-refinement]].
- 2026-10-27, user on OLE 6_hetero_tree's 1,456-1,986-particle leaves: "those ... states are really good and seem to sample the
  conformational landscape really well. The priority is not the number of states but a set of states which represents the
  landscape. This could be compositional or conformational." So small leaves are WANTED; do not impose the research note's
  "OLE >= 2,000 per state" floor or skip sub-2*minpop nodes in a way that drops them. Judge a split/merge rule by how well
  its state set covers the landscape (conformational and compositional), not by K.
- 2026-11-08 user, on the robustness redesign: "remember we are trying to produce the most diverse set of conformations possible".
  - Robustness fixes must ADD states at uncertain decisions (joint seeding over tied sites, keeping both candidates), never trade
    diversity for a single safe choice.
  - Judge a fix on two things: the closed-state hit rate across seeds, and the number of distinct, honestly different
    conformations (span and extremes).
