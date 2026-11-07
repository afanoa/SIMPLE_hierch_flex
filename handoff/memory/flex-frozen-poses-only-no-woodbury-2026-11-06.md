---
name: flex-frozen-poses-only-no-woodbury-2026-11-06
description: "STANDING HARD GATE (user 2026-11-06): the latent model (flex) works only on FROZEN poses; no recalculation of the latent model inside refinement and no Woodbury / SOLVAR-type marginal algebra in refinement, unless rigorously proven useful"
metadata:
  node_type: memory
  type: feedback
  originSessionId: cc23038c-fea9-4b5f-b321-d58eb2942350
  modified: 2026-11-06T19:16:15.565Z
---

User, 2026-11-06, on my flex + posterior-machinery proposal (doc/implementation_notes/planned/flex_posterior_integration.md),
whose layer L2 put flex's z-marginal Woodbury term into every refined row entry and the latent motion into the insertion:
"I do not want to recalculate or do any woodbury algebra like SOLVAR. A hard gate (unless rigorously proven useful) is that
the model works only on frozen poses."

**Why:** every marginal-likelihood pose test of ours hurt or did little (flex_pca shift polish 2.55 vs 0.56 px, flex2D
9.4 -> 15.9 px, toyJ ~15 %), and the user wants flex and refinement coupled without flex inside the pose search.
**How to apply:** integration proposals keep flex outside refinement. Flex fits at the poses refinement committed, and
refinement only reads what flex publishes (the state weight set, labels, seeds, masks), at stage or node boundaries. No
latent term in the row, the in-plane, shift or state draws, or the insertion. Lifting the gate needs a rigorous
demonstration first. Related: [[solvar-marginalizes-z-during-fitting]], [[joint-pose-heterogeneity-toy-verdict]],
[[tree-as-full-mppca-assessment-2026-11-06]].
