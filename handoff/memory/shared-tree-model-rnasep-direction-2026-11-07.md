---
name: shared-tree-model-rnasep-direction-2026-11-07
description: "USER DIRECTION 2026-11-07 night: pursue the shared tree model (one sparse site basis, per-site forest, fission-honest tests) until it works on RNase P; working on IgM-like regimes is already 'a great result'. The RNase P blocker is per-particle misfit (amplitude / fit quality), which is being modelled"
metadata:
  node_type: memory
  type: project
  originSessionId: cc23038c-fea9-4b5f-b321-d58eb2942350
  modified: 2026-11-08T04:28:26.677Z
---

User, after the shared latent tree failed RNase P: "We need to pursue this to work on the RNase P. If it works on IgM like
regimes that is a great result." Earlier the same evening: "our current approach lacks elegance which is why I am interested
in a shared model (of course refinements per subset are necessary)", and "port the shared tree latent if it's useful yet".

**Why:** the shared model is the user's preferred design for the method and the paper. It already beats the production tree on
the IgM-like fixtures and holds every null, but on RNase P its splits track per-particle misfit rather than the closed state
([[tree-latent-shared-basis-verdict-2026-11-07]]).
**How to apply:**
- Keep the shared model as the main methods thread.
- The RNase P route runs through the per-particle observation model: the rnasep-misfit session characterises the misfit and
  tests per-particle band scale / noise variance / interpolation-aware noise against the solvent-control floor. The tree-latent
  session runs the misfit-removed rerun (amplitude + nuisances profiled, stratified sat5).
- Port to Fortran once RNase P shows the closed state honestly (r131 corner or the 692k closed lineage), using the forest design in
  `handoff/tree-latent/PORT_DESIGN.md`.
- Present IgM-like success as a result in its own right.
- 2026-11-08: the noise model is found and works ([[rnasep-misfit-noise-spectrum-verdict-2026-11-08]]). With it, the 692k root
  split stops tracking noise but splits along amplitude (a 10 % low-amplitude side), with closed not at level 1. Deeper growth,
  honest maps of both sides, the r131 tree and the native port were running that night.
- 2026-11-08 ~05:30: at r131 neither the latent forest nor the site chain's full-image round isolates closed (7 % at d' ~.9; the
  seeds halve the node, and the halves' maps read the node mean). The production tree finds the corner only at depth 2, after a
  halving carve that enriches closed 1.5x. Open design question: allow an exploratory carve, then an honest test one level down.
Related: [[site-coupling-factorial-verdict-2026-11-07]], [[sparse-ppca-pursuit-verdict-2026-11-07]],
[[rnasep-closed-state-must-be-retained]].
