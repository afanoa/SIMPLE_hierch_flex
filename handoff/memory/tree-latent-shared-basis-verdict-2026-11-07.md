---
name: tree-latent-shared-basis-verdict-2026-11-07
description: "VERDICT 2026-11-07: a tree grown in latent space over ONE tree-wide sparse site basis (soft Bishop-Tipping hierarchy with XD, growth by the fission saturation test, frozen poses) passes every null incl. the IgM 10 % twin, isolates qtoyflex closed (F1 .558 vs .39/.47) with the arm whole, IgM 10 % five-arm .26 (tree .042, lever 3 .338), qtoy2 .76 (pp_on .811), IgM 5 % one leaf"
metadata:
  node_type: memory
  type: project
  originSessionId: cc23038c-fea9-4b5f-b321-d58eb2942350
  modified: 2026-11-08T02:03:49.454Z
---

tree-latent session (user: "a cohesive model over the entire tree"); `hetero_tree_20261028/handoff/tree-latent/{VERDICT,STATUS}.md`;
runs `math_reworks/tree_latent/`. Python on sparse PPCA's bases; nothing in SIMPLE.
- Model: the sparse site basis B shared by the whole tree and deconvolved node Gaussians (XD with the posterior precision A,
  soft responsibilities). Growth: a site's K = 2 XD proposal is accepted by sat5 (fission saturation).
  - v3 = proposals and sat5 input labels on y, node memberships and child weights on fission copy y1. Membership is a selection
    too, so it must be fission-protected.
- Results (labels on y):
  - qtoyflex: 3 leaves (arm whole, holo, closed). Closed is its own depth-2 leaf (F1 .558, purity .63, recall .50), holo .713,
    kind NMI .399. Gate trees: closed .39 / .47, holo .62, kind .255 / .265.
  - qtoy2: 8-9 leaves, NMI .752-.764, closed .84-.90, holo .95-.96. It loses the conformer 1 | 2 pair (no K = 2 proposal
    separates it). pp_on .811.
  - IgM 10 % fix: five-arm product, joint .260 / forest .258; the 2-state arms reach lever-3 quality (a1 .394, a3 .372, a4
    .340).
  - IgM 5 % fix: 1 leaf (root z_mid 2.2-2.8; half the information is too little at d' ~1).
  - Nulls: 1 leaf everywhere (IgM 10 % twin 0 / 5 DISCRETE; the hetero tree splits it).
- Forest vs joint: equal. The forest never splits a 3-state arm, because the balanced K = 2 split (state 1 vs 2 + 3) straddles
  and sat5 calls it continuous; lopsided splits that isolate one state pass.
- No null accept in any variant, so no multiplicity correction was needed.

**Why:** the user wants one elegant model over the whole tree; this is its first prototype.
**How to apply:** the shared-basis latent tree is the candidate cohesive model.
- Next (approved): K = 3 proposals with a pairwise saturation test in forest form, for 3-state arms and close pairs.
- Its limits: power at d' ~1 (fission halves information), and frozen poses. The real tree must still refine poses per subset
  (block-coordinate ascent).
Related: [[sparse-ppca-pursuit-verdict-2026-11-07]], [[saturation-fission-site-test-verdict-2026-11-07]],
[[per-site-discreteness-test-gapw-verdict-2026-11-07]], [[lever3-latent-refs-seeding-verdict-2026-11-07]],
[[tree-as-full-mppca-assessment-2026-11-06]].
