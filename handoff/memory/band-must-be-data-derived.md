---
name: band-must-be-data-derived
description: "STANDING (user 2026-11-08): the resolution band (lp cap) of site / latent fits must be data-derived, not a fixed constant (not LPF 15, not 20); choose it by held-out gain per node or per block (NU-filter style)"
metadata:
  node_type: memory
  type: feedback
  originSessionId: cc23038c-fea9-4b5f-b321-d58eb2942350
  modified: 2026-11-08T05:12:40.122Z
---

User, after sparse PPCA found all IgM arms at lp 20 but none at the tree's LPF 15 (seed geometry ruled out): "really the band
should be data derived".

**Why:** a fixed band is a data-specific knob. IgM's arm signal lives below ~20 A while qtoy and RNase P work at 15 A, so any one
constant fails some molecule. This is the same principle as the knob-free rule ([[flex-pca-must-be-knob-free]]).
**How to apply:** for sparse sites (and any site or latent fit), pick the band cap from the data:
- the cap that maximises the leading seeds' held-out gain per node; or
- per block, each block's modes taking the band with the best held-out gain (the nonuniform-filter analogy).
Use a fixed geometric ladder of caps and 2 SE as the only constants. Never propose "use lp 20 for IgM" as the fix.
- 2026-11-08, user, on the ladder's 1.6-4.7x cost: "couldnt you bake this into the latent variable model itself ... we know how
  particles are related at this particular frequency and lets stick with that" (like sparse PPCA's regions of high variance).
  - Preferred form: the band learned INSIDE one fit, not an outer search over caps.
  - Either a per-site per-shell prior variance by marginal likelihood (ARD over shells), or the per-shell half-set agreement
    of the modes as Wiener weights.
  - Briefed to sparse-ppca2.
- 2026-11-08 afternoon: four in-model forms all lost the IgM 10 % fix ([[sparse-sites-inband-band-verdict-2026-11-08]]).
  - A seed needs its band fixed before its first iteration, and a node-wide probe is too coarse.
  - Last try (briefed): per-shell even/odd covariance agreement measured inside the seed balls, pooled into one band before
    phase 1.
  - Otherwise the ladder is final, with its two cost cuts (concurrent first caps, a round-1 null screen).
