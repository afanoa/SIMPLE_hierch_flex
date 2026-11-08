---
name: site-coupling-factorial-verdict-2026-11-07
description: "VERDICT 2026-11-07: sparse sites are factorial across sites, joint within a site. The IgM product of per-site states beats a free 72-cell table by 45 held-out nats and a joint mixture by ~2,000; planted couplings are found (z 15.6); a covariance screen alone false-calls 6/90 null pairs; TRAP: integrating nuisances under their fitted prior invents couplings, so profile them"
metadata:
  node_type: memory
  type: project
  originSessionId: cc23038c-fea9-4b5f-b321-d58eb2942350
  modified: 2026-11-08T03:52:13.403Z
---

site-coupling session (Yan & Pollefeys-style subspace coupling, for the shared tree model's structure);
`hetero_tree_20261028/handoff/site-coupling/{VERDICT,STATUS,results}.md` + scripts/; runs `math_reworks/site_coupling/`.
Method: XD with the posterior precision A on sparse PPCA's per-particle site coordinates (each particle's full G_i, cross-site
blocks included). Held-out likelihood: fit on even, score on odd.
- IgM product-form subset: per-site K comes out 2, 2, 2, 3, 3, the truth's 72 cells. The five-arm product (89 params) beats a
  free 72-cell weight table by 45 nats (z -5.3) and a joint mixture (K 64) by 2,057.
- Full sets: the product beats the joint mixture by 2,030 nats (fix twin) and 1,651 (null twin).
- Planted coupling (subsampled a1-a3, V .32): free table +425 nats, pair test +489 (z 15.6). Null twin, coupled angles
  (theta 1): +180 (z 9.6). 0 false calls on 65 independent pairs; only the weak plants (V .026) are missed.
- qtoyflex (one site): its blocks are coupled under both nulls (free table +373, all-free +545), so it is joint within the site.
- The covariance screen alone fails the null twin: 6 / 90 negative pairs called at 2 SE (4 under a parametric bootstrap), 4 of
  them involving site a5. None reaches the Sidak level, so the screen can shortlist pairs but not decide.
- TRAP: integrating the nuisances under their fitted Gaussian prior gets 5-6 of 10 independent pairs called, and the calls
  survive the bootstrap. Profile the nuisances.

- Final (all jobs done): the null twin's 10 pair tests all favour the product (free table -4.1 to -55.1 nats, z -1.1 to -7.5);
  0 false calls on 70 independent pairs. CAVEAT: refits from another start of the site mixtures keep the decision but move the
  gain by up to 5 paired SE (e.g. a3-a4 -4.1 vs -42.7), so fit site mixtures with restarts before the merge test. The
  covariance screen catches 2 of 4 same-arm IgM block pairs (a2 [4]-[7] z 4.3, [4]-[8] 2.7); joint-within-site stands.

**Why:** decides the shape of the cohesive model the user wants ("a cohesive model over the entire tree"): a forest of per-site
trees, not one joint tree. That also explains the equivalence merge-up: a joint tree duplicates a site's states under other
branches.
**How to apply:** one tree per site (a particle's leaf = the product of its site leaves), a joint tree over each site's blocks;
merge two sites only when their free table beats the product by 2 SE held out; profile nuisances; keep each particle's full G_i.
Related: [[tree-latent-shared-basis-verdict-2026-11-07]], [[sparse-ppca-pursuit-verdict-2026-11-07]],
[[lever3-latent-refs-seeding-verdict-2026-11-07]].
