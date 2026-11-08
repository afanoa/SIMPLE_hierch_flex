---
name: separability-proposals-verdict-2026-11-08
description: "VERDICT: projection-pursuit split proposals (homoscedastic K=2 XD on y1, best of 2d directions by held-out gain) find pure tails of the missing states but not the states; qtoy2 1|2, IgM a5 40/40 and r131 closed stay unresolved; nulls held by sat5's side test"
metadata:
  node_type: memory
  type: project
  originSessionId: 1218bb2b-a21a-4d25-9e48-c40ddd89ceff
  modified: 2026-11-08T09:55:15.823Z
---

tree-latent, 2026-11-08 (batch2 `handoff/tree-latent/STATUS.md` addendum 5). The proposal added at every node and site:
- a homoscedastic K=2 XD, i.e. a 1D two-component split along the Fisher direction;
- searched on fission copy y1 over 2d directions;
- ranked by held-out gain and proposed above 2 SE;
- tested by sat5 on y1 labels, with Sidak counting 2d per site.
Scripts: `tl_tree10.py`; `tl_mftree10.py` for r131 under rnasep-misfit's nvarB.

The result:
- The fixtures are unchanged and the IgM twin stays one leaf. Its proposals reach z_mid 5-10 and are rejected only by sat5's side
  test.
- At the miss nodes the best direction is the right one, but it hardens only pure tails (qtoy2 r12: 877 at 74 % conformer 2; IgM a5
  r1: 1,562 at 97 % state 1). sat5 makes no call on them.
- r131 under nvarB stays one leaf. Its best direction still tracks log-likelihood (.33).

**Why:** the coordinator hoped separability ranking would fix the misses that per-subtree basis refinement could not
([[nested-basis-refinement-verdict-2026-11-08]]). It does not: the missing pairs overlap along every single direction.

**How to apply:**
- Don't add 1D separability proposals to the port. The open problem is splitting two overlapping states as wholes.
- On nulls, sat5's side test, not the multiplicity bar, is what rejects projection-pursuit splits.

Related: [[tree-latent-verdict-2026-11-07]], [[real-data-fission-needs-misfit-model-2026-11-08]].
