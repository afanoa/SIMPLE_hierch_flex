---
name: site-chain-verdict-2026-11-08
description: "VERDICT: site chain (sparse sites -> K=2 tied seed -> one masked cross-half round -> GAPw) recovers IgM's balanced arms (.348 / .148, nulls 0) but does NOT split RNase P r131's 7 % closed minority: no block DISCRETE, near-half groups, maps = node mean .89 / .59; v8 per-site table built (commit waits on the band rule)"
metadata:
  node_type: memory
  type: project
  originSessionId: d63181a2-43e0-44aa-a87c-6f7c55606053
  modified: 2026-11-08T10:36:58.052Z
---

site-chain session, 2026-11-07/08. Handoff `hetero_tree_20261028/handoff/site-chain/{STATUS,VERDICT}.md`, runs
`math_reworks/site_chain/`.

- What carries it: the seed, not the round. sklearn's tied recipe (k-means++ start, stop at |dLL| < 1e-3) keeps the state split;
  converged fits lose it (5 % a3 .164 -> .005). refine3D's state_eval round adds nothing over site_split's own round.
- IgM: the worker on lp-20 fits accepts 5 / 5 arms at 10 % (product .348 vs lever 3 .338) and 4 / 5 at 5 % (.148). The nulls
  accept 0.
- In the trees at lp 15 the sparse fit finds no IgM site. The cause is the band, not the seed geometry; the commit waits on
  sparse-ppca2's data-derived band rule ([[band-must-be-data-derived]]).
- v8 adds a per-site factorial table beside the tree: per site, the lowest-GAPw DISCRETE block's labels and maps, in
  `<root>/factorial_chain/node_<id>/`. Its parity with v7 is exact.
- The commit candidate is v8 + the coords-read fix (b10 on 5b17de846), which gives .349 / .145
  ([[simple-nlines-skips-comment-lines-trap]]).
- r131 closed (the HARD GATE state, [[rnasep-closed-state-must-be-retained]]), b10 worker on rnasep-misfit's native fits with the
  noise model on (9 blocks) and off (18):
  - no block is DISCRETE;
  - every split is a near-half (17-28k of 44k);
  - the purest group holds .115 closed (1.7x the .069 share; .100 with the noise model off);
  - its maps read .887-.896 / .59, the node mean, so own data already fails the corner.

- Two levels (those groups as child nodes, own mean, noise model on): the sparse fits find no surviving site, so the chain has
  nothing to test.
- Where the production tree's r131 corner comes from (gates pp_off_r131 / nv_r131, node r1): it crosses .95 in stage A, the
  iterated 2-state refine3D (prob_neigh, lp 15, 3 its) seeded by the fflex groups.
  - Its r1_s2 labels, rebuilt at FIXED root poses, still read .972 / .963 at only .085 / .106 reference-row purity.
  - So the partition (iterated likelihood selection) is the ingredient, not pose re-search or PCG.
  - Reference-row purity does not predict the readout.
  - On own data the corner is mostly imprint ([[extremes-fission-verdict-2026-11-07]]: a corner leaf at .93 honestly).

**Why:** the K=2 tied seed halves a node whose minority state has no second mode, and the full-image round keeps the halves. The
hybrid "shared basis -> seeds -> round" resolves balanced binary states, not a 7 % tail.

**How to apply:** for RNase P's closed state, don't expect the chain or any balanced K=2 seed to isolate it. A minority needs a tail
or unbalanced proposal, and its honest maps need a fission re-read ([[extremes-fission-verdict-2026-11-07]]). Keep the chain for
balanced site states; batch2 does not enforce memory, so size its worker from MaxRSS ([[batch2-memory-not-enforced]]).
- 2026-11-08 afternoon, b11 = b10 + band3 (guide 7.10 thread 3):
  - The gate battery passes 14 / 14, and the r131 corner holds on own data. The chain fired on 7 gate nodes, keeping 0 sites each
    time: the battery shows no harm, not benefit.
  - IgM stock twins: 10 % fix 2 leaves (combo NMI .034), 5 % fix .001, nulls 1 leaf.
  - The 4 IgM chain trees are still running; `sc11_igm_score.sbatch` scores them. Commit only if IgM gains in the tree.
