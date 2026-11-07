---
name: hetero-tree-paper-plan-2026-11-06
description: user started a 3D hetero_tree paper (2026-11-06); integration of latent model + Bayesian refinement is the named methods gap; positioning vs 3DVA/RECOVAR
metadata:
  node_type: memory
  type: project
  originSessionId: dc17548d-1b94-4fe8-ad02-93f812b3f2a7
  modified: 2026-11-06T17:13:01.686Z
---

User, 2026-11-06: "Can we start to write a paper? Please look at 3DVA and RECOVAR papers. What do each of our test molecules
add to the paper. At this point it seems we have something that works, but the continuous latent variable model and the
bayesian refinement need to be more tightly integrated." Then: "propose some more molecules to demonstrate on and we can
pursue those as well." An outline was given in chat; decisions pending at that point: integration inside this paper or as
future work, which new molecules to download.

- 2D (flex2D) already has its own IUCr-format draft in the MAIN checkout's untracked doc/manuscripts/flex2D_cls_split_2026_10_03
  (RNase P arms, 10076, 10180, NOT, TRPM4, PfCRT, 10028); the new paper is the 3D tree paper and only points to it.
- Node coupling as built (claude/tree-lever-floor simple_hetero_tree.f90): flex_pca K = 1 at fixed poses seeds states; stage A/B
  refine3D prob_neigh with prob_athres 0.01 (greedy: best pose, best state), each state one map (no latent).
- 3DVA (JSB 2021) and RECOVAR (PNAS 2025) both fit at FIXED consensus poses; neither has a calibrated test of how many states;
  RECOVAR names an alternating pose/conformation loop as future work (its Sec. 2.3). Shared datasets: 10076 and 10180.
- Novelty must not rest on the per-node latent model: flex_pca follows RECOVAR's supplement in places
  ([[recovar-supplement-unit-conventions]]) and SOLVAR already marginalises z with pose refinement
  ([[solvar-marginalizes-z-during-fitting]]).

**Why:** the paper is the user's next deliverable; these facts frame every section.
**How to apply:** draft from the outline; read [[joint-pose-heterogeneity-toy-verdict]], [[pose-state-coupling-toy-verdict]],
[[soft-multistate-fixed-pose-verdict-2026-10-30]] and [[ppinfo-lever-verdict-2026-11-02]] before designing the integrated node.
Headline stays RNase P ([[flex-pca-authorship-and-headline-result]]).
