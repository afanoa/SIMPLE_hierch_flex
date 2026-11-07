---
name: poles-are-minority-vs-rest-k-includes-rest-2026-10-22
description: "Pole design (user-approved 2026-10-22 ~10:50): the two-class site EM's fixed points are MINORITY-vs-rest states, so a 70/30 switch gives ONE pole; K = accepted poles + a rest state; seed a pole only at posterior >= 0.5 (not argmax); rescue mutual-best twins; poles are per site (another site's difference is invisible inside the mask)"
metadata:
  node_type: memory
  type: project
  originSessionId: 25d47be4-7d1d-48d8-abdf-97b1c3c946ab
  modified: 2026-10-22T14:53:33.010Z
---

Port session's toy calibration of the halves-based pole acceptance (site_rank build 41, `stree_port/halves_toy2`, truth in factorial
`toys_for_port`):
- **Poles are right.** simsites closed arm (8 %): every start reaches the closed class (AUC 0.945, oracle 0.95); twins merge at about 0 SE.
- **Argmax seeding fails.** For 91 % of particles both pole posteriors are < 0.5, so the open bulk lands in the closed pole (precision
  0.12).
- **A 70/30 domain switch gives P = 1** (every start finds the 30 % minority), so K = 1 and the node becomes a leaf; the switch is missed.
- **Lost extreme.** IgM arm 4: twin endpoints miss the merge by a hair and are both rejected.

User approved, now the defaults:
- (a) **Rest state:** a particle joins a pole only at posterior >= 0.5, the rest form a state, so K = P + 1 (P = 1 -> K = 2; P = 0 -> leaf).
  simsites A K = 2 at 0.95 precision; C closed seed precision 0.63; RNase P r1 root K 3 -> 4 (rest 8-10 %).
- (b) **Twin rescue:** a rejected endpoint whose best raw partner is the same other endpoint in both directions joins it, then
  acceptance is re-tested.

Joint K-class EM seeding inside one site mask is NOT the fix (`repro_r1/jointem`): it merges open A/B, whose difference lies outside the
1f1 mask. **Poles are per site;** a difference in another site belongs to a child node with its own discovered site.

Open:
- open A vs B at the root sits 0.7-2.1 SE apart (closed twins about 0 SE, closed vs open 3-21 SE). Whether to require strong evidence
  to split, deferring weak splits to the child, is for the user;
- the rest state's reproducibility: RESOLVED (12:00, build 43). The r1 root gives K = 4 in all three runs, with all four stage-B
  states within noise (0.95-1.01). The rest (19.1-19.6k) reads both arms at about half occupancy (0.57 / 0.56) and is NOT the census's
  middle corner (Jaccard 0.03-0.05). That corner sits inside the open A-open B pair, so middle-band starts are not needed.

**CORRECTION (2026-10-22 ~14:00, port builds 44-45):** the site EM was not converged at the fixed 10 iterations. Run to convergence
(19-42 iterations at r1), the r1 root is P = 2 (closed vs open) in all three runs: open A / B merge, the rest vanishes, and the 3 SE
rule agrees. The K = 4 root was an early-stopping artefact. A converged pole-vs-rest EM finds the node's dominant BIPARTITION, so the
knob-free design is a binary tree of reproducible splits; further structure (open A / B) comes from the child's own site.

**Why:** the site EM models "this population vs everyone else"; reading its fixed points as all the states undercounts.
**How to apply:** count states as accepted poles + rest; never seed by argmax over separately normalised pole posteriors.

Related: [[pole-merge-line-inside-closed-pair-2026-10-22]], [[landscape-must-be-reproducible-and-splits-data-derived]].

**ADOPTED (user 2026-10-23 morning): the two-test pole rule.** An endpoint is a pole only if (a) its hard class (posterior >= 0.5)
holds < 1/2 of the node AND (b) its hard-class seed maps pass the existing cross-half 3-SE test. Otherwise: the next site, then a
2-state stage A / B from fflex, then a leaf. It fixes the empty (toy) and collapsed (RNase P r1 child) modes, which both pass the
halves test on soft maps. The user also approved per-particle CONTRAST in the site EM and stage A / B ("shouldnt that be handled by
contrast": the same thing flex_pca fits) and symmetry support.
