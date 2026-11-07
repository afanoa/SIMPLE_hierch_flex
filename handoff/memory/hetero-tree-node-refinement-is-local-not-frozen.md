---
name: hetero-tree-node-refinement-is-local-not-frozen
description: "TRAP (user caught me 2026-10-23): hetero_tree DOES re-align particles in every node (stA0 / stage A / stage B = refine3D prob_neigh, local around the parent's poses, lp 15 / 15 / 10 at box_crop 64 / 64 / 88, 3 iterations each, pcg_focus against the consensus background); prob_athres=0.01 makes the prob sampling GREEDY, it does not freeze angles"
metadata:
  node_type: memory
  type: project
  originSessionId: 25d47be4-7d1d-48d8-abdf-97b1c3c946ab
  modified: 2026-10-23T17:14:45.549Z
---

I told the user "the tree never re-aligns particles" (relaying a peer's "stage B freezes the angles"). Both were wrong. In
`simple_hetero_tree.f90` (refine_keys; ITA = ITB = 3, BOXA = 64, BOXB = 88, LPA = 15, LPB = 10, NSPA = 2500, NSPB = 5000, ATHR = 0.01),
every node runs three multi-state refine3D stages with one reference per state:
- stA0 and stage A: lp 15, box_crop 64;
- stage B: lp 10, box_crop 88;
- all with refine=prob_neigh (geom neighbourhood around the current pose), center=no, rec_backend=pcg, pcg_focus=yes with the node
  mask.
`prob_athres` (default 10) is the angle threshold for the probabilistic SAMPLING. At 0.01 the sampler takes the best candidate (greedy);
the neighbourhood search itself still runs.

**What is really limited:**
- the search is local around the parent's poses, so a large consistent offset would not be undone. (CAHRA's "16 deg twist" turned out
  NOT to be an error: the two simulated models differ by +16.2 deg in frame, and the open poses are right in their own frame,
  1.3 deg median; see [[cahra-13991-pose-entanglement-run-2026-10-23]]);
- references share the consensus outside the site mask (pcg_focus), so alignment is mostly core-driven;
- alignment stops at 10-15 A.
The motion session's test: whole-molecule re-refinement at +-10 deg moved the leaf poses only 0.03-0.05 deg, so they are converged.

**How to apply:** before claiming what a SIMPLE stage does, read its keys in the code (refine_keys / common_keys). Never relay a peer's
mechanism claim as fact without checking it.
