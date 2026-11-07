---
name: stored-weights-pipeline-toys-then-not
description: "2026-10-05 user direction for the flex 3D soft-EM / M-estimation work on fwl-c142339: test the 'store better upstream weights' pipeline on toys first; the real-data target is NOT (well separated micro-states) once its refinement is good, not the spliceosome"
metadata:
  node_type: memory
  type: project
  originSessionId: 2a150d2d-ec74-4838-81a9-912e171235b3
  modified: 2026-10-05T18:57:11.880Z
---

User, 2026-10-05: "on toys for now. but NOT is a good example once I get a good refinement because it has well
separated micro states."

Context: the pipeline under test = flex latent -> lowest-inertia k-means (over-clustered start) -> hard states ->
M-estimation maps -> one cross-half image-likelihood E-step (score_states) -> calibrated posteriors stored as the
project's flex weights -> refine3D_states flex=no flex_weights=yes (+ consensus/deviation references) -> per-state
refine3D_auto. Development box fwl-c142339: `~/SIMPLE_softem` (branch flex-softem), `~/softem/` (toySE, scripts
arm.sh / r3s.sh / r3s_from.sh / est_place.sh / latent_test.py / kmeans_labels.py / patch_states.py / score.py).

**Why:** NOT's refine3D_states states were smeared and its consensus view-biased (cFAR 0.07-0.21); the user is fixing
the sampling/refinement first (see [[not-sampling-fix-reruns-2026-10-04]]); the spliceosome's five states differ only
at >= 20 A, so it cannot show a resolution gain.
**How to apply:** keep real-data runs off until the user says the NOT refinement is good; then take NOT's consensus
project, run the stored-weights pipeline there (distributed on batch2), score by state distinctness (deviation and
pair-difference reproducibility at 20 A and finer), per-state resolution vs consensus, and the view/defocus null.
Related: [[flex-placement-not-latent-loses-small-states]], [[flex3d-soft-em-toy-results-2026-10-05]],
[[not-ribosome-flexpcg-plan-2026-10-02]].
