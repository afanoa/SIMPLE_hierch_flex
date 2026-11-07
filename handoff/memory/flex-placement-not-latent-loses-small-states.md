---
name: flex-placement-not-latent-loses-small-states
description: "2026-10-05 toySE latent-vs-placement test: flex's GMM macro-cluster fold loses the 10 % state at every pose set (13 purity 0.20-0.27) while the SAME latent holds it (LDA/kNN recall 0.96-1.00); lowest-inertia k-means K=4 on the latent + M-est maps + one-shot image-likelihood E-step = 0.997 at consensus poses; better poses help the latent modestly (Fisher 1.4/2.9/4.4) but are not the lever"
metadata:
  node_type: memory
  type: project
  originSessionId: 2a150d2d-ec74-4838-81a9-912e171235b3
  modified: 2026-10-05T18:43:01.508Z
---

Question (user, 2026-10-05): is there evidence that better poses improve state separation in the heterogeneity model, as the
basis for a pose/state loop? And isn't the non-monotone toySE refit result (acc 0.785 consensus / 0.833 loop / 0.722 true)
a problem in the model?

Test on fwl-c142339 (`~/softem/`, build ~/SIMPLE_softem flex-softem 697d5ccb0; data ~/softem/toySE copied from batch2;
scripts latent_test.py, est_place.sh, kmeans_labels.py, patch_states.py; results logs/est_place2.out, logs/est_km.out).
toySE = toy2 states 14/15/17/13 at 8k/6k/4k/2k, SNR 0.1, box 96 @ 2.44; flex delivery K=4 missed state 13.

| poses | latent: 5-fold LDA / 15-NN acc (state-13 recall) | Fisher | lin R2(theta) | flex GMM fold K=4 acc (13 purity) | E-step vs M-est maps of that delivery |
|---|---|---|---|---|---|
| consensus | 0.973 / 0.985 (0.96 / 0.97) | 1.39 | 0.888 | 0.831 (0.20) | 0.872 (0.20) |
| loop (r3sx) | 0.987 / 0.996 (0.96 / 0.99) | 2.87 | 0.930 | 0.865 (0.24) | 0.892 (0.35) |
| true | 0.998 / 0.999 (1.00 / 1.00) | 4.39 | 0.952 | 0.854 (0.27) | 0.897 (0.34) |
| true, K=6 requested | same latent | | | 0.951 (0.99) | 0.996 (1.00) |

Remedy at CONSENSUS poses (no pose change): lowest-inertia k-means (10 seeds, unsupervised) on flex's own latent -> hard
states -> reconstruct3D nstates=K -> score_states (fixed pose, lp 6, trs 2): K=6 majority-map 0.963 (13 purity 0.999);
**K=4 0.997, every state >= 0.994 pure**. The one-shot E-step only reproduces the references it is given (it cannot invent a
missing state) and is near-perfect once a reference exists.

**Reading:** the latent carries every state at every pose set; poses improve it monotonically but modestly; the loss is in
flex's state delivery (deconvolved-mixture K=10 -> macro-cluster fold), which merges the small state at K=4 regardless of
poses. The "non-monotone with poses" result was this fold, not the model's latent. Pose/state coupling is not the lever;
placement is, and the image-likelihood E-step finishes the job.

**How to apply:** replace or back up the macro-cluster fold with a lowest-inertia k-means (or run the E-step from an
over-clustered start) before building M-estimation references; judge "state separation" on the latent (LDA/kNN) and on
the delivery separately. Related: [[flex3d-soft-em-toy-results-2026-10-05]], [[pose-state-coupling-toy-verdict]],
[[flex-pca-gmm-undoes-equal-mass-placement]], [[flex-pca-state-placement-needs-motion-frame]].
