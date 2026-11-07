---
name: flex-pca-bagging-beats-single-fit
description: "Pooling two independently-fitted eigenbases beats the best single fit at the same rank on IgG, and closes the deformation-encoder question"
metadata: 
  node_type: memory
  type: project
  originSessionId: 98926916-3af3-4bb7-9587-824545bc5ef2
  modified: 2026-08-05T19:56:18.804Z
---

2026-08-05. Nominally equivalent linear fits vary **5x** in explanatory power
([[flex-pca-basis-rank-over-specified]]), so basis VARIANCE, not representational
power, is the binding problem. Bagging tests the cheap fix.

N=2 independent held-out fits (basisrep basis A on eo=0, basis B on eo=1), scored
on three GT bin pairs with the cross-halfset estimator, Vbar deflated:

| pair (ceiling) | eigA | eigB | mean | bag-pred | bag-sub5 (rank-matched) |
|---|---|---|---|---|---|
| b00->b80 (.235) | +0.256 | +0.160 | +0.208 | +0.250 | **+0.265** |
| h00a->h80a (.337) | +0.254 | +0.229 | +0.242 | +0.267 | +0.248 |
| h00b->h80b (.103) | +0.166 | +0.129 | +0.147 | +0.195 | **+0.203** |

**Bagging beats the mean of the individual fits on all three pairs (+27 %, +10 %,
+38 %), and the RANK-MATCHED pooled basis (`bag-sub5` = top-5 singular directions
of the pooled [A|B]) beats the BEST single fit on 2 of 3** (bag-sub10 takes the
third). A better basis at the same rank, from two fits, at zero modelling cost.

**This closes the deformation-encoder question.** Bagged rank-5
(+0.265/+0.248/+0.203) beats rigid-domain at rank 18 (+0.149/+0.211/+0.046)
decisively. The generated basis's last advantage — zero fit-to-fit variance — is
obtained more cheaply by averaging linear fits without surrendering
representational power. See [[flex-deformation-model-falsified-on-igg]].

**To implement:** fit on K disjoint particle subsets, pool the K*k basis volumes,
keep the top k singular directions. `heldout=yes` already fits two bases on
disjoint halfsets, so K=2 needs only the pooling step.

**Caveats before committing:** N=2 is the MINIMUM possible bag, so these gains are
a lower bound and the saturating K is unmeasured; cost is linear in K (each fit is
a full column stage); the 38 % basis-error term
([[flex-pca-deconvolved-mixture-works]]) is the ceiling on what bagging can return;
and all three pairs are one dataset and one contrast, so confirm elsewhere before
trusting the magnitude.
