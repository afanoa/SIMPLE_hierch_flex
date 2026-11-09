---
name: textbook-flex-weighting-production-level-2026-11-09
description: Bishop-Tipping M-step-only particle weights in flex_pca (nodelevers b6cb5f258, install5) validated (no shrinkage); at r13 -> r131 the soft child gives x1.09-1.20 (production level, 2 of 3 seeds agree), the tempered x1.4 was an artefact; confident seeding = hard exactly (k-means sets the partition)
metadata:
  type: project
---

node_levers SOFT_LOOP.md (2026-11-09): five arms at r131: hard x1.00; sqrt(w)-tempered x1.37-1.42 (artefact); members-at-1 x1.27-1.34
(8 states kept: the state count is not the lever); confident seeding (SIMPLE_FLEX_PLACE_ROWS, a86b61e20) identical to hard to the last
label because GMM AUTO's k-means over all members re-forms the regions; textbook weighting x1.20 / 1.11 / 1.09, Jaccard .87 (s2/s3)
.69 (s1), latent AUROC .714 unchanged; qtoyflex r12 one leaf in every arm. Validations on install5: unset untouched (2e-5), weights 1
= unset (1e-8), hard subset and members-1/extras-0 reproduce the child fit (|corr| .96-.98), corr(w,|z|) .02. A labelled
zero-weight row makes the weight-set publish refuse (w = 0 corner; harness floors at .01).

**Why:** the user wanted the soft child's extras weighted the textbook way before judging soft vs hard.

**How to apply:** weighted soft children are not a lever at r131 (production-level at best); do not pursue further without a new
idea. Related: [[sqrt-weight-on-planes-tempers-flex-estep-2026-11-09]], [[node-levers-replay-verdict-2026-11-08]].
