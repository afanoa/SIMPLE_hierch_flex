---
name: klem-develop-on-50k-subsets
description: "Feedback 2026-10-15: 'are these really slow tests worth it, can they not be expedited' -> develop the latent model on 50k-particle subsets (--nsub) and short schedules; full 200k runs only as confirmation"
metadata:
  type: feedback
---

User (2026-10-15, after the full-scale klem runs had been running 15 min without output): "are these really slow
tests worth it. can they not be expedited somehow".

**Why:** the model question (does the latent separate the conformers) is answered by 40-50k particles as well as by
200k (the oracle used 40k); the per-particle SNR does not change, only the per-cell statistics. Full-scale runs cost
an hour per variant and hide bugs for 15+ minutes (block-buffered logs).

**How to apply:** every development run of klem / windows / E-steps goes first on `--nsub 50000` (or a 50k toy),
schedules of 6-9 cheap iterations, line-buffered logs (`grep --line-buffered`), and visible first-iteration
checkpoints; launch the 200k version only after the fast one says the variant works. Related:
[[compile-and-test-means-batch2-fast-gate]], [[keep-test-loops-minimal]], [[klem-kernel-latent-em-2026-10-15]].
