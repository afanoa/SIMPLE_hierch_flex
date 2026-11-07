---
name: flex-pca-em-speedup-plan
description: The 2026-08-27 EM speedup review reopened convergence+SQUAREM (the 08-18 closure premise fell) and set the GPU target as a single-process resident loop
metadata: 
  node_type: memory
  type: project
  originSessionId: 0ab31dc7-e373-49b7-ba89-3571851ff6e4
  modified: 2026-08-27T15:10:00.120Z
---

Plan lives at `doc/for_developers/ideas/flex_pca_em_speedup_review.md` (2026-08-27). Key facts
that reverse or anchor earlier conclusions:

- The 2026-08-18 "do not build SQUAREM or a convergence criterion" verdict rested on a
  5-iteration operating point. That premise fell: the winning ribosome recipe runs 25 probe
  iterations, the moment path is retired (EM is the only estimator), and the Gamma-floor rank
  criterion requires convergence that no run has reached. The acceleration note's target
  architecture (two independent half-fits + cross-fit softdim stop + pooled delivery, ~+11%
  cost since each fit sees N/2) is revived as Phase 1.
- Amdahl anchor: at 10180/nparts=10 the master+qsys section is 6.4 s of a 21 s iteration (66%+
  at small scale); worker time is 78% M-scatter, 18% E-step. Worker cost is strictly linear in
  window size; the master floor is N-independent — capping windows saturates at ~2x unless the
  serial section is deleted (single-process GPU loop, G8).
- ONLINE windows currently stand down the fused GPU E-step and banked adjoint (measured
  representation collapse 2026-08-19, device state keyed to a fixed particle set); the polar
  stage-2 keying pattern (per-stage/per-iteration/per-batch state only) is the required fix.
- GPU strategy is Blackwell-or-CPU (Turing measured 1.7x slower than CPU on FP64); nothing GPU
  has ever been benchmarked end-to-end under the current recipe — that benchmark is Phase 0
  item 1, and all §5 speedups (~10-25x/iteration) are arithmetic until it exists.

Related: [[flex-pca-ribosome-winning-recipe]], [[flex-pca-probe-iters-and-whitening]],
[[remote-build-recipe-fwl-c143766]], [[flex-pca-must-be-knob-free]]
