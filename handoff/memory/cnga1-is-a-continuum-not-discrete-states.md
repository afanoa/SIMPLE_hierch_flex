---
name: cnga1-is-a-continuum-not-discrete-states
description: "2026-09-14: cnga1 has NO discrete states -- dip statistic exactly 0.00000 over 21,652 directions at BOTH d=10 and d=16; the C4 orbit-degeneracy explanation is REFUTED (subspace is not C4-closed, max sv 0.886); box 128 / lp 8 is the measured optimum"
metadata:
  type: project
---

**Multimodality: negative, twice.** `pp4.py <run> 20000 30000 20 20` on the real 245,424-particle
latent: best dip over 20000+1652 directions = **0.00000** at d=10 *and* at d=16 (raw and
deconvolved), valley required to split >=2 % of mass per side. The winning 2-component EM is a
50/50 symmetric split of one cloud (standardised separation 0.504 raw / 0.764 deconv) -- not a
minority cluster. K=3 macro-clusters appear at rank 10, 16 and 24, but they are regions of a
continuum, not islands; only the ~23 k small cluster is stable across ranks
(22902 / 17713 / 24123).

**The C4 explanation is refuted, not merely unproven.** `pgrp` is never passed to flex_pca (it is
command-line-only, `simple_parameters_parse.f90:159`, default `c1` at `simple_parameters.f90:328`),
so nsym=1 and the replication loops at `simple_flex_pca_pcg.f90:463,487,579,599` are inert --
there was never anything to switch off. The subspace is **not** C4-closed: svs of Q^T R90z Q are
0.886 ... 0.018, none >= 0.95, so orbit degeneracy cannot explain the flat spectrum, and C4 irrep
magnitudes are not register-invariant (|lambda1| = 0.738). A hidden register-scrambled asymmetric
mode would show a near-degenerate C4-closed block at the top of the spectrum; three independent
fits show none. The two "conflicting" symmetry numbers are the SAME statistic:
recapture = mean(sv^2) = 4.6417/10 = 0.4641.
**NEVER pass `pgrp` to flex_pca**: the cmdline allowlist is global, so `pgrp=c4` parses silently,
sets nsym=4, and changes every accumulate path with no log line. SIMPLE has no symmetry expansion.
What *is* true: refine3D ran pgrp=c4, so 99.97 % of poses are folded into e1 in [0,90) and
protomer-specific states are scrambled at the data level -- but the concerted symmetric modes carry
3.26 of 10 units of basis power and are themselves a continuum.

**Configuration space is closed** (all at rank 10, |cos|>=0.7 in both halves / mean min |cos|):
box 96 -> 6, 0.774 | **box 128 -> 8, 0.798** | box 160 -> 6, 0.759 | box 160 + lp 6 -> 4, 0.704.
Finer band is strictly worse (cluster count rises to 4 while reproducibility falls: noise).
The `het-resolution(FSC 0.143)` readout is worthless -- it tracks the SOLVER BUDGET, not the data
(10.49 A at maxits_pcg=0, 8.04 at 3 and 8, 7.54 at 2 warm), is initialised to the band edge
(`simple_flex_pca_em_mstep.f90:189,210-220`), and its companion count is right-censored
(10 of 10, 16 of 16, 24 of 24). Neither half can set lp or neigs. See
[[cnga1-continuum-is-the-micelle-halo]] for what the continuum actually is.
