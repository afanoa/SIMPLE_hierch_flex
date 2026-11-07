---
name: cnga1-c4-pose-folding-scrambles-protomers
description: "2026-09-14: C4 pose folding DOES scramble protomer-specific states in cnga1, but it does NOT explain the continuum -- the orbit-degeneracy inference was REFUTED; replaces the note that claimed the continuum is a C4 artifact"
metadata:
  type: project
---

**This note corrects and replaces `cnga1-c4-symmetry-blocks-discrete-states`, whose central claim
was wrong.** That note argued the continuum is a C4 artifact because the eigen-subspace is not
closed under the 90-degree rotation (52-77 % of a rotated mode recaptured). **The inference runs
the wrong way:** orbit degeneracy *requires* closure. If a strong asymmetric mode had been
orbit-copied into four lambda/4 partners, the top of the spectrum would carry a near-degenerate
C4-closed block -- and three independent fits show none (max singular value of Q^T R90z Q = 0.886,
none >= 0.95; ||A^4 - I||_F/sqrt(10) = 0.808). Non-closure is evidence *against* the story it was
offered as evidence for. It also follows that C4 irrep magnitudes |z_m| are not register-invariant
(|lambda_1| = 0.738), so that readout does not exist.
The two "conflicting" symmetry numbers are one statistic: recapture = mean(sv^2) = 4.6417/10 =
0.4641 (matched null 0.041-0.045, z ~ +60; disjoint half fits 0.4291 / 0.4415). The register
genuinely costs rank -- it just is not the cause of the continuum.

**What survives, and is measured:**
- The C4 axis is **z**: a delivered state map self-correlates **+0.961** at 90 deg and +0.948 at
  180 deg about z, vs ~0.52 about x/y (exact `np.rot90`, no interpolation). The consensus
  `rec_final_state01_pproc.mrc` is C4-symmetric to **0.99973**.
- `refine3D_auto` and `abinitio3D` both ran `pgrp=c4`, so **99.97 % of poses are folded into the
  asymmetric unit** (245,349 of 245,424 have e1 in [0,90); 45-deg histogram
  [115737, 129612, 14, 13, 14, 10, 14, 10]), per `simple_sym.f90:168-170` (`eullims(1,2)=360/ncsym`).
  Four molecular orientations scored identically, so **the recorded pose carries no information
  about which protomer is which.**
- Consequence, correctly scoped: **protomer-specific / asymmetric states** (one subunit displaced,
  k-of-4 sites occupied) are scrambled at the data level before any basis is fitted. That is a real
  blind spot. But **concerted C4-symmetric modes are immune** and carry 3.26 of 10 units of basis
  power (pc1 0.712, pc2 0.760, pc5 0.623 symmetric fraction vs 0.25 chance) -- **and they are a
  continuum too.** So symmetry does not explain the result away.
- **Symmetry expansion does not exist in SIMPLE**: `rg 'sym_expand|symmetry_expan|symexpand'` over
  src/ and production/ finds nothing; the only sym programs are `symaxis_search`, `symmetrize_map`,
  `symmetry_test` (`simple_exec_sym.f90:25-29`), all map-level, and `sym%symrandomize`
  (`simple_sym.f90:724-740`) randomises in place. The old note's "route" was not available.

**Hard rule: NEVER pass `pgrp` to flex_pca.** `params%pgrp` defaults to `c1`
(`simple_parameters.f90:328`) and is command-line-only (`simple_parameters_parse.f90:159`); no path
reads it from the project, so `nsym=1` and the replication loops at `simple_flex_pca_pcg.f90:463`,
`:487`, `:579`, `:599` are inert -- replication was never active and there is nothing to switch off.
But the cmdline allowlist is **global, not per-program** (`simple_cmdline.f90:447-451`), so
`pgrp=c4` parses fine, sets `nsym=4`, turns on replication in every accumulate path, and prints no
log line saying so. The banked-adjoint fast path is also gated to `get_nsym()==1`
(`simple_flex_pca_em_iter.f90:231-242`), i.e. production assumes C1 anyway.

For what the continuum actually is, see [[cnga1-continuum-is-the-micelle-halo]]; for the
multimodality and configuration results see [[cnga1-is-a-continuum-not-discrete-states]].
Related: [[cnga1-has-reproducible-continuum]], [[pfcrt-poses-are-the-blocker]].
