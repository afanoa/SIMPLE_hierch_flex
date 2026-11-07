---
name: transport-dominated-mor-is-the-right-framing
description: "The flex_pca rank problem is the Kolmogorov n-width barrier from transport-dominated model order reduction, and that literature says to COMPOSE deformation with the linear basis, not replace it"
metadata: 
  node_type: memory
  type: reference
  originSessionId: 98926916-3af3-4bb7-9587-824545bc5ef2
  modified: 2026-08-05T20:30:48.063Z
---

Found 2026-08-05 while searching applied math for better deformation approaches.
This problem has a name outside cryo-EM.

**The barrier.** For transport-dominated problems the Kolmogorov n-width decays as
**O(1/sqrt(n))** instead of exponentially (Greif & Urban 2019). A linear basis
provably cannot compress moving features. That is exactly the rank
over-specification in [[flex-pca-basis-rank-over-specified]], with a lower bound
attached.

**The fix is COMPOSITION, not replacement.** The methods developed there:
- **shifted POD** (Reiss/Schulze/Sesterhenn, arXiv 1512.01985): decompose in a
  CO-MOVING frame; shifting to compensate transport makes the singular spectrum
  decay rapidly.
- **registration-based MOR** (Taddei): learn a registration map aligning snapshots,
  then POD in the registered frame — "registration improves linear compressibility
  of the solution manifold".
- transported subspaces (arXiv 1912.13024); quadratic approximation manifolds
  (Barnett & Farhat, arXiv 2204.02462).

**Every deformation test in the flex_pca investigation pitted deformation AGAINST
the linear basis — the one role this literature says it should not play.** The
prescribed architecture is `deformation o linear`: register to a co-moving frame,
THEN PCA. It also explains why the growth/amplitude term beat the warp term
([[flex-deformation-model-falsified-on-igg]]): growth is not transport, so it was
already linearly compressible; only the transport part needs the co-moving frame.

The cryo-EM literature (RECOVAR included) treats linear-subspace, deep-learning and
deformation-map methods as ALTERNATIVES. Register-then-reduce looks underexplored
there — a possible contribution rather than a catch-up.

**TESTED on IgG 2026-08-05 — NO BENEFIT, as predicted.** Composition at matched
total DOF (warp(Kw)+linear(Kl) jointly orthonormalised vs linear(Kw+Kl)) gave ONE
favourable cell in sixteen (+0.284 vs +0.255 at K=10 on one pair), not replicated
on the second; spending any DOF on warp at K=5 strictly hurts.

**Why, and the condition for it ever to pay:** the n-width theorem binds only when
variation is TRANSPORT-DOMINATED, and three independent measurements say IgG's is
not — the growth/amplitude term beats the warp term at a third of the rank, a
properly fitted nonlinear diffeomorphism explains only ~half the consensus->state
difference, and pseudo-atom DISPLACEMENT columns carry nothing while AMPLITUDE
columns carry it all. At lp=15 A in averaged reconstructions IgG's heterogeneity is
OCCUPANCY CHANGE, not mass moving, so a co-moving frame has nothing to remove.

**Do NOT conclude registration-based MOR fails for cryo-EM — conclude this dataset
cannot test it.** The condition for it to pay is a large clean domain motion
resolved well enough that grad Vbar has internal structure: ribosome ratchet,
rotary ATPase, tens of angstroms at 4-6 A. IgG at 15 A is near the worst case —
small amplitude, blob-like density, occupancy-dominated heterogeneity.

**The original experiment, still not run in its full nonlinear form:** register the GT state maps to a
common frame (SimpleITK demons works — see
[[flex-deformation-model-falsified-on-igg]] §28, and note a hand-written demons
diverged), then check whether the residual differences need FEWER linear
components. If the singular spectrum steepens, build the composition; if not,
transport is not what limits rank here and the n-width framing does not apply.

Also relevant and already in the repo but NOT in memory before: an existing
Python harness at `scripts/flex_models/` implements unbalanced transport (WFR),
metamorphosis, articulated domains, B-spline deformation and elastic modes, with
a proposal doc at
`doc/for_developers/ideas/python_first_flex_heterogeneity_models.md`.
