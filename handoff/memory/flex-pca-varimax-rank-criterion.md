---
name: flex-pca-varimax-rank-criterion
description: "mean matched varimax |cos| across disjoint-half refits picks rank (argmax), but the mean has an unmeasured low-rank bias; proposal + Phase-0 battery in doc/for_developers/ideas/flex_pca_halfset_em_varimax_rank_proposal.md"
metadata: 
  node_type: memory
  type: project
  originSessionId: 74eda56d-80ed-4b21-b788-b0776c91d800
  modified: 2026-08-31T16:32:31.919Z
---

2026-08-29 overnight round (R5–R27 in `/Users/atifao/ribo_analysis/FINDINGS.md`), proposal
written 2026-08-31: `doc/for_developers/ideas/flex_pca_halfset_em_varimax_rank_proposal.md`
(supersedes §3 of the independent-halfsets proposal and Change 2 of the staged schedule plan).

**The criterion:** varimax-rotate each disjoint-half basis independently (orthonormalized span,
0.45-box mask, unit-norm, no sqrt-eig weighting), Hungarian-match on |cos|, take the mean;
REFIT at every rank (R22: rank-10's 5th dim reproduces at 0.510, rank-5's at 0.075 — truncating
a generous fit is refuted). Argmax over the ladder picked rank 4 on BOTH 10028 (0.858 vs 0.610,
margin 0.248 — but only two rungs measured) and 10180 (0.691 peak over five rungs, margin only
0.045). Root cause it fixes: eigenvalues 2..k are near-degenerate (~1.3x band), so ordering and
in-span rotation are arbitrary across fits. **Rotation buys reproducibility, not
identification** (0.806→0.985 cross-half; but on full data varimax COSTS ~0.05 identification,
R20: 0.905→0.853). Varimax is low-rank only (destroys at n10/n20). Spatial ICA is
seed-unstable (knob) — rejected; latent ICA is deterministic at K≤5 and finds bimodal axes
(RACK1 0.752) but not continuum ones — readout only, never estimation.

**The unresolved defect:** the mean's implicit threshold is the running mean (pc1-dominated),
and R4's earlier rank-2 numbers (raw 0.987/0.922, mean ~0.955 > 0.858) suggest plain argmax
selects K=2 on 10028 and loses RACK1. Rule B alternative (largest K whose WORST matched pair
clears a floor: 0.771 at n4 vs 0.007 at n5) fixes the low end but needs per-rank floors.
Phase 0 of the proposal (offline, before any Fortran) decides: 10028 rungs 1–3, three mod-4
pairings on 10180, random-rotation + no-signal nulls, matched-axis rerun floors.

**Instrumentation caveat on ALL ladder evidence:** `simple_flex_pca_em_iter.f90:956` keys the
ONLINE stochastic-EM window off worker-local `it` (always 1 under nparts>1), so the window
never advances — the 2026-08-29 runs (ONLINE=1 BATCH=40000 nparts=10) fitted a fixed ~76% of
each shard. One-token fix (`it` → `it_eff`); re-run one rung before trusting §1's numbers.

**Final-stage decision (2026-08-31, user-set constraint):** NO fresh full-N EM anywhere in
production — the paired fits ARE the fit, following the refine3D idiom (alignments/basis never
restart; "final" products are linear replays: `sum_pair_into_sum_rec` accumulator merge,
`combine_eo` = one extra pass). Delivery = frame-align B→A, rotate B's M-step accumulators
(Y·R, RᵀρR), sum four quarter-sets, ONE joint solve with cross-fit-FSC Wiener, then one all-N
embed; escalation is 1–3 warm-started iterations, never a refit. The honest cross-fit FSC is a
first-class versioned artifact (per-shell curves) because future SSNR/τ² shrinkage and shared
frequency marching hang off it — both have exact refine3D precedents (`add_invtausq2rho`,
`set_bp_range3D`). Also: abinitio3D itself turns gold-standard matching OFF
(`GOLD_STD_STAGE=TURNED_OFF` — eo-averaged refs during lp-set marching, eo split kept for
reconstruction/FSC only); refine3D proper routes each particle to its own half's reference.

**Why:** this is the first working GT-free rank rule after eigenvalue order, split-half rho,
the ratio rule (R23, retracted), and in-band FSC counts (saturates k of k) all failed.
**How to apply:** never compare half-set bases by index or signed FSC — use the §5 contract
(mask, unit, SVD-frame start, varimax, Hungarian, span-cosine backstop); read the proposal
before touching rank/halfset code. See [[flex-pca-basis-rank-over-specified]],
[[flex-pca-em-speedup-plan]], [[flex-pca-must-be-knob-free]],
[[flex-pca-halfset-fits-lose-fsc-regularization]].
