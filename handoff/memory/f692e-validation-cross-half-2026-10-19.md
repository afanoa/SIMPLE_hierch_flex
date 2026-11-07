---
name: f692e-validation-cross-half-2026-10-19
description: "user 2026-10-19: validate the f692e states with the cross-half pass but RETAIN the closed state; splitbench6 (identity pairing, --protect labels never merge/prune with unprotected, frozen membership in the cross-half EM) on candidates = f692e leaves + closed sub-states mapped onto the rnasep692 halves"
metadata:
  type: project
---

Design (2026-10-19 23:5x): discovery stays the f692e route on all particles (boards at family scale found the closed
state; half-size trees do not: [[klem-state-graph-umap-and-real-halves-2026-10-19]]); validation = the cross-half pass
on the independent halves. toy55/f692e_candidates.py: labels per half row = f692e consolidated_12_fullw_region leaf
(63) except zoom (consolidated_12_fullw_zoom) state 1 -> 64 (pure closed, 0.94, 4,896 ptcls: 2,629 A / 2,267 B) and
zoom 143/77/30 -> 65 (closed-near 0.83-0.87, 7,123); leaves < 1500 per half -> 0 (23 of 63, left out: the bench loads
labelled rows only), so 40 leaves + 2 closed groups, ~280k per half. splitbench6.py = splitbench5 + --pair identity +
--protect 64:65 (protected merge only with each other, never pruned or folded into, membership frozen in the EM,
never drained; logs "would MERGE / survives without protection" per protected label and "protected group -> final
state N"). f692e_validate.sbatch (PROTECT, XS): jobs 54009100 (xem_prot) and 54009101 (xem_noprot, the data-alone
answer), out rnasep_klem/dtree_rnasep692_f692e/<XS>; then state_recs.sbatch MODE=rnasep692 TAG=f692e FIN=<XS> for the
closed-arm readout. CAVEAT: f692e labels came from boards on BOTH halves, so the first merge test is biased toward
"reproducible" (shared selection noise); the cross-half EM and the re-test after it are the stricter decision.
RESULT (consensus poses, 2026-10-20 01:15): pre-EM test: no merges among the 42 groups (pair stats 0.24-0.39 vs
threshold 0.13); closed 64 vs nearest unprotected 0.48, 65 at 0.38 -> "survives without protection" (biased, see caveat).
Cross-half EM changes 73 % of labels, 18 small groups drain -> 24 states. UNPROTECTED: 64/65 are NOT drained, they GROW
(weight 0.009 -> 0.032, ~8.5k per half) and turn intermediate (closed-arm 0.71/0.73; known-closed proxy 19.7 % / 9.3 %,
vs 47 % in the frozen group) -> the EM pulls near-closed particles in. PROTECTED: frozen state 24 = 64 (2,629/2,267,
proxy 47 %), state 22 = 65 (proxy 13 %); re-test 0.66 / 0.45 vs 0.13. Family-pose runs: 54011181/182 (xem_*_famposes).
Post-EM re-test is inflated by the EM's cross-half coupling (stats 0.74-0.84 after vs 0.24-0.39 before); splitbench7
--retest_beyond_em tests only shells beyond the EM band. Null (null_xem.sbatch: one conformer in 4 random groups):
pre-EM K* 1, forced EM collapses to 1 group -> the EM does not invent states on homogeneous data.
FAMILY POSES (f692e family refine3D_states poses, 93 % of particles; 54011181/182, 02:23/02:34): 28 states (vs 24 at
consensus); protected 64 -> state 24, closed-arm 0.79 (vs 0.73), proxy 47 % / 43 %; 65 -> 21 (13 % / 19 %); UNPROTECTED:
closed particles gather in state 18 at 25 % / 25 % proxy (vs 19.7 % / 9.3 % split over two states at consensus poses),
closed-arm 0.76. Pre-EM: 64 vs nearest 0.53 (threshold 0.14). -> refined poses help retain the closed state; without
protection it is still diluted (~4x the frozen group's size), never a pure closed map in the bench readout.
