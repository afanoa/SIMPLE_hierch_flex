---
name: softem-toyoc-minority-recovery-2026-10-08
description: "toyOC (RNase P region-1 swing, 13% minority, 20k, SNR 0.1): the cross-half soft-EM loop recovers the minority EXACTLY (purity 0.47 -> 1.00, fraction 0.28 -> 0.130) in 2 E-steps from a half-contaminated lp-20 flex seed where the linear latent could not separate it; frozen M-est cannot move; BUT the loop fills redundant states with a VIEW split (NMI 0.23 vs null 0.001, e/o 0.99), so K must be pruned and the view-NMI control is mandatory"
metadata:
  type: project
---

Chain `/mnt/beegfs/elmlund/afan/softem/toyoc.sbatch` (env AMP, FRAC, TAG, N, SNR), dirs `softem/toyOC30`, `toyOC20`,
`toyOC20null`; scorer `softem/toyoc_score.py <toy> <flex8|flex20> [frozen|sx|sx/final_mest]`; store writer
`softem/softem_store.py <flexdir> <armdir> <K> [restore]` (writes the last softem_iterNNN posteriors into the flex-softem
store: int64 header(6), digest = row layout only, so overwriting weights is safe); `softem/toyoc_final.sbatch` (env D F K).
Toy: S2 consensus (rnasep_hires/S2) with G10 region 1 swung by AMP deg (toy2_make_vols, levels 3), box 96 @ 2.44, centre
state 10 (87%, 17,400) + swung state 19 (13%, 2,600), SNR 0.1, true poses -> 3 prob refine3D its at lp 8. Flex fits with the
SIMPLE_softem build (the bandsort build writes NO weight store; the softem build's flex_pca must deliver K states that the
arms then use: nstates=K read from the delivery, not 2).

| start | minority recall / purity / est. fraction (true 0.130) |
|---|---|
| 30 deg, flex lp 8 | 1.000 / 1.000 / 0.130 (latent alone solves it) |
| 20 deg, flex lp 8 | 1.000 / 0.227 / 0.572 (swallowed) |
| 20 deg, flex lp 20 | 0.999 / 0.471 / 0.276 (coarse band halves the contamination) |
| 20 deg, lp 20 + frozen M-est | unchanged by construction |
| 20 deg, lp 20 + sx loop it 1 / 2 / 3 | 0.990/0.893/0.144 -> 0.997/1.000/0.130 -> 1.000/1.000/0.130 |

Maps (FSC of minority-minus-centre vs the true difference, 40-20 / 20-12 / 12-8 A): 30 deg frozen .897/.767/.552, sx
.905/.785/.566. 20 deg frozen .663/.554/.390; sx final maps (shrunk refs) .650/.266/.235 -- NOT a loop failure: the loop
split the 17,400 centre particles between the two redundant states BY VIEW (NMI view cones 0.235, defocus 0.009, null
0.001; state1-state2 difference e/o 0.99/0.99/0.94), so the pooled centre is poisoned; minority minus state 1 alone =
.757/.587/.417, better than frozen in every band. "minority minus the TRUE centre" is useless (filter mismatch, 0.10-0.25
even at 30 deg).

**How to apply:** (1) the loop is a classifier: it cannot invent a state absent from its references but it purifies a
2:1 contaminated seed to exact in two cross-half E-steps; (2) seed at the band where the motion is quasi-linear (lp 20
here), not the fit band; (3) redundant states get filled with view splits that reproduce across halves -- e/o
reproducibility is NOT a validity test; run NMI(state, view cones) against a permutation null on every delivered state
and merge/prune before trusting maps; (4) final maps from the final weights with an unshrunk M-est (softem_store.py +
reconstruct3D flex_weights=yes), difference against the best single reference state. Pending: toyOC20 lp-8 arms, toyOC30
lp-20 arms, the homogeneous null. Related: [[flex-placement-not-latent-loses-small-states]],
[[flex3d-soft-em-toy-results-2026-10-05]], [[pose-state-coupling-toy-verdict]], [[rnasep-s4-four-arm-support-test-2026-10-08]].

**30 deg from the contaminated lp-20 seed (purity 0.42, fraction 0.31):** sx loop -> 0.98/0.132 at it 1, exact (1.000/0.130) from it 2;
maps minority-minus-pooled-centre vs truth .898/.772/.556 vs frozen .821/.709/.480. Here the redundant centre state DRAINED
(459 of 17,400) instead of view-splitting (NMI view 0.026 all, 0.057 centre-only, null 0.002). So the loop's treatment of a
redundant state is unstable between drain and nuisance split; prune K to islands+1 before the loop, or merge after.

**Homogeneous null (`toyOC20null`, centre state only, same N/SNR/pipeline):** the flex lp-20 start on a homogeneous sample is
already a VIEW partition (3 clusters 3653/8925/7422, NMI view 0.228, defocus 0.001); the sx loop drains one state (0 by it 3)
and keeps a 2-way view split (11199/8801, NMI 0.152, null 0.001). Its state differences reproduce across halves at
0.95-0.99 (40-20 A), 0.67-0.99 (20-12), 0.27-0.96 (12-8) for BOTH frozen and sx arms. => even/odd reproducibility of a state
difference carries no information about conformational validity on this pipeline; the real minority (view NMI 0.03) and the
null's view split (0.15-0.23) are separated ONLY by the view/defocus NMI against a permutation null. Make that control
part of every delivery.

**Swallowed seed (20 deg, lp 8 start: purity 0.23, fraction 0.57), rerun 2026-10-09:** the loop improves monotonically but
slowly and refine3D's convergence stop ended it after 3 E-steps: purity 0.31 -> 0.36 -> 0.40, fraction 0.41 -> 0.36 -> 0.32,
recall stays >= 0.97; maps .50/.36/.19 vs frozen .41/.31/.11. So the reach of the loop is between a 2:1 seed (exact in 2
E-steps) and a 4:1 seed (not converged in 3); for bad seeds force more iterations (minits) and watch the trend.
**Null, lp 8 start:** the loop drains the spare states on a homogeneous sample (4326/13388/2286 -> 0/19415/585 by it 5,
view NMI 0.24 -> 0.08): on a null it collapses toward one state.

**Spare-state merge rule (2026-10-09, `softem/merge_by_view.py <toy> <armdir>`):** greedy: for every pair of surviving states
the NMI between the pair's binary split and 30-deg view cones over a 5-permutation null; merge the most view-dependent pair
(ratio > 10 and NMI > 0.02, or a near-empty member), relabel, recompute, stop. Pairwise testing without recomputation
fails (a view-split product makes every pair involving it look view-dependent, ratios 16-44 even vs the real minority);
greedy + recompute gives: 20-deg lp20 K 3->2 (17400/2600, minority recall 1.000 purity 1.000), 30-deg lp20 K 3->2 exact,
both nulls K->1. On the unconverged swallowed-seed arm (purity 0.40) it merges everything, as it should for a split that
is still half view. Real-split signature: view NMI ~0.002 (ratio ~1); view-split signature: 0.05-0.24 (ratio 18-206).
Real-data run (rnasep_softem, 742k, nparts 10) was cancelled after one 2-h iteration: develop on the toy first (user).

**Forced iterations on the swallowed seed (minits=6, `softem/toyoc_arm.sbatch`):** purity 0.31/0.36/0.40/0.46/0.57/0.72 over
6 E-steps, fraction 0.41 -> 0.18, recall >= 0.97 throughout, maps .705/.508/.342 (3 E-steps gave .50/.36/.19). A 4:1
seed converges, just slowly (~8-10 E-steps); the default refine3D stop rule is wrong for the loop -> always set minits.

**END-TO-END PASS (2026-10-09, `toyOC25e2e`, fresh 25-deg toy, one job: `toyoc.sbatch AMP=25 TAG=e2e LPS=20 STS=flex20
ARMS=sx MINITS=6 MAXITS=10 E2E=1`):** lp-20 seed purity 0.25 / fraction 0.52 (a 4:1 seed) -> sx loop purity 0.42 / 0.76 /
0.999 / 1.000 by E-steps 1-4, fraction 0.130; spare state drained to 0; merge rule K 3 -> 2 (17400 / 2600); final hard
M-est from merged labels: minority-minus-others vs truth .901/.802/.622, e/o .986/.915/.647. Swallowed lp-8 seed with
minits=10: purity 1.000, fraction 0.129, recall 0.995 by E-step 8. VALIDATED RECIPE: coarse-band flex seed (lp 20) ->
S0 (reconstruct3D flex_weights=yes) -> refine3D nstates=K SIMPLE_SOFTEM=1 SIMPLE_SOFTEM_REFS=sx flex_weights=yes refine=prob
minits 8-10 -> merge_by_view.py -> selection infile=merged_labels.txt -> reconstruct3D nstates=K' -> e2e_score.py.
