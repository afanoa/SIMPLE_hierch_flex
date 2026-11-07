---
name: rnasep-softem-real-data-2026-10-13
description: "RNase P soft-EM loop on 741,592 particles (S1/S2/S4 rows, junk 0), 2 seeds: from the S1/S2/S4 labels the partition settles by it 6-8 at 0.72 agreement with refine3D_states; from flex's lp-20 seed (0.38 start) it settles at 0.41; cross-seed agreement 0.67 (chance 0.34): the S1-like state is seed-robust (~80% mutual), S2/S4 is a continuum cut differently per seed; view NMI 0.007, defocus scrubbed 0.021 -> 0.003; ~3 h per iteration at nparts 10"
metadata:
  type: project
---

Chain `~/autostates/rnasep_softem.sbatch` -> `/mnt/beegfs/elmlund/afan/rnasep_softem/` (p.simple = 25_refine3D_states
relabelled 1,2,4 -> 1,2,3 via `selection infile=states.txt`, junk 3,5 -> 0; pc.simple = merged copy for consensus AND for
flex_pca, which refuses multi-label projects; flex lp 20 npreimages 3 on the softem build -> K=3 store; arm A = hard labels
written into the store with `softem_store.py labels:`; arm B = flex's own delivery; refine3D nstates=3 refine=prob
SIMPLE_SOFTEM=1 SIMPLE_SOFTEM_REFS=sx flex_weights=yes minits 8 maxits 10, nparts 10 x nthr 8). Scorer
`softem/rnasep_softem_score.py W armdir K` (needs W/states.txt, W/poses.txt, W/automask3D.mrc); finalize
`~/autostates/rnasep_finalize.sbatch` (merge_by_view with junk rows kept at 0 + hard M-est + cross-arm map FSC).

**Timing:** iteration 1 = prob_tab over 742k x 3 states, 1 h 38 min on 10 x 8 cores, then the search; ~3 h per
iteration; 8 iterations ~ 26 h; the 2-day walltime killed arm B after 7 iterations (before its merge/final).

**Arm A (S1/S2/S4 seed):** agreement with the refine3D_states labels 0.849 (it 1) -> 0.750 (it 4) -> 0.717 (it 8), populations
233k/253k/255k stable from it 6, median max posterior 0.97, NMI view 0.010 -> 0.007, defocus 0.0002, null 0.0001; loop-map
state differences e/o 0.999/0.99/0.85-0.88. **Arm B (flex seed, start agreement 0.384):** settles at 0.414, populations
211k/360k/170k, posteriors sharpen 0.61 -> 0.95, defocus NMI 0.021 -> 0.003 (the loop scrubs the nuisance the flex seed
carried). **Cross-seed:** A1 <-> B1 174,914 shared (of 233k / 211k); A3 -> B2 almost entirely (197,908); A2 splits B2/B3
(116k/124k); best-permutation agreement 0.671, NMI 0.28. => one seed-robust state (S1-like), S2/S4 = a continuum cut
differently per seed (consistent with the within-state ladders).

**Traps:** (1) merge_by_view's argmax labelled zero-posterior (junk) rows as state 1 -> the first "final hard M-est" of arm
A contained the 190,904 junk particles; fixed (rows with no mass stay 0), finals rerun by rnasep_finalize. (2) refine=neigh with SIMPLE_SOFTEM=1 never finished ONE iteration on the 20k toy (13 h TIMEOUT, job 53804533): do not use neigh
in the loop; prob is the only validated search. (3) the chain's 2-day limit cannot hold two 8-iteration
arms; run arms as separate jobs. Related: [[softem-toyoc-minority-recovery-2026-10-08]],
[[rnasep-s4-four-arm-support-test-2026-10-08]], [[rnasep-final-poses-make-states-separable]].
