---
name: multi-motion-fixture-site-product-2026-10-29
description: "mm3q multi-motion fixture (3 independent 2-state domains on RNase P conf01): the tree fragments (1/8 combinations a leaf, motions re-split in 3 branches) but the site layer finds each motion once; product of root site poles = combination F1 .81 vs leaves .65"
metadata:
  node_type: memory
  type: project
  originSessionId: eb6a15a5-b0c0-4ff3-86a2-d4fd7a97f327
  modified: 2026-10-30T00:13:29.802Z
---

2026-10-29, the independent-motions investigation (coordinator local_97df86ed; notes tree_independent_motions_2026_10_29.md).

Fixture (batch2 routes/toy55):
- `make_mm_maps.py`: base RNase P conf01 plus 3 rigid textured domains on the rim. Each domain has 2 states (+-25 deg hinge) with
  independent P = .5 / .65 / .8, so 8 combinations (combo id = 3 digits, e.g. 122).
- `run_mm_sim.sbatch`: the qtoy simulator. QT=mm3q at SNR0 0.0048 = qtoy2's noise (the primary fixture); simmm3 at 0.0012 = the null's
  noise.
- Scorers:
  - `mm_score.py`: leaves, repeats, factor states, merge-up, site layer.
  - `mm_node_score.py`: per node while the tree runs.
  - `mm_lattice.py`: split difference-map cosines.
  - `mm_fassign.sbatch` + `mm_fassign_score.py`: a one-pass stage-B re-assignment, and its score.
- TRAP: qtoy2 was simulated at SNR0 0.0048, not the script's 0.0012 default (only the null fixture used 0.0012).

Result (recipe17d floor-free, 69 core-h):
- 16 leaves; 1 of 8 combinations is a clean leaf. f1 and f3 are each re-split in 3 disjoint branches.
- 5 accepted splits follow no motion (not quality); the merge-up fused none.
- The site layer kept 3 independent root sites, one per motion (pole F1 .82-.92).
- Their pole product: combination F1 mean .81 (8/8 at >= .7) vs leaves .65 (1/8). Plus one stage-B pass: NMI .780 vs .647.
- Leaf-map ICA: .801. A redundant 4th site drops it to .700.

**Why:** the user asked whether a tree is the wrong structure for several independent motions. The partition is, but the site layer
already holds a factorial description.

**How to apply:**
- Proposal: an opt-in product step over the root-level sites (independence filter, then product labels, then an optional pass).
- Never on RNase P, which must stay byte-identical (coordinator). See [[no-region-basis-for-rnasep]]: leaf-map spatial ICA is a toy
  diagnostic only.
- The leaf-map cross-half rank is broken by selection imprint (82/82 on IgM); the factor count needs a particle-level test.
- Related: [[igm-choose-vs-em-replays-2026-10-28]], [[stability-cut-framing-2026-10-28]].
- CONTROLS (same day; site_product_readout.py over every non-RNase-P site-layer tree). Without a gate, 7 null runs get a fake factor:
  their kept root site has 0 hard members.
  - GATE: >= minpop hard members and decisiveness (hard / soft pole size) >= .5. Real sites .81-3.29, fakes <= .40. Set after seeing the
    controls.
  - Results: nulls 0; continua 1; micelle and rare-state 0; qtoy2 1 (quality and junk poles rejected; product NMI .241 << leaves .706,
    closed/holo are leaf-only); IgM 8 real arm angles; OLE 2 motions/seed, reproduced (|r| .89 / 1.00).
  - OLE cross-seed NMI: product .536 ~ root children .523 > leaves .33.
  - One-pass re-assignment HURTS on OLE (.536 -> .462, minority pole 18 -> 43 %). The step should use the site poles with no pass.
- HELD-OUT (43 new sites on trees not used for the .5 cut): 25 real pass, 14 fake rejected, 3 fake pass, 1 real reject.
  - 2 of the false passes are weak sites that merge into real motions by dependence.
  - 1 is a decisive junk pole (IgM r14 r4_s3, dec .71, 94 % junk): add an in-tree junk rule using the root junk screen's site_junk.txt
    (pruned at the end of the run, so read it before).
  - The miss is a weaker copy of a captured motion.
- Tools: site_product_readout.py, gate_heldout.py (routes/toy55).
