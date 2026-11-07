---
name: junk-free-component-needs-scale-fair-maps-2026-10-21
description: "junk free components at fixed poses work (precision 0.99, blob 1.00) ONLY with a per-shell data-fitted reference scale and free maps reconstructed exactly like the states; shrink-to-consensus free maps eat 86 % of real particles"
metadata:
  node_type: memory
  type: project
  originSessionId: 25d47be4-7d1d-48d8-abdf-97b1c3c946ab
  modified: 2026-10-21T17:46:51.434Z
---

Prototype junk_em3.py (routes/toy55; runs rnasep_klem/junk_em{,2,3}), null_50k / two_50k, fixed refine3D poses, each half scored
against the other half's references, klem start (free 1 = consensus lp 35 A, free 2 = half-contrast consensus, 2 % of the mass, used in
iteration 1 only), band to 12 A, per-shell reference scale beta_k = sum Re(P conj X) / sum CTF^2 |P|^2 over the real-state members:
- truth junk hidden in random real states: free comps 87-89 % recall at precision 0.99 (blob 1.00, broken half-molecule 0.80-0.83) on
  null_50k; two_50k blob 0.99, broken 0.28 + a separate 93 %-broken STATE (partial composition, as splitbench15's SSNR rule says);
  removing the junk lifts NMI(conformer; state) 0.40 -> 0.54.
- refined states as they are (h0): free comps stay EMPTY, the junk already sits in junk-majority states (0.94-0.99) -> the junk
  decision then belongs to the unit level (low SSNR).
Traps found on the way: (1) a zero-map "noise" component catches nothing (junk correlates with the molecule at low res);
(2) free maps shrunk toward the consensus have ~0.75 contrast vs ~0.5 for the regularised states and win EVERY particle (86 %,
precision 0.07); (3) without beta the high shells (pose error) penalise sharp references: recall 0.11; (4) em_lim 20 leaves too
little high band for junk to lose on.
**Why:** the user remembered klem's junk component as robust and wants it in the node classification.
**How to apply:** any extra component competing in a likelihood (junk, minority) must have the same regulariser/contrast as the
states and a fitted per-shell scale; see [[klem-six-conformer-700k-and-junk-flag-2026-10-19]], [[klem-integrated-hierarchical-model-vision]].
