---
name: data-fission-merge-test-verdict-2026-11-07
description: "VERDICT (2026-11-07): data fission (X1 = X + tau Z, X2 = X - Z/tau; labels from X1, test maps from X2) removes the self-selection that failed lever 4's map-level cross-half merge test: qtoynull merges 4 -> 1 in both bands at tau 1 and 0.6, qtoy2 keeps its real states; cost = labelling SNR (-.08 accuracy at tau 1, -.02 at tau 0.6). TRAP: a purity split (same majority structure, 5 % different conformer admixture) is kept by any map-level test; finer nulls and corr-trimmed cores do not merge it"
metadata:
  type: project
---

Data fission at stage B (fixed poses, state_eval build), 2026-11-07, `handoff/fission/VERDICT.md` (scripts `fis_*.py/.sbatch`,
runs `math_reworks/latent/fission/<key>/`, stacks `fission/stacks/<toy>_tau<t>/X{1,2}/all.mrcs`):
- The toys' noise is white Gaussian (radial power flat to 0.1 % from 10 A to Nyquist), so Z is white with each image's own sd
  (from its 3.3 A-Nyquist power); on real data Z must be shaped by the committed sigma2 group spectrum. Validated: cross-copy
  noise correlation +-.0004 in noise-only shells, corr(X, X1) = corr(X, X2) = .707 at tau 1 (.857 / .513 at tau 0.6).
- The null passes: qtoynull root, X1 labels (random or mixture seeds, 5 rounds) on X2 maps merge 4 -> 1 at 39-10 A (z <= 2.2
  at tau 1, <= 2.7 at tau 0.6) and at 9.4-5.6 A (z <= 0.8); the same labels on X1 maps give lever 4's signature (4 -> 4, z 4-13).
  qtoy2 root keeps every state of different composition (z >= 6.9 / 4.1) and merges a clean duplicate (28 | 28, z 1.5).
- Cost: labelling on X1 loses class accuracy at equal passes: tau 1 .695 vs .772 (closed F1 .80 vs .95), tau 0.6 .751 vs .772
  (closed .90 vs .95), the test copy paying at tau 0.6 (one-state even/odd FSC .92 vs .96). Storage 2x the stack; fission 3 min
  per 60k x 160^2; one fission at the root serves every node.
- TRAP (purity split): two holo states (82 % / 90 % holo, differing by 5 % misassigned open conformers) are kept at z 8.9 by
  the X2 test; finer shuffle nulls (480, 1,920 cells) do not move it (the view cone of the impure state is a symptom); a
  corr-trimmed core removes the junk but not the conformers (z 9.1); only the true-holo cores are equivalent (z 2.1-2.5). No
  map-level test of whole states can call that merge; the merge-up needs per-particle state margins or the tree's
  half-independent re-labelling. A sharp view-cone cut of a state gives a missing-cone reconstruction (even-map sd 30, cFAR
  .005): never cut controls that sharply.
- The hier-ppca frequency split (test band 9.4-5.6 A above the assignment lowpass) also passes the null on X1 but has no power
  on small / low-SNR states (half-map correlation .02-.05): it folds qtoy2's closed state into a junk lump (closed F1 .80 -> .42).

**Why:** the parent asked whether fission beats the self-selection that made the map-level merge test keep five states of a
homogeneous node; it does, at a labelling-SNR cost that tau 0.6 mostly recovers, and the remaining false keep is a different
problem (purity), not noise.

**How to apply:** fission once at the root (two stacks, two projects via update_project ptcl_old_root/ptcl_new_root); every
selection (poses, states, junk reference, leaf labels) on X1; every split-acceptance and merge-up map test on X2 over the full
39-10 A band with the view x defocus shuffle null; final maps from X. Decide purity splits from state margins or half-independent
re-labelling, not from whole-state maps. Related: [[oversplit-junk-ref-merge-verdict-2026-11-07]],
[[cross-half-test-needs-band-above-assignment-2026-11-07]], [[null-acceptance-majority-leaf-2026-10-27]].
