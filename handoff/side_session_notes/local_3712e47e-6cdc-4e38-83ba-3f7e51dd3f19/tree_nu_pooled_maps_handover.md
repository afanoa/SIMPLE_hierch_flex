# Non-uniform tree-pooled leaf maps (session tree-nu-pooled, worktree tree-nu-pooled, branch claude/tree-nu-pooled off claude/tree-dist 9005a8c77)

Verdict 2026-11-02 ~01:00 (batch2 time): NO-GO. No Fortran module, no key; the note is
doc/implementation_notes/tree_nu_pooled_maps_2026_11.md on the branch, the prototype scripts/hetero_tree_nu_pooled_proto.py (+ .sbatch).
Work dir on batch2: /mnt/beegfs/elmlund/afan/hetero_tree_20261028/nu_pooled (scripts/, runs/{raw,lp10,lp10sqnc,lp10sq}, logs/).

- Prototype on the saved 692k r17d plain leaves: state_tree port rebuilds the plain-map leaf tree (same sibling pairs as the saved
  stage-B tree), pooled even/odd at all 133 nodes, per leaf the nu_filter unary (whitened Huber cross-half error, tent window 15 A)
  per level of its path, argmin, nu_filter's ordered-label Potts, composition from the chosen pool Wiener-weighted by its FSC.
- Held-out FSC (even-built vs own odd, 40-10 A, median / worst over 67 leaves): plain .943/.759, tree-Wiener emp2 .958/.808, NU
  .946/.774 (10 A-band objective), .928/.727 (raw band), .945/.774 (squared loss + full-map correction). NU never beats emp2.
- Closed corner r1311_s2: plain 1.004/0.577, emp2 1.030/0.594, NU 0.92/0.60 in every variant (99-100 % of the closed-arm mask pools
  with r1311_s1, closed arm 0.77); the own level is chosen on 0 % of the envelope for 64 of 67 leaves. The sibling difference is
  below the 10 A-band half-map noise in a 15 A window, so a local MSE-optimal choice pools it; the tree-Wiener's global per-shell
  weighting keeps it.
- Cost 3-5 CPU-h per variant (12 cores, 26 GB, 22-34 min). Label maps per leaf and run: runs/<run>/labels/unit_NN_level.mrc.
- If a non-uniform final map is wanted later: a local tree-Wiener (its covariance fit per window), not a level selection.
