---
name: per-node-junk-rule-prototype-2026-10-22
description: "per-node junk rule PROTOTYPED + validated (junk_score.py xd3): best-state amplitude NUM/DEN vs every converged state, 3-component heteroscedastic mixture, junk = >= 2 sd below the heaviest; per-defocus-bin mixture REJECTED (splits real-data tails, 8-34 % false flags); real r1/r2/r4/r6 have no junk, f692st root 1.4 %"
metadata:
  node_type: memory
  type: project
  originSessionId: 1ca67e82-4d15-4725-b95d-0ee5719be9f4
  modified: 2026-10-22T19:02:06.785Z
---

Prototype for the port's site_rank seeding step. Code: worktree laughing-carson-4e8ca1, `scripts/cls_expansion_proto/toy55/junk_score.py`
(rule = `junk_flags_xd3.txt`), `junk_test2.sbatch`, `junk_recut.py`. Note: `doc/implementation_notes/per_node_junk_detection_2026_10_22.md`.
Outputs: rnasep_klem/factorial/junk_test/.
- Score: every node particle against every state, cross-half refs, a = NUM/DEN at the best state (genuine ~1, no signal 0), v = 1/DEN.
- Cut: a THREE-component heteroscedastic mixture (junk + genuine core + broad tail); genuine = the heaviest component; junk = >= 2 sd below it,
  weight >= 0.5 %.
- Toys (20 % junk): precision 0.96-0.99, recall 0.65-0.79 (half-molecules ~0.45); the rare closed 8 % is flagged at the base rate.
- At convergence: junk sits inside the decisive classes and is still found.
- TRAPS:
  - the per-defocus-bin 2-Gaussian mixture splits unimodal real distributions (r1, f692st r2/r4: 8-34 % "junk", reproducible across
    runs, rebuilt = a normal map);
  - a 2-component fit misses junk when genuine amplitudes are heavy-tailed (f692st root with its 8 states).
- Real data:
  - r1 at convergence (conv_r1 p1-p3), r2, r4, r6: no junk population;
  - f692st root: 9,137 flagged (1.4 %), whose rebuild has no half-map FSC signal (consensus correlation 0.03); an equal-size random
    unflagged set gives 18 A and 0.62.
- Halves (fit on one half, apply to all): Jaccard 0.92-0.99; the user approved the rule as DEFAULT ON, the port implements it.
- The 107k "neither position" particles fit their own states (0 flagged at their nodes, 2.85 % at the root): a junk STATE, not junk particles.

**How to apply:** use the xd3 rule for junk particles; review junk STATES via the per-state amplitude quantiles, never auto-remove.
Related: [[junk-detected-per-node-not-only-root-2026-10-22]], [[pole-census-synthetic-validation-2026-10-22]], [[never-overwrite-scripts-under-running-jobs]].
