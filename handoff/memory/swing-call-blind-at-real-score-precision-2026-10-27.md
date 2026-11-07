---
name: swing-call-blind-at-real-score-precision-2026-10-27
description: "WHY real swings call switch: node_motion nb falls ~3x per 0.2 drop in score rho (swing30 0.23/0.07/0.036/0.020 at rho 1/.8/.7/.6); real rho 0.5-0.7 -> switch range; deliver by axis separability (mid share) + refined-evidence minpop poles; raw tails fail the RNase P closed state"
metadata:
  node_type: memory
  type: project
  originSessionId: 4ec828f0-469c-44cc-b5b4-4bfed9c0e765
  modified: 2026-10-27T19:35:29.272Z
---

2026-10-27, session swing-switch-extremes-58027e (write-up: that worktree's
doc/implementation_notes/planned/hetero_tree_swing_switch_extremes_2026_10_27.md; cluster swing_switch_20261027/).

- Diagnosis: the blend test's non-blend share nb needs a precise per-particle score.
  - Noise-free maps160/swing30 (30 A tip): nb 0.23 / 0.07 / 0.036 / 0.020 / 0.010 at rho 1 / 0.8 / 0.7 / 0.6 / 0.5.
  - swing60: 0.034 at rho 0.6.
  - Real rho is 0.5-0.7, so real swings land in the switch range. A noisy switch and a noisy swing give the same series and step
    share. Only IgM-size motions call swing.
- Separability: mid = share of node rows with |post_hi - post_lo| < 0.5 on the split axis.
  - qtoy2 0.12-0.27, OLE root 0.36, cdisc root 0.31-0.34: cuts reproduce.
  - OLE depth-1 0.53-0.78: cuts do not reproduce.
  - IgM continuous arms are also 0.13-0.29: mid measures separability, not discreteness.
- Poles: fixed minpop tails of the stage-B evidence difference (prob_state_evidence.txt).
  - Free refinement drifts continuum poles inward (OLE r2: -1.17 -> -0.75 sd; one pole grew 2.6k -> 10.4k).
  - Raw site-posterior tails overlap RNase P r1211's closed leaf only 0.63-0.78 (hard-gate risk).
  - Refined tails are purer on qtoy2 classes (closed 0.96-1.00, conf43 0.99) and about equal on IgM.
  - On IgM truth a node's pole pair is about as extreme as the best leaf below it (median 0.57-0.59 vs 0.56 of the oracle): not
    further, but without the slicing.
  - OLE minpop poles reproduce across seeds (Jaccard 0.46-0.70) where the cuts don't.
  - RNase P full692 r22 (graded, mid 0.92): a 5-10k open pole reads s1 -6.9 to -7.0 A own-pixel vs -5.3 for the most open leaf.
    The band-split held-out check shows that gap is IMPRINT: held out, pole -2.91 = leaf r221 -2.87 (random -2.60 +- 0.22).
  - OLE r2 replay (tree 6 rows, seeds 1/2): poles reproduce at Jaccard 0.91 (an upper bound, same rows). Own-pixel pole pA reads
    -0.67 on the r2_s4 axis vs saved leaves -0.60 / -0.53; not held-out.
**Why:** this is the user's switch-vs-swing request ([[discrete-splits-capture-the-span-states-that-matter]]). It corrects the
reading of [[at-node-switch-swing-call-2026-10-24]]: "switch" there means "blend-like".
**How to apply:** don't tune NB_SW. Decide delivery by map-call swing OR mid >= 0.5. Build poles from refined evidence at fixed
size.
- Native A/B (key swing_delivery, uncommitted in worktree swing-switch-extremes-58027e; cluster swing_switch_20261027/) was
  cancelled before any tree finished.
  - The rule acted at OLE r1/r2 and qtoy2 r211; every IgM node read switch (0.15-0.35); fallback nodes have no axis.
  - qtoy2 r211 poles were impure (holo 0.20; the other pole was a quality end), and the key-off refined cut isolated holo better
    (0.33).
  - Nulls: the swing key cannot hold the bar.
  - Before any rerun: a minority-vs-rest guard, and 2-CPU controllers (about 450 of the 1,670 core-h was idle controllers).
- 2026-10-27 night: RNase P r22 open pole FAILS the band-split held-out check (pole -2.74..-2.91 A vs leaf r221 -2.80..-2.93,
  random 5k -2.60 +- 0.22); its own-band -7 vs -5 A gap was imprint. qtoy2 r211 swing poles impure (one = best-quality end of
  the majority conformer). Native swing A/B trees cancelled unfinished; OLE r2 single-node replays (2 seeds) pending.
