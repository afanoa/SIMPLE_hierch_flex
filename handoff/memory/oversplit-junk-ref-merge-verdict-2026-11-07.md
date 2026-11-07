---
name: oversplit-junk-ref-merge-verdict-2026-11-07
description: "VERDICT lever 4 (over-split + junk reference + merge, stage B at fixed poses): junk ref works, over-split recovers states K lacked (qtoy2 root .674 -> .803 with cross-half re-assignment) but not close pairs (r1 fixed point, r131 closed); the map-level cross-half difference merge test fails the null by self-selection"
metadata:
  node_type: memory
  type: project
  originSessionId: 5482e66b-b4b2-4b69-b534-a7b318bf3360
  modified: 2026-11-07T17:30:02.875Z
---

Lever 4 at stage B (fixed stage-B poses, state_eval build), 2026-11-07, handoff/oversplit/VERDICT.md:
- A data-only junk reference (robust-z sum of -corr and -contrast > 3) works wherever junk exists (junk F1 .14 -> .73 r1, .23 -> .58 root) at no cost;
  on real data (r131) it grows into a low-contrast sink (1.2k -> 7.7k).
- Over-splitting (2K seeds from the node's deconvolved flex mixture, agglomerated closest-first by Bhattacharyya) recovers the
  states stage B's K lacked (qtoy2 root class accuracy .674 -> .785 after 7 passes, .803 with cross-half re-assignment; closed F1
  .55 -> .95; qtoyflex closed .17 -> .45) but never the close pairs: r1 stays at its 80/20 fixed point (.757 vs .756), 43|49 never
  separate, r131's closed corner is neither found nor lost. Seeds are the limit (seed NMI .53 root, .21 r1, .004 r131); K+2 seeds
  fuse the big classes and end below stage B.
- TRAP: a map-level cross-half difference test (even-half diff vs odd-half diff of a state pair, against label shuffles within
  view x defocus strata) keeps every pair everywhere and FAILS THE NULL: five states of one homogeneous structure (qtoynull, s1
  constant) score x .10-.15 (z 8-13) with view overlaps .9. Cause: self-selection (a particle's label comes from its own image, so
  its noise is correlated with its label in both halves); shuffled labels cannot mimic it. Tilt-band splits are a second blind
  spot (an artificial tilt split of one class scores x .34, z 7). One cross-half re-assignment inherits the signature (null agreement
  .60 vs chance .20, x rises to .41); only a half-independent EM with state matching across halves (the tree's merge-up design)
  can pass the null.

**Why:** the parent asked whether over-split-then-merge with a junk reference fixes stage B's reference limit; the gain is real
but the merge step as built cannot be trusted, and the null is the proof.

**How to apply:** never judge a data-driven split by a cross-half test of the resulting maps unless the halves were labelled
independently; port over-split + junk reference only with the merge-up's half-independent rule; do not retry the flex mixture as
a seed for close pairs. Related: [[stageb-oracle-references-limit-2026-11-07]], [[flex-to-refinement-coupling-verdicts-2026-11-06]],
[[null-acceptance-majority-leaf-2026-10-27]].
