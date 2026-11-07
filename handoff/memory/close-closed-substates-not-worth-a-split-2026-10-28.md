---
name: close-closed-substates-not-worth-a-split-2026-10-28
description: "user judged recipe14's r111-r1 closed sub-state pair (1.04 vs 0.97 closed-arm, 14.8k/10.7k) 'pretty close, not really swinging, but indeed closed' -> statistically real != worth a state; resolution-aware merge being replayed"
metadata:
  type: feedback
---

2026-10-28 ~19:40, after viewing ~/Downloads/rnasep_r111r1_closedA_s02_14840.mrc and closedB_s01_10688.mrc (10 A, unsharpened):
"those states are pretty close. Not really swinging but they are indeed closed."

**Why:** the pair passed every statistical test (re-test z 8.4, rho 0.049, localisation L 0.78) and the equivalence merge-up keeps it
apart (cannot-link), yet the user sees no meaningful motion; splitting 25.5k closed particles into 14.8k + 10.7k costs resolution.
Consistent with [[discrete-splits-capture-the-span-states-that-matter]] (extremes at good resolution; more splits lose resolution).
**How to apply:** don't present statistically distinct close sub-states as wins; prefer a resolution-aware merge (difference must be
resolvable at the resolution the split costs) - replay queued with adaptive-split-merge; the closed STATE itself is the gate
([[rnasep-closed-state-must-be-retained]]), not a closed sub-pair.
- Later the same evening the user added: "those two closed states differ in one of the other helices" - a real, local second
  variation inside closed (diff_centroid: one compact blob 78 A from the closed arm, more density in s01). Deliver-as-two vs merge is
  then a resolution trade-off, offered as a depth-3 vs depth-4 cap choice, not a verdict of over-division.
