---
name: coarse-retest-band-leaks-in-em-band-2026-10-29
description: "TRAP: moving state_split's re-test into the cross-half EM band (40-20 A, where real splits live) leaks: all 8 null roots split at T 0.92-0.99, z 15-29; a coarse band needs a held-out (fit one half, score the other) design"
metadata:
  type: project
---

2026-10-29, adaptive split/merge session (route9/band, band_task.sh; doc section 10.19).
- The cross-half EM labels each half's particles from their own data, with the other half's maps as references, over shells
  <= 20 A. Inside that band the groups carry selected noise that correlates across halves.
- A 40-20 A re-test (lp=20 + negative SIMPLE_SPLIT_KTGAP) splits every qtoy2-matched and OLE-matched null (T near 1, while the
  shuffle null stays 0.02-0.24). The "power" on m33 r2121212 (z 0.98 -> 16.8) and the closed pair (3.66 -> 11.7) is the same leak.
- The production re-test (about 15-10 A) is outside the EM band; that is why the earlier leak test was clean.
**How to apply:** never test a split in a band the EM labelled from; for a coarse band, fit the K groups on one half and score
held-out per-particle log-likelihood gain on the other (proposed, not built).
Related: [[contrast-switch-kills-real-rnasep-splits-2026-10-27]], [[null-split-diffuse-real-split-localized-2026-10-28]].
