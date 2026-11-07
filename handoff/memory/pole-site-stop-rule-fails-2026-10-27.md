---
name: pole-site-stop-rule-fails-2026-10-27
description: "Replay (11 saved trees + truth toys): 'continue splitting only if a site keeps a (novel) pole' FAILS as a stop rule - drops every 1.5-2k OLE leaf (0/15) and merges dozens of truth-supported IgM/qtoy splits; 'no kept pole' is a power limit"
metadata:
  type: project
---

2026-10-27, user asked "cant we use the pole or sites to figure out if we need to continue to split?"; the adaptive split/merge
session replayed it on saved trees (ole_12707/split_merge_replay/pole_stop/, code pole_stop.py; ~16 s jobs).
- Most nodes split via the fflex fallback, not a kept pole: OLE 6_ 5 pole / 32 fallback, 7_ 6 / 22, RNase P m33 20 / 99.
- P1 (a kept pole must exist): OLE 6_ 31 leaves -> 6 states, 0 of 15 small (1.5-2k) leaves kept; 7_ -> 7, 0/12.
- Truth: IgM m52M fallback splits 47 supported / 7 weak / 0 none (P1 would merge 50 truth-supported splits); qtoy2_bp P2
  (novel pole) merges 11 supported splits. Only qtoy2_a's fallback splits are weak.
- Novelty (|rho| < 0.5 vs ancestor axes) removes real nested sub-states of one motion.
**Why:** the pole rule's >= 3 SE halves test loses power at small nodes like the split test; absence of a pole is not absence
of structure. **How to apply:** don't propose pole/site presence as a stopping rule; poles stay for extremes and swings.
Related: [[discrete-splits-capture-the-span-states-that-matter]], [[gate-fakes-always-find-structure-2026-10-25]].
