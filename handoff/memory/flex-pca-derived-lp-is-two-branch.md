---
name: flex-pca-derived-lp-is-two-branch
description: "Derived-lp is a TWO-BRANCH rule (run one rung finer: chases the band = compositional, pins = that's your lp), and derived-neigs is an information BOUND, not a setter"
metadata:
  type: project
---

Measured across the 2026-08-18/19 representation round (10076 + 10049). The single-rule
derived-`lp` was REFUTED and salvaged as a branch test.

**derived lp** = heterogeneity resolution: mean even/odd FSC of the top components, first
crossing below 0.143. It must come from the HETEROGENEITY FSC, never the map FSC — 10049
pins at 17.8-19.4 A while its consensus map is 3.58 A, a 4-5x gap.

**The branch test:** run one rung finer and watch the crossing.
- crossing CHASES the band (10076: crossing ~= band edge at lp=20/16/12) => band-filling
  compositional data; the rule is circular there and `lp` must come from the science ladder
  (measured interior optimum 16 A @ box_crop=80).
- crossing PINS (10049: 17.8-17.9 A at both lp=20 and lp=16; 19.35 A interior at lp=16/bc=60
  with DEFLATE_BG) => band-limited continuum data, and that pin IS the derived lp.

**derived neigs** = count of components with in-band FSC > 0.143 — but on both datasets this
came back RIGHT-CENSORED (40/40 on 10076, 32/32 and 40/40 on 10049). So it is an information
BOUND, not a setter: the delivered optimum sits BELOW it exactly while the state stage is the
binding constraint (see [[flex-pca-supervised-ceiling-instrument]]). Auto-neigs cannot be a
pure information rule until the state stage stops binding.

**Why:** this is the honest status of the knob-elimination program in
[[flex-pca-must-be-knob-free]] — lp/neigs are now DERIVED (a measurement exists) but not yet
AUTOMATIC (the measurement does not uniquely determine the setting). Treat any claim that
flex_pca is knob-free as premature.

**How to apply:** on a new dataset, run two band rungs and read the BAND/RANK line before
committing to `lp`. Do not read "top components" by index on the EM arm — the EM basis is NOT
variance-sorted (components 39/40 carried the 2nd/3rd largest sd on 10076), so rank by in-band
mean FSC instead. Operating point as of 2026-08-20: neigs=40, lp=16 @ box_crop=80 on 10076;
lp=20 discovery band on 10049.
