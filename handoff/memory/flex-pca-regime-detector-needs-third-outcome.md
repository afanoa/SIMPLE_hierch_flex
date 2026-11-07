---
name: flex-pca-regime-detector-needs-third-outcome
description: "The merge pair-ratio regime detector has no 'insufficient signal' verdict; EMPIAR-10049 gives the calibration point — ratios above 1.0 and a negative median"
metadata:
  type: project
---

The handoff's regime detector reads the merge's disattenuated pair-ratio distribution and separates
discrete from continuum:

| dataset | pairs | min | median | max |
|---|---|---|---|---|
| Ribosembly (discrete) | 105 | 0.293 | 0.670 | 0.999 |
| 10076 (discrete) | 120 | 0.518 | 0.790 | 0.993 |
| IgG (continuum) | 171 | 0.815 | 0.878 | 0.963 |
| **10049 RAG (2026-08-17)** | 120 | **-1.022** | **-0.930** | **1.026** |

RAG is off-scale against all three: every calibration set is positive and bounded below 1.0, while
RAG has a negative median and a max ABOVE 1.0. The merge said why in its own words — nine pairs
"indistinguishable within their own noise" (ratios 0.989-1.026) — and collapsed 16 -> 6 states.

**Why:** a ratio above 1.0 means cross-state correlation exceeds the within-state noise floor, i.e.
the states are not separable at their own noise. That is a THIRD regime the detector was never
calibrated for: not discrete, not continuum, but **insufficient signal**. Without it the pipeline
silently delivers 6 or 16 states that mean nothing — exactly what happened here.

**How to apply:** for the knob-free design ([[flex-pca-must-be-knob-free]]), gate on this BEFORE
choosing a readout: max ratio >= 1.0 or median < 0 => report "no separable heterogeneity" instead of
delivering states. Also note `state_axis=1` is NOT a continuum readout — it changes target placement
but still emits DISCRETE-STATE maps, so regime-switching needs a trajectory output that does not yet
exist. See [[empiar-10049-rag-flex-pca-null]].

## Second detector added 2026-08-19: deviation-FSC vs band (label-free)

A distinct, cheaper instrument from the same round: compare the state-deviation even/odd FSC
crossing against the assignment band.

- deviation-FSC **>>** band => COMPOSITIONAL (10076: crosses ~10 A against a 16 A band)
- deviation-FSC **~=** band => CONTINUUM (10049: 21.5 A against a 20 A band)

Mechanism: a missing subunit puts full structural detail into the state difference, far beyond
the band where *assignment* stops improving; continuum flexibility yields neighborhood-average
states whose differences die at the band. Needs no labels. Runs alongside merge contraction as
the second regime signature (10076 contracts toward truth; 10049 merges nothing).

**The third outcome is STILL OPEN as of 2026-08-20.** Neither detector reports "insufficient
signal" — and note the pair-ratio calibration above came from the pre-native-pose 10049 run whose
null is now superseded ([[empiar-10049-rag-flex-pca-null]]), so that off-scale row may be measuring
pose error rather than a genuine third regime. Recalibrate on the native-pose project before
building the gate on it.
