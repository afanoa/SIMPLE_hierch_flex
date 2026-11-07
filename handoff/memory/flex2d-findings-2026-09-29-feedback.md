---
name: flex2d-deliver-results-first
description: when the user asks to see results, give the fastest working output first and explain before any digging; no multi-step root-cause chases
metadata:
  type: feedback
---
On 2026-09-29 the user asked for hemoglobin class averages and a viewable directory; I spent ~30 min on wrong
diagnoses (header stats, averager change) and they rejected three tool calls ("wild goose chase", "what are you even doing").
**Why:** they need to look at results to decide; a broken output with three explanations is worse than a quick workaround.
**How to apply:** when output is broken, say in two lines what is wrong and the quickest way to get a viewable result,
get a yes, produce it, THEN fix the code. Verify a diagnosis before acting on it (e.g., test e2display in the browser mode
the user uses) and never present a guess as the cause. Related: [[flex2d-findings-2026-09-21]].
