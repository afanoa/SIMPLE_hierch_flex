---
name: numerically-close-not-bitwise
description: "STANDING (user 2026-10-24, 'again'): don't chase bit-identical results; numerically close is adequate for validating changes (speed levers, node_proj, subtree jobs, determinism)"
metadata:
  node_type: memory
  type: feedback
  originSessionId: 25d47be4-7d1d-48d8-abdf-97b1c3c946ab
  modified: 2026-10-24T05:27:43.265Z
---

User, 2026-10-24 about 01:20, while the overnight work was verifying flex determinism and same-seed identity byte by byte: "again i
dont care about being bit identical. numerically close is adequate". The "again" means it was said before.

**Why:** reproducibility that matters is the science (the same states, close maps), not the last float bits. Bitwise checks cost core-h
and dev time, and they turned node_proj's float-level noise (maps 0.9998) into a false "failure" that blocked a default the user had
already chosen.

**How to apply:**
- Validate result-safe changes by numerical closeness: maps corr >= about 0.999, the same tree shape and decisions, labels close (ARI
  at about the level of two same-seed runs).
- Don't gate defaults on byte equality, and don't spend runs or code to make outputs bitwise.
- Make flex determinism or ordered sums optional niceties, never blockers: flex_pca at nthr > 1 is not bitwise reproducible
  (dynamic-schedule sums), so a bitwise bar only measures float noise. The patch .claude/patches/flex_determinism_on_056895a5c.patch
  is for the user's review only. About 160 core-h of bitwise verification was cancelled on this rule.
- Never queue bitwise-identity jobs; compare against the same-seed run-to-run floor instead.
- Related: [[landscape-must-be-reproducible-and-splits-data-derived]] (reproducibility = the same states across runs, not bits),
  [[flex-pca-ab-needs-simple-seed]], [[hetero-tree-4p5m-scaling-2026-10-23]].
