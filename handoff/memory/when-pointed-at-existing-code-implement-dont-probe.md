---
name: when-pointed-at-existing-code-implement-dont-probe
description: "when the user points at an existing SIMPLE pattern (\"look at the other code\", \"a cache like elsewhere\"), read it once and implement; no rounds of file-count/format/size probing jobs"
metadata:
  node_type: memory
  type: feedback
  originSessionId: b8bc4431-3743-4b17-aafd-6d31e354ea19
  modified: 2026-10-09T17:08:17.026Z
---

2026-10-09, cls_expansion reader race: after the user said "do a cache like we do elsewhere" I kept probing (reader history,
both caches, stack counts, MRC modes, free space via extra sbatch jobs) and the user interrupted twice: "why are you taking
so long". **Why:** the direction was already decided; the probing added nothing the implementation needed.
**How to apply:** read the referenced pattern once (API + build/read routine), say the plan in two lines, write the code,
build, run one test. Sizes/timings come from the test run itself. Related: [[do-the-requested-job-not-diagnostic-sprawl]],
[[act-on-the-literal-ask]].
