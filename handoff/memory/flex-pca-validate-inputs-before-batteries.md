---
name: flex-pca-validate-inputs-before-batteries
description: Gate every flex_pca campaign on a cheap input-quality check (K=1 pair) before committing multi-hour batteries; prefer one generous-K CV pair over per-K ladders
metadata: 
  node_type: memory
  type: feedback
  originSessionId: fb1d9e33-0137-47d8-8d28-e3a32330f274
  modified: 2026-09-01T15:18:10.079Z
---

On 2026-09-01 a full K-ladder battery + T5 was launched on SIMPLE-refined 10028
poses before checking pose quality. K=1 cross-half cosine came back 0.146 (vs
0.965 on cryoDRGN consensus poses) — every queued arm on those poses was dead on
arrival, hours were burned, and the user called it out twice ("wasting time",
"be more principled").

**Why:** flex_pca results are extremely pose-sensitive (even a mild 4-iter
re-refinement from good poses degrades every rung). Any battery on a new
project/pose set is worthless until the input is validated.

**How to apply:**
1. Before ANY multi-run campaign on a new pose set: run ONE K=1 (or small-K)
   half pair (~25 min) and check cross-half cosine against the known-good
   reference (~0.96 on 10028 cryoDRGN). Below ~0.8 → stop, fix poses first.
2. Rank selection: never queue a per-K ladder. One paired job at generous
   neigs (10) gives the whole CV curve (per-component matched cosines) in one
   fit-pair; k = #components above floor. (User-prescribed protocol.)
3. When the user says a better input exists (cryoDRGN poses), switch ALL queued
   work to it immediately — don't let doomed arms keep draining the box.
4. Cheap decisive experiment first, expensive confirmation second — always
   state the cost and what each run decides before launching.
5. Remote orchestration: ONE flat sequential script with idempotent stages
   (DONE markers); on any plan change kill everything, edit, relaunch (costs
   seconds). Never layer background waiters polling logs for sentinel lines —
   they add poll latency, pkill -f self-match kills, and silent stalls when a
   producer dies (all three happened 2026-09-01; user: "stupid computing
   controls").
6. When new code lands, RE-EXAMINE the queued jobs: a step queued for an old
   reason (rec3D to feed the project-FSC filter) becomes superfluous the moment
   the reason is superseded (per-state eo-FSC filter) — prune it, don't run it.
7. Monitors must cover CRASH and HANG, not just chain sentinels: watch the
   job's own log for ERROR/Backtrace AND process liveness. A crashed worker +
   hanging master writes no sentinel — silence looked like progress for 50 min
   (2026-09-01 rec3D objfun=euclid crash; user: "you just wasted 50 mins").
   Also: imported/bare projects need objfun=cc — euclid requires sigma2 files.
8. EVERY job monitor uses the full template, no exceptions: (a) error patterns
   (ERROR|THROW|Error termination) in the job's run.log AND in
   SIMPLE_SUBPROC_OUTPUT (worker crashes appear ONLY there while the master
   hangs silently), (b) stall detection (log sizes frozen across 3 checks with
   processes alive = hang alert), (c) process liveness vs completion sentinel.
   A sentinel-only monitor missed a worker crash + master hang AGAIN on
   2026-09-02 (distributed polish round-1 missing-basis crash) after the same
   lesson was recorded the day before. Recording a lesson is not applying it —
   reuse the template verbatim.
