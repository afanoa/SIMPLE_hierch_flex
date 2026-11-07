---
name: simple-run-parallelism-preference
description: "Default SIMPLE runs to nparts=10 nthr=8 whenever RAM/VRAM allows, instead of a single fat nthr"
metadata:
  type: feedback
---

Standing preference given 2026-08-17: run SIMPLE jobs as **`nparts=10 nthr=8`** rather than one
process with a large `nthr`, provided RAM/VRAM is not the binding constraint.

**Why:** 10x8=80 workers uses a big box (112 cores on
[[remote-build-recipe-fwl-c143766]]) far better than a single nthr=48 process, and `nparts=1`
degenerates some flex_pca code paths (see [[flex-pca-ab-testing-traps]]).

**How to apply:** pass `nparts=10 nthr=8` by default. Check the memory budget FIRST when the GPU
path is on — the handoff measured the banked M-step buffer at 9.4 GB per worker at 24 components
/ box 256, which OOM'd 4 workers on a 32 GB card. Ten workers sharing the 48.9 GB Blackwell card
leaves only ~4.9 GB each, so for GPU flex_pca either lower the component count, rely on the
byte-budget blocking (default 768 MB/worker), or cut nparts. RAM is not a concern on that box
(251 GB).

**Reinforced 2026-09-02 ("this is not how the code works"):** default EVERY
SIMPLE run to multiple parts x <=8 threads (nparts=10 nthr=8 on an 80-thread
box); 1 part x 48 threads is wrong — OpenMP saturates and the code is
architected for distributed parts (which also adds the window-striping
regularization benefit). The ONLY current exceptions, forced by code guards:
(a) paired runs with SIMPLE_COV_MERGE_POLISH (polish v1 is shared-memory only:
its distributed probe handoff is hardwired to the fit-A namespace — lift by
threading the polished namespace through save/load_probe_state), and
(b) legacy note about SIMPLE_ env validation under nparts did NOT bite for the
SIMPLE_COV_* set (production20 ran nparts=10 with the full env config fine).

**Hard rule after third violation (2026-09-02, "half the machine is idle
because of these inane settings"):** before launching ANY flex/SIMPLE job,
grep the launch command for nthr>8 without nparts — that combination is a
BUG, full stop. And VERIFY every remote script edit by printing the changed
line AND by checking the RUNNING process's cmdline (pgrep -af | grep -oE
"nparts=... nthr=...") — two remote in-place edits (sed -i, python heredoc
over ssh) failed SILENTLY this day and launched 48-thr jobs anyway. The only
reliable remote-edit path is: Write the file locally, scp it, verify on disk,
then verify the live process.

2026-09-24 update: this applies to flex2D/abinitio2D test runs too (synthetic suite and real 10076). I had used
nparts=8 nthr=24 (192 cores per real run) and 4 x 24 for synthetic, reached 640 running cores and was called
out. All suite scripts (synth2D/scripts/*.sh, 10076_flex2D/f200u/f200u.sh template) are now nparts=10 nthr=8.
Never run more than one real-data 10076 run at a time.
