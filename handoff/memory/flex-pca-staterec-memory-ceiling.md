---
name: flex-pca-staterec-memory-ceiling
description: flex_pca prints a machine-wide reconstructor memory estimate at startup; npreimages sets the state count that drives it, and exceeding RAM OOMs the whole box
metadata:
  type: project
---

flex_pca announces its own memory need in the first few log lines:

    >>> FLEX_PCA states=20 at reconstruction box 360 -> approx 23.73 GB of
        reconstructors per process x 11 processes = 260.99 GB machine-wide

**Read that line before walking away.** On 2026-08-28 a run printed 260.99 GB on
a 251 GB box (fwl-c142339), ran 35 minutes, finished the probe and all worker
partials, then was SIGKILLed (rc=137) during the master's merge. Available RAM
went 23 GB -> 239 GB in one 60 s sample. It also left 106 GB of partial
reconstructions behind on a filesystem that was already 98% full.

Cost scales as `states x box_rec^3` per process, times `nparts+1` processes. It
does NOT scale down with nparts — the master holds every state's padded grid at
once, so cutting parallelism does not rescue an over-large box_rec.

`npreimages` sets the state count: `nstates = max(MIN_NSTATES, params%npreimages)`
(simple_flex_pca_model.f90). It is a provisioning CEILING, not a target — the
merge stage can collapse states below it. Measured points on that box:

| states | box_rec | machine-wide |
|--------|---------|--------------|
| 6      | 128     | 3.90 GB      |
| 6      | 320     | 55.39 GB     |
| 20     | 360     | 260.99 GB (OOM) |

6 states is the usual working value. When queueing unattended runs, add a
watchdog that kills the job's process GROUP if available memory drops near zero
— see [[never-signal-simple-processes-by-name]].
