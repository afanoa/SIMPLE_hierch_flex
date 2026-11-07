---
name: prototype-in-python-before-fortran
description: new flex2D model ideas get a Python toy prototype that shows a clear win before any Fortran change or batch2 build
metadata:
  node_type: memory
  type: feedback
  originSessionId: 95245f2d-90ec-4676-978e-28a5533c3f54
  modified: 2026-10-03T22:10:37.490Z
---

2026-10-03, after my FLEX2D_CLS_MARG tilt-marginal went straight into the E-step (two commits, two batch2 builds) and lost on
the hbK15k toy in every variant: "You can still work on the out of plane tilt, but i would try to get this to work on a toy
first and with python scripts before commiting to a large build".

**Why:** a Fortran change costs a build, a snapshot and distributed toy runs per iteration of the idea, and an idea that has
not been shown to work in principle burns those cycles.

**How to apply:** for a new modelling idea (tilt subspaces, two-level partitions, new class evidences), first write a
self-contained numpy/scipy toy with known truth (synthetic projections of the toy maps, controlled views/states/SNR), compare
against the mean-only baseline and an oracle, and only port to flex2D once the toy shows a clear, explainable gain. Related:
[[flex2d-findings-2026-09-21]] (ep).
