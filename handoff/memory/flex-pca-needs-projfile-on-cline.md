---
name: flex-pca-needs-projfile-on-cline
description: "distributed flex_pca without projfile= dies with \"mean-scale band mismatch\" (em_mean.f90:569) — master runs at the native box, workers at the derived box_crop"
metadata:
  node_type: memory
  type: project
  originSessionId: 846f1d9d-8f6f-4d4a-b48d-b023b38c86e4
  modified: 2026-09-23T17:50:58.252Z
---

Reproduced 2026-09-23 on origin/master 622d0acb4 (IgG, nparts=4 nthr=20, fwl-c142339).
`apply_flex_pca_defaults` runs BEFORE params%new and derives box_crop/lp/box_rec (and the
project consensus vol1) only when `projfile` is ON THE COMMAND LINE. Without it the master
skips all of that and params%new falls back to box_crop=box (native). params%new then
auto-discovers the cwd *.simple and sets projfile on the cline, so the worker job_descr carries
projfile but no box_crop/lp — each worker re-derives box_crop from the default lp=16 -> nyq
disagrees with the master's mean-scale file -> THROW. Code dates from 8b9922b5c (2026-09-17);
not a regression from later pulls.

**Why:** the user hit this with `simple_exec prg=flex_pca ... nparts=4` run inside a project
dir with no projfile=; with projfile= the same run passes (workers get box_crop=64, nyq=32).

**How to apply:** always pass projfile= to flex_pca. Even shared-memory runs without it
silently fit at the native box (no crash). Fix direction (not applied, user said no changes):
resolve the cwd project before the derivations, or stamp the master's box_crop/lp/box_rec into
job_descr. Local-qsys runs on the box also need `SIMPLE_PATH` exported or qsys_env%new throws.
Related: [[flex-pca-ab-testing-traps]], [[remote-compute-box]].
