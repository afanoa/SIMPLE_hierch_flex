---
name: batch2-new-source-needs-cmake-rerun
description: adding a NEW .f90 under src/main needs `cmake .` in the build dir before make (the main/ GLOB_RECURSE has no CONFIGURE_DEPENDS); and a background build's done-file must have a fresh name per build, or the waiter reads the previous build's exit code
metadata:
  type: feedback
---

2026-09-20: `simple_flex_pca_varimax.f90` was added; `make -j 24 install` in `~/SIMPLE_region2/build` failed
at the first `use simple_flex_pca_varimax` ("Cannot open module file") because the source list is a
configure-time `file(GLOB_RECURSE SIMPLE_src "main/*.[Ff][90][08]")` (src/CMakeLists.txt:81). Only the UI
defaults glob carries CONFIGURE_DEPENDS. My waiter also reported BUILD_EXIT=0 from the previous build's
`make.done`.

**Why:** the glob is evaluated when cmake configures; make never re-globs.

**How to apply:** for a new file run `cmake . > cmake.log 2>&1 && make -j 24 install > make.log 2>&1` in
the build dir (the cache keeps Release/ARCHOPT/OpenMP/prefix); name the done-file per build
(`make2.done`) and confirm with `ls -la bin/simple_exec` (timestamp) plus `grep -c "Error [0-9]" make.log`.
Related: [[build-on-batch2-login-node]], [[flex-region-pca-branch-state]].
