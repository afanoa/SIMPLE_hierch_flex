---
name: box-build-needs-gcc16-ld-library-path
description: "On fwl-c142339 a `make install` fails at 93% with `print_ui_json failed with exit code 1` unless LD_LIBRARY_PATH holds the gcc-16.1.0 lib64: the UI-JSON validation step RUNS the freshly built simple_exec, which needs GLIBCXX_3.4.32/GFORTRAN_10; the Fortran itself had compiled fine"
metadata:
  node_type: memory
  type: project
  originSessionId: 7e382bce-970c-4cdb-bc32-4a1a1c1801a9
  modified: 2026-09-08T18:45:11.599Z
---

Every build/run script on the box needs BOTH before `make`:
`export PATH=/mnt/nasapps/development/modules_elmlund/gcc/16.1.0/bin:$PATH
LD_LIBRARY_PATH=/mnt/nasapps/development/modules_elmlund/gcc/16.1.0/lib64`. Symptom without the
second: make_*.log ends with `RuntimeError: print_ui_json failed with exit code 1` from
production/validate_complete_ui_json.py, and a few lines above `simple_exec: /lib64/libstdc++.so.6:
version GLIBCXX_3.4.29 not found`. Cost me one queue cycle on 2026-09-08 (queue9).
Also: a fresh full build at -j32 takes ~4 min there; an incremental flex rebuild ~2 min.
Related: [[flex-pca-em-split-and-syntax-check-recipe]], [[remote-compute-box]].
