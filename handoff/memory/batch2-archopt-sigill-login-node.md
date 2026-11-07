---
name: batch2-archopt-sigill-login-node
description: "On batch2.ncifcrf.gov a USE_ARCHOPT build compiled on a csbdevel compute node SIGILLs on the login node: build on Ice Lake (x6342), login is Skylake 6150, so every simple_exec that touches vectorized code must go through srun"
metadata:
  type: reference
---

Verified 2026-09-14. `USE_ARCHOPT:BOOL=ON` is the default in the batch2 CMakeCache, so `-march=native`
bakes in the *build node's* ISA.

- login/xfer node `fsitgl-xfer03p` = Intel Xeon Gold **6150** (Skylake-SP)
- `csbdevel` is **heterogeneous**: `x6342` (Ice Lake, 48 cores, 256 GB) and `x6740P`
  (Granite Rapids, 96 cores, 2 TB)

A build that lands on an Ice Lake node runs fine on every csbdevel node (Granite Rapids is a
superset) but dies on the Skylake login node with `Program received signal SIGILL: Illegal
instruction.` — one backtrace per OpenMP thread, no other diagnostic.

Light paths (`print_project_info`, `info_image`, `script=yes` generation, `update_project`)
happen to survive on the login node; anything numeric (e.g. `automask2D`) does not. Do not read
a login-node success as proof the binary is portable.

- `norm` (2026-10-25 sinfo) holds 32 x6150 Skylake nodes (36 cores) beside 58 x6342 and 7 x6740P. A tree whose step jobs may land on
  norm (SIMPLE_QSYS_PARTITION=csbdevel,norm) must export `SBATCH_CONSTRAINT='x6342|x6740P'` in its launch script, before
  simple_exec, so every sbatch it calls avoids them (stree_port/stree_f692t.sbatch does).

**How to apply:** wrap interactive SIMPLE calls on batch2 in
`srun -p csbdevel -c 8 --mem=16G -t 20 bash -c "cd <dir> && simple_exec ..."`.
Note `--pty=false` is rejected (`must be numeric file descriptor`) — just omit it.
Real work goes through `script=yes` + `sbatch` anyway, which always lands on a compute node.
Related: [[simple-script-yes-slurm-submission]], [[simple-args-generator-lib64-trap]].
