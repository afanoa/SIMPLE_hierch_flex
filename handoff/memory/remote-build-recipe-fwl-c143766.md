---
name: remote-build-recipe-fwl-c143766
description: "fwl-c143766 is the Blackwell GPU box; the flex_pca CUDA build needs cuda/13.0 (NOT 12.8) and arch 120, with a split Fortran/CUDA compiler pair"
metadata:
  type: reference
---

`ssh atifao@fwl-c143766.ncifcrf.gov` — 112 cores, 251 GB RAM, **NVIDIA RTX PRO 5000 Blackwell,
compute capability 12.0, 48.9 GB VRAM, driver 580.105.08**. Distinct from
[[remote-compute-box]] (fwl-c142339, 80 cores, Quadro RTX 5000 / sm_75). `/home` is local XFS per
box, NOT NFS-shared, so ssh keys must be installed on each host separately.

Working GPU build of the `flex-polar-gpu` line (verified 2026-08-17, BUILD_EXIT=0):

```
module load cuda/13.0 cmake/3.31.0
GCC16=/mnt/nasapps/development/modules_elmlund/gcc/16.1.0   # Fortran
GCC14=/mnt/nasapps/production/gcc/14.3.0                    # C/C++ and nvcc host
export LD_LIBRARY_PATH=$GCC16/lib64:$GCC14/lib64:$LD_LIBRARY_PATH
cmake -S . -B build-gpu -DCMAKE_BUILD_TYPE=Release -DUSE_FLEX_CUDA=ON \
  -DCMAKE_CUDA_ARCHITECTURES=120 -DCMAKE_Fortran_COMPILER=$GCC16/bin/gfortran \
  -DCMAKE_C_COMPILER=$GCC14/bin/gcc -DCMAKE_CXX_COMPILER=$GCC14/bin/g++ \
  -DCMAKE_CUDA_HOST_COMPILER=$GCC14/bin/g++ -DUSE_OPENMP=ON
```

Three traps, each of which fails the build outright:

1. **cuda/12.8 does not work.** `simple_flex_gpu_kernels.cu` calls `cublasSetEmulationStrategy` /
   `CUBLAS_FP64_EMULATED_FIXEDPOINT_MATH` (the FP64-emulation path for the fourth-power term) with
   NO version guard. That API first appears in cuBLAS 12.9.
2. **Architecture is 120, not the 75 in the 2026-08-17 handoff** — that handoff was built for a
   Turing card. The repo default in `cmake/Dependencies.cmake` is already 120.
3. **nvcc refuses any host compiler above GCC 14** (cuda/12.8 ceiling 14, cuda/13.0 ceiling 15),
   but SIMPLE's Fortran needs gfortran 16.1.0 — hence the split pair above.

Binaries land in `build-gpu/production/`, NOT `build-gpu/bin/`. Everything is user-level: module
loads plus a build dir in `$HOME`, no sudo. Note `/mnt/nasapps/development/modules_elmlund` IS
group-writable — never install there, it is shared with the whole lab.

**Distributed runs need `bin/` populated.** The build puts binaries in `build-gpu/production/`,
but the distributed launcher looks for `$SIMPLE_PATH/bin/simple_private_exec`. With `bin/` empty,
every worker dies instantly with "No such file or directory" and the PARENT KEEPS RUNNING — the
job looks alive in `ps` and makes zero progress, and the message does not contain the string
"ERROR!" so it slips past log greps. Fix once per build:

```
ln -sf $SIMPLE_PATH/production/simple_private_exec $SIMPLE_PATH/bin/
```
(same for simple_exec, single_exec, simple_stream, simple_test_exec, simple_persistent_worker).
`nparts=1` masks this entirely because it runs in-process.
