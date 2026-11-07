---
name: flex-cuda-build-option
description: flex GPU was REMOVED from master 2026-09-28 (uncommitted in ~/SIMPLE_master); older branches still carry USE_FLEX_CUDA — arch must precede enable_language(CUDA)
metadata:
  type: project
---

**2026-09-28: the whole flex GPU implementation was cut on master** (worktree
`/Users/atifao/SIMPLE_master`, left uncommitted for the user): the .cu, `simple_flex_gpu.f90`,
`compile_flex_gpu.sh`, the `USE_FLEX_CUDA` option/CMake block, the `flex_gpu` test entry, and every
GPU branch in flex (em_iter banked adjoint/forward, fused E-step, polar device E-step, rec3D device
insert, latent_ops `flex_dev_prep_hook`, `align_halfplane_inplane`). The OpenMP-offload GPU stack
(`compile_gpu.sh`, `simple_gpu_utils`, `USE_OPENMP_OFFLOAD`) is separate and was kept. Flex branches
(flex2D etc.) still have the GPU code and need the same cut when merged. VERIFIED CPU-neutral 2026-09-30: fast
gate 13/13 (job 53399766) and a seeded nthr=1 flex_pca A/B vs unmodified master, all 60 outputs byte-identical
(job 53488470, [[flex-pca-ab-needs-simple-seed]]); batch2 trees ~/SIMPLE_nogpu (patched) and ~/SIMPLE_nogpu_base.
COMMITTED 2026-09-30 as 6642aaefa 'removing flex gpu kernels and compile instructions' on local master in
~/SIMPLE_master, rebased onto origin/master 32457e631, PUSHED to origin/master 2026-09-30 (32457e631..6642aaefa). compile_clean.sh + compile_debug.sh both rc=0 with
gate 13/13 on the pre-amend d189f28fb (jobs 53489971/53489972); the amend only dropped the unused em_iter
`nll_mix_add` (single-file recompile clean). User: no full rebuilds for one-line amends.
The notes below describe the pre-removal state.

The flex GPU kernels are gated by one CMake option, `USE_FLEX_CUDA` (default
OFF). It is independent of `USE_OPENMP_OFFLOAD` / `USE_OPENACC` — verified
2026-08-28 by configure matrix: with any of those ON but `USE_FLEX_CUDA` unset,
nvcc is absent from the cache, the `.cu` is absent from the source list, and the
linked binary has zero CUDA symbols.

Design: all `#ifdef USE_FLEX_CUDA` live in `src/main/flex/simple_flex_gpu.f90`
only; the `#ifdef` wraps the `bind(c)` interface block while module state is
declared unconditionally. Callers carry no guards — they branch at runtime on
`flex_gpu_available()`, which is a literal `.false.` without the define. So a
CPU build compiles the facade and every entry point is a `THROW_HARD` stub.

**The trap:** since CMake 3.18 `enable_language(CUDA)` itself defines
`CMAKE_CUDA_ARCHITECTURES` (to nvcc's default, 52). So the common idiom

    enable_language(CUDA)
    if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)   # never fires
        set(CMAKE_CUDA_ARCHITECTURES <x>)

is dead code. Set the architecture BEFORE `enable_language`. This is why the
tree's original `set(CMAKE_CUDA_ARCHITECTURES 120)` never took effect — builds
were silently getting 52, and 75 only ever came from an explicit `-D`.

**Not Blackwell:** fwl-c142339 has a Quadro RTX 5000, compute capability 7.5.
The arch-120 default in the original code was written for a different box; see
[[remote-build-recipe-fwl-c143766]]. Box build recipe: cuda/12.8 nvcc,
gcc-14.3.0 g++ as CUDA host compiler, gfortran 16.1.0, `-DCMAKE_CUDA_ARCHITECTURES=75`.
