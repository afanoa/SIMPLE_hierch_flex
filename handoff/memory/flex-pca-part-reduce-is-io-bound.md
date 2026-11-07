---
name: flex-pca-part-reduce-is-io-bound
description: "2026-09-14: the flex_pca master part reduce is I/O-bound, not compute-bound -- ~16 GB of STREAM reads at ~86 MB/s; threading the elementwise folds bought only 4.6% because the arithmetic was never the cost"
metadata:
  type: reference
---

**Measured anatomy of one paired iteration** (cnga1, 245,424 ptcls, box_crop 128, lp 8, rank 10,
nparts=10, nthr=16, batch2/beegfs; from `sweep_fold10/fold10.53202057.out`):

| stage | seconds | share |
| worker qsys round | 100-112 | 24 % |
| **master part reduce** (`reduced v5 paired`) | **185-196** | **45 %** |
| 4 x PCG kernel fold | 4 x 11 = 44 | 10 % |
| 4 x PCG M-step solve | 4 x 19.6 = 78 | 18 % |
| total `PAIRED ITER` | 427-453 | |

`paired_reduce_parts_v5` (`simple_flex_pca_em_estep.f90:566`) is a strictly serial
read -> add -> read -> add chain: `access='STREAM'`, one part at a time, and every
`!$omp parallel workshare` add waits for its blocking read. Bytes per reduce at rank 10
~= 10 parts x 2 fits x (10 comps x (2 x complex + 2 x real) + 2 x npairs=55 real) on a
~1.2 M-voxel crop grid **~= 16 GB, i.e. ~86 MB/s** -- under a tenth of what beegfs can do.

**So threading the folds was optimising the wrong term**: opt10 (`8e081f9d9`) vs opt8 moved the
reduce only 196.8 -> 187.8 s mean (4.6 %, n=4 each, p~0.06). The instrumented PCG kernel fold is
**11 s**, not 190.
It is also NOT bit-identical (max rel 9e-3 in eigenvalues, median 8e-6) -- but the science is
unchanged: 8 of 10 reproducible both ways, mean |cos| 0.796 vs 0.798, same K=3 with the same
occupancies to 0.03 % (22977/158379/64068 vs 64155/158367/22902, labels permuted).

**The real levers, in order:** (1) the accumulators are ~95 % zeros -- only shells <= kstop=30 of
64 are populated at lp 8, so shipping a band-limited payload is a ~20x byte cut and *exactly*
bit-identical; (2) read the P parts concurrently into P buffers and combine in fixed part order --
bit-identical to serial, costs P x 1.6 GB of RAM; (3) the write side emits ONE record per
component but the read side does FOUR reads of it, so 4x the syscalls on that payload.

**Trap hit:** arms are pinned to a build tree by `SIMPLE_PATH` in their own `run.sh`, but WORKERS
take the binary from the projfile's compenv -- `sweep_n24`'s workers ran opt8 while its master ran
opt8 too, and `mk_sweep.sh` does `cp -f $B/proj.simple`, so re-running it would overwrite a
projfile that live workers are reading. Fresh directories only.
Related: [[flex-pca-speed-levers-2026-09-08]], [[flex-pca-ab-testing-traps]].

---

**CORRECTED 2026-09-14 (measured on a real part file; the earlier numbers here were wrong).**

A probe part file at box_crop=128 / rank 10 / nparts=10 is **11.83 GB**, so the reduce moves
**~118 GB per iteration, not 16 GB**, at **622 MB/s, not 86 MB/s**. 622 MB/s is a respectable
beegfs single-stream rate, so the reduce is slow because of VOLUME, not inefficiency -- which kills
the "read the parts in parallel" lever (there is no 10x bandwidth gap to reclaim).

**Where the bytes actually are** (header dims: ncomp=10, lattice 70x133x133, npairs=55, kmix=16):

| payload | lattice | size | share |
| cmat_e/o + rho_ex/ox, 10 comps x 2 fits | crop 1.24 M | 0.58 GB | 5 % |
| rho_e/o, 55 pairs x 2 fits | crop | 1.06 GB | 9 % |
| **kpk_e/o + rpk_e/o (PCG kernels)** | **doubled 8.45 M** | **~10.1 GB** | **~86 %** |

So packing cmat/rho alone buys ~14 % of the bytes (190 s -> 165 s, worthless). **The PCG kernel
payload is the whole prize.**

**THE TRAP: nonzero FRACTION is not bounding-box VOLUME.** Both payloads are ~93 % zeros
(crop 6.51 %, kernels 6.89 % nonzero) but the doubled lattice is in **physical FFT addressing**
(`self%wrap` puts real mass on the opposite face), so an axis-aligned bounding box around that
content covers **~50 % of the lattice**. A bounding-box implementation therefore yields only
**2.2x blended**, not the 15x the sparsity suggests. I built the box version first on exactly this
confusion and had to redo it. An explicit **index list** of populated lattice points is
addressing-agnostic (centred ball, wrapped octants, anything) and recovers the full ~15x.

**Bit-identity is structural, not asserted:** build the list as the UNION of the nonzeros of exactly
the arrays shipped, so every omitted voxel is an exact zero and folding it is a no-op. Accumulators
start at exact +0.0 and deposits are gated on `h^2+k^2 <= nyq_eff(nyq_eff+1)`.
Must also validate the LATTICE DIMS on read (`subhdr(2:4)` / `header(4:6)` were dropped) -- a part
from another lattice otherwise yields in-range indices at wrong addresses = silently misplaced mass.
`write_probe_part_kernels`/`fold_probe_part_kernels` are SHARED by the v4/v8 and v5 paths, so both
version constants must move together (now PROBE_PART_VERSION=11, PROBE_PART_VERSION5=10).
The box scan / list build runs on the parallel WORKERS, off the reduce's critical path.

**Per-iteration budget at maxits_pcg=3** (worker round 105 s, reduce 190 s, M-step 32 s, fold 44 s
= 371 s): band-box -> 242 s; index-list -> ~157 s, i.e. 20 iterations in ~52 min vs 2.1 h.
**There is no particle-subsample knob in flex_pca** (no nsample/update_frac/ptcl cap), so the 105 s
worker round cannot be cut for cheap exploration without adding one.
