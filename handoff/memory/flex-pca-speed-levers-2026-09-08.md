---
name: flex-pca-speed-levers-2026-09-08
description: "Measured 2026-09-08 anatomy of a flex_pca iteration on the 80-core box: the coupled K x K M-step INSERT is the CPU hotspot (37-47 of 56 s per 20k shmem iteration), the distributed slowness was the part reduce over beegfs (85-90 s), raw-plane residency is worthless (5 MB planes), and the insertion loop is barrier-bound so nothing else may run on the box during a timing arm"
metadata:
  node_type: memory
  type: project
  originSessionId: 7e382bce-970c-4cdb-bc32-4a1a1c1801a9
  modified: 2026-09-08T18:44:27.268Z
---

20k harness (10028 rows 1-20000, v3 paired recipe, neigs=20, 6 half-set iterations), fwl-c142339:

| arm | iteration | wall | where it goes |
|---|---:|---:|---|
| 10 parts x 8 thr, raw reads, parts on beegfs | 178 s | 1092 s | round 83 s + part reduce 85-90 s |
| 10 x 8, cache=yes, parts on local scratch | 28 s | 275 s | (numerics differed: see [[flex-pca-cache-taper-mismatch]]) |
| 1 process x 80 thr, raw | 56 s | 713 s | read 1.6, prep 9, project+solve 1.2, INSERT 37-47 |
| 1 x 80, raw, resident planes | 60-70 s | 572 s | prep 0 but copying 5 MB planes costs 7-10 s |

Full N (105k) production: paired stage 1312 s + states 820 s = 36 min at 10 x 8.

**Corrections to earlier beliefs.** (1) The E-step proper is ~1 s; the cost is the M-step insertion
into the K(K+1)/2 + K coupled accumulators (231 volumes at K=20). (2) "Shared memory beats
distributed" was the reduce reading ~1 GB parts per worker over beegfs; with parts on scratch the
distributed iteration is 2x faster than 80 threads in one process, as SIMPLE's architecture expects.
What is unlike the rest of SIMPLE is the part SIZE (K x K volumes vs an oritab + 2 volumes). (3) A
raw-prepped plane lives on the 720-pixel padded grid (5.3 MB); the cache serves it on the 128 grid
(~130 KB), so residency (`simple_flex_pca_planes.f90`, box branch flex-pca-resident) is tied to
cache=yes. (4) The insertion (`insert_planes_oversamp_coupled_batch_scaled`) parallelises over the
ROWS of one particle at a time with a barrier per colour class: ~22 rows of work for 80 threads and
100k+ barriers per pass. Any other process on the box (my own ssh log reads included) stalls it:
two iterations went 37 -> 161-181 s while I was reading logs. Never touch the box during a timing arm.

**Levers, in order.** (1) cache reads with the taper fixed; (2) parts always on local scratch,
cropped to the band (~3x smaller), threaded reduce; (3) re-decompose the insertion over particles
with thread-owned accumulation targets (up to ~3x on the paired stage, both paths); (4) box_rec 240
for the states stage; (5) `SIMPLE_COV_PROBE_MAX` capping the halves (user said NOT yet).
Related: [[flex-pca-distribution-changes-results]], [[flex-pca-conventions-branch-state]].

**2026-09-13 measured on cnga1 box 128 / lp 8 / 30 shells (job 53201259, csbdevel):** the master
tail is the wall-clock dominator and it is effectively SERIAL. Passive `/proc` sampling of the
master: 31 threads exist, **2-3 running, 29 sleeping**, and slurm's own accounting gives
AveCPU 19:31 against 20:16 elapsed = 0.96 cores. One probe-iteration master tail ran >9 min with
no log line while the distributed round over 10 parts took only 103.7 s. **system time is 26% of
master CPU** (utime 118103 / stime 41185 ticks), which points at allocation churn or page faults,
not arithmetic.
Consequences: (1) the MASTER THREAD BOOST is chasing the wrong thing -- the stage does not use the
threads it already has, so neither 16 nor 160 changes it, and the boost's write-through to
params%nthr is a live bug (see [[simple-script-yes-compenv-trap]]); (2) PCG solver work is NOT the
bottleneck at box 128, so warm starts / block CG / a mask-aware preconditioner all optimise a
minority of the wall; parallelising the coupled K x K insert in the master tail outranks them.
PCG conditioning numbers for when that is worth doing: box-64 M-step stops on maxits with
resid 0.527 after 8 iterations -> contraction 0.926/iter -> kappa ~= 680 -> rtol=1e-3 needs ~90
iterations, so 8 is ~9% of the requested solve. The preconditioner is `M(k)=rho(k)+floor(shell)`
(Fourier diagonal, sampling density, `build_precond` in simple_reconstructor_pcg.f90:1628) and is
**blind to the real-space mask**; the system solved is the projected `(P H P) u = P b` with P hard
(window>0, P^2=P), so ~91% of the spectrum is exactly zero at a 9% support.

**CORRECTION, later on 2026-09-13.** The paragraph above concluded "parallelising the coupled
K x K insert in the master tail outranks them". **That was wrong**, and so was my first reading of
the 9-minute log gap. Measuring instead of inferring gave:
1. The gap was NOT the pair-kernel finalize. It was the PCG solve: `PCG MSTEP it=1 even iters=8
   seconds=523.2` (odd 536.8, next fit 552.7) at box_crop=128, against **3.4 s at box_crop=64** --
   154x for an 8x volume.
2. The coupled K x K insert is already `!$omp parallel do collapse(3)`. Not the problem.
3. A real defect: `image%wthreads` defaults to `.false.` and `flex_pcg_t` built all three of its
   images through the 3-argument `image%new`, so **every FFT in the flex PCG path was
   single-threaded** (simple_reconstructor_pcg has always passed `fft_nthreads` through; see
   simple_reconstructor_pcg.f90:279 and :611). Fixing it changed NOTHING at box 64 (3.4 s) and a
   full M-step did not finish in 18 min at box 128 -- i.e. threading the transforms alone is
   neutral-to-worse. Do not assume this one is a win.
4. **The actual binding constraint is allocation churn.** `apply_operator` runs once per CG
   iteration and allocated ~2.8 GB per call: `cq` (ncomp x cdim complex = 681 MB at box 128 /
   rank 10), `acc` (68 MB), `work`, plus a `get_cmat` result (134 MB) and a
   `get_rmat`+`center_crop` pair (75 MB) PER COMPONENT. Measured with /proc sampling of the right
   pid: **1.05 of 16 cores busy, 31% of CPU in the kernel**, 3-6 threads flickering in R state.
   The solve was mmap/page-fault bound, not arithmetic bound.
**Method note that cost me a wrong conclusion:** `pgrep -f simple_exec | head -1` grabbed the
OTHER co-located job's master and gave "0.43 cores". Always disambiguate by
`readlink /proc/<pid>/cwd`.
**Also measured:** the sigma2 bootstrap subsample is NOT deterministic -- two runs of identical
configuration produced different `sigma2_state.bin` (md5 1e730fc0 vs 70aa994d) and noise scales
(3.361e-1 vs 3.346e-1, 0.4% apart) despite identical `SAMPLED E/O: 12500/12500` accounting. So no
flex_pca A/B can be bitwise; pin `sigma2_state.bin` if you need an exact comparison. Downstream
that is worth ~1% on eigenvalues, ~1% on GMM cluster membership, and up to 0.1 in |cos| on the
marginal components (1-3 stay identical to 3 dp). See [[flex-pca-ab-testing-traps]].
**Operational:** SIMPLE reused the `1_flex_pca` execution directory for a second run instead of
creating `2_flex_pca`, overwriting the first run's artifacts in place.
