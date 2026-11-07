---
name: flex-pcg-affordability-breakthrough
description: "2026-09-14: the flex PCG M-step at box 128 was 520 s/solve and is now affordable -- a SCALED warm start from the gridding solution at maxits=1-2 beats 8 cold iterations on residual at 31-40% of the cost; the cost is (maxits+1) operator applications and the per-application cost is transform-bound"
metadata:
  type: project
---

**The measured ladder** (cnga1, 245,424 ptcls, box_crop=128, lp=8, ncomp=10, nthr=16, csbdevel,
one solve = even or odd of one M-step):

| config | iters | init resid | final resid | seconds |
| gridding (`maxits_pcg=0`, **unmasked basis**) | 0 | - | - | **0.3** |
| warm+scaled | 1 | 0.634 | 0.534 | 161.3 |
| warm+scaled | 2 | 0.634 | 0.524 | 205.4 |
| warm unscaled | 8 | 0.959 | 0.474 | 514.4 |
| cold (the old default) | 8 | 1.000 | 0.594 | 518.5 |

**ONE scaled-warm iteration beats EIGHT cold ones** (0.534 vs 0.594) at 31% of the cost. Extra
iterations have sharply diminishing returns (0.534 / 0.524 / 0.474 at 1 / 2 / 8).

**What made it work.** `solve` already computed the gridding (block-Jacobi) solution and then threw
it away with `x = 0.0`, calling a cold start "the production convention". It is not: BOTH entry
points of simple_reconstructor_pcg warm start via `window_div(x)`, with three guards learned the
hard way -- a division floor `PCG_SUPPORT_DIV_MIN=0.1`, a free initial-residual rejection
(`start_max_rel_resid`), and a curvature cold restart. Flex had the cold-restart machinery only.
Then: `cg_core` COMPUTED `start_scale = <b,Bx0>/<Bx0,Bx0>` and never used it. Measured
`corr(b,Bx0)=0.7736`, `scale=0.5186` -- the gridding iterate points the right way and is ~2x too
large. Applying the scale is free (hp = B x0 is already in hand) and takes the initial residual
from 0.959 to exactly `sqrt(1-corr^2)=0.6337` (measured init 0.6337).

**Cost model:** a solve costs (maxits + 1) operator applications, ~60-80 s each at box 128, because
the start residual needs one. For a COLD start that application computes B*0 = 0 and is pure waste
(~11% of a cold solve); production documents the shortcut, flex now does it too.

**Per-application cost is transform-bound, and the inventory was bigger than it looks:** per CG
iteration, 2*ncomp transforms on the boxpd=2*box lattice PLUS 4*ncomp on the native lattice,
because `bandlimit` does two transforms and is called twice per component; `apply_precond` adds
2*ncomp more. All were serial through ONE shared `self%wimg`/`self%nimg` while OpenMP was applied
only to inner voxel loops. Fix = per-thread image pool over COMPONENTS, copying
`image%construct_thread_safe_tmp_imgs` (simple_image_core.f90:126), which builds every pooled image
with `wthreads=.false.`.
**DO NOT thread inside a transform.** Passing `wthreads=.true.` was measured a net harm: all four
gdb stack samples landed in `do_futex_wait -> __new_sem_wait_slow -> os_sem_down` in
libfftw3f_threads (160 transforms per M-step x a 16-thread pool = semaphore bound). Reverted.

**Dead ends, measured, do not retry:** threading inside FFTW (above); removing 2.8 GB/call of
allocation churn from `apply_operator` (real, bitwise-neutral, kernel time 26%->3%, but <10% of wall
because the solve was never allocation-limited); the gridding/PCG hybrid schedule where only the
delivered basis is mask-constrained -- **rejected on science**, the probe iterations are where the
basis forms, so an unconstrained M-step lets solvent/micelle variance into the components and then
into the E-step (cf. [[flex-pca-solvent-mask-dominates-rna-basis]]).
**`maxits_pcg=0` ships an UNMASKED basis:** `put_back` does `set_rmat(x*self%env)`, the envelope not
the support, so the 0.3 s figure is a floor, not a usable option. A masked cheap path is not a
one-liner because the half-set FSC that sets the Wiener filter is taken on UNMASKED halves.

Related: [[flex-pca-speed-levers-2026-09-08]], [[cnga1-has-reproducible-continuum]],
[[flex-pcg-cold-start-and-fsc-prior]], [[simple-script-yes-compenv-trap]].

**2026-09-14 VALIDATED ON SCIENCE, not just speed.** Rank 10, box 128, lp 8, same data, original
build vs the optimised one (pool + scaled warm start + all transform loops pooled). Cross-half
matched |cos| requiring BOTH halves >= 0.7:
  original:  1 .995  2 .913  3 .872  4 .894  5 .777  6 .690  7 .664  8 .644  9 .732  10 .615 -> 6 pass
  optimised: 1 .989  2 .945  3 .757  4 .780  5 .740  6 .813  7 .852  8 .678  9 .843  10 .581 -> 8 pass
Components 6 and 7 went from failing to passing (.690->.813, .664->.852). The prediction held: a
better-converged M-step lifts the MARGINAL components, which were solve-limited, and leaves the
leading ones alone. Individual numbers scatter both ways because sigma2 is redrawn per run.

**CORRECTED COST BREAKDOWN (this is the number to plan from), rank 10 / maxits_pcg=8 / nparts=10:**
| stage | per run | share of fit |
| part reduce (`paired_reduce_parts_v5`, 4 calls) | 788 s | 47% |
| distributed rounds (5) | 500 s | 30% |
| PCG solves (18 @ ~22 s, = 9 operator applications each) | 396 s | 24% |
Per operator application 2.4 s. maxits=2 gives 7.5 s/solve, maxits=8 gives 22 s/solve.
**The part reduce is now the dominant stage and was never measured before.** It is
`paired_reduce_parts_v5` in simple_flex_pca_em_estep.f90: a SERIAL loop over parts, one part
resident at a time, folding ~116 GB of full-size accumulators off beegfs per reduce at ~0.59 GB/s.
It also allocates/frees ~37 GB of read buffers per reduce inside `fold_probe_part_v5_fit`
(cbuf/rbuf/rcbuf; rcbuf alone is npairs x the expanded lattice, ~1.87 GB, twenty times per reduce)
-- the same churn pattern as the solver had, NOT yet fixed, and two routines share the pattern.
**MORE PARTS IS HARMFUL, not helpful:** each part writes the FULL accumulator set, so reduce I/O is
nparts x accumulator size. At 20 parts the round halves (~52 s) but the reduce doubles (~394 s):
446 s/iteration vs 300 at 10 parts. This is what the recorded "shmem beats distributed was the
beegfs part reduce" note meant. nparts=10 is near optimal.
**Overlap is NOT a local change:** the round's wait lives in simple_qsys_ctrl's scheduling loop
(jobs_done array, returns only when all parts finish), so by the time the reduce runs every file
already exists. Overlapping needs a submit/wait split in the controller -- the layer where a bug
hangs every future run on a file that never arrives.

---

**PER-COMMIT ATTRIBUTION (2026-09-14), box_crop=128, rank 10, cnga1 — measured at IDENTICAL budget:**

| build | commit added | iters=2 (warm) | iters=8 (cold) |
| pre-optimization | -- | 205.4 s | 518.5 s |
| opt (`af106f7c9`+`4bacc96f3`) | per-thread image pool; single-threaded FFTW plans | ~58.8 s | ~73 s |
| opt7 | preconditioner transforms pooled | 55.3 s | 70.6 s |
| **opt8 (`3dbe7f86e`)** | **finalize, finalize_rhs and the shell ridge pooled** | **~7.6 s** | 21.2 s |
| opt9 (`d1cac887b`) | single-threaded Y plans + voxel-outer coupled multiply | -- | 18.8 s |
| opt10 (`8e081f9d9`) | threaded master-side part fold | -- | 19.5 s |

**The last pooling commit was the biggest win (7.3x), not the first (3.5x).** `finalize`,
`finalize_rhs` and the shell ridge were SERIAL transform loops that run once per solve regardless
of the iteration budget -- i.e. they were the fixed cost, independently visible in the warm ladder
as ~111 s of setup inside a 161 s single-iteration solve (161.3 s at iters=1, 205.4 at 2, 512.8 at
8 -> ~50 s marginal per iteration on ~111 s fixed). **Lesson: on this solver, profile fixed
per-solve cost before per-iteration cost.**

**The honest overall factor is ~27x, and it reproduces at two budgets** (205.4 -> 7.6 at iters=2;
518.5 -> ~19-21 at iters=8). **Do NOT quote "69x"** -- that divided an 8-iteration cold solve by a
2-iteration optimised one and smuggled a 4x budget cut into a code-speedup claim.

**Numerically inert, as required:** residual at iters=2 is 0.5238 / 0.5229 / 0.5231 across
pre-opt / opt7 / opt8. Warm start's gain is in the INITIAL residual (0.634 vs 1.000 cold), and
convergence is non-monotone -- 0.634 -> 0.534 (1 it) -> 0.524 (2) -> 0.382 (8) -- so iterations
1->2 buy ~0.010 while 2->8 buy 0.142; that argues against the very short budgets independently of
cost. The part reduce is the remaining dominant stage: [[flex-pca-part-reduce-is-io-bound]].
