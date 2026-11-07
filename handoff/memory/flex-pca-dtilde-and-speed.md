---
name: flex-pca-dtilde-and-speed
description: "d_tilde WAS the dominant cost; on the packed+CG path the solve is now 0.2 s and d_tilde governs memory, not time (remeasured 2026-08-16). Tree default is 128, not 64. No analytic criterion can pick it, and why"
metadata: 
  node_type: memory
  type: project
  originSessionId: e31f34db-3e71-4d09-ab25-71f10f48fa4c
  modified: 2026-08-06T19:42:40.205Z
---

**SUPERSEDED IN PART, 2026-08-16.** Remeasured on IgG 20k / 11 cores / shared memory / current
defaults, release build: the REDUCED SOLVE is **0.2 s including packing** at d_tilde=128
(`packed+CG, 0.545 GB`), i.e. <1% of a 144.6 s run. The solve term is constant in N, so the small
particle count does not explain it -- the packed-svec + matrix-free CG path collapsed the constant.
d_tilde now governs MEMORY, not time. The "48% of the run" framing below is history.
Also: `COV_DEFAULT_DTILDE` in the tree is **128**, not the 64 this note claims.
Where the time actually goes now: probe EM 40.8 s (28%, 2 of 5 iters after early stop at cos 0.979),
state reconstruction 33.2 s (23%), orthonormalize 12.5 s, column accumulation 12.9 s. Probe EM is
the LARGEST stage and beats state reconstruction 1.23x even when it stops early; at the full 5 iters
it would be ~3.1x. See doc/for_developers/flex_pca_method.md "Where the time goes".

Original note, State as of 2026-08-06, IgG 100k, nparts=10 nthr=8, box_rec=128. **893 s -> 237 s** and better
state separation, from two independent things:

- **1.32x** from the Aug-6 optimizations, measured cleanly at matched d_tilde=177 (561.8 -> 425.7 s).
  Per stage: SNR 25.1 -> 8.0 s, state reconstruction 100.2 -> 68.0 s, trajectory 51.5 -> 39.5 s.
- **The rest from d_tilde.** The reduced solve costs ~[d(d+1)/2]^2, so d is the most expensive
  number in the program. It had been set by the MEMORY budget, which made it depend on solver path
  and machine (packed -> 250, dense -> 177) rather than on the specimen.

## The d_tilde curve (spread = mean pair corr of the most-spread state triple, LOWER is better)

    d=250  893 s  16.3 GB  0.628      d=90   242 s  2.24 GB  0.594
    d=177  426 s   5.25 GB 0.624      d=64   237 s  2.24 GB  0.542  <-- best
    d=128  279 s   2.27 GB 0.605      d=45   225 s  2.24 GB  0.587
The 2026-08-05 handoff volumes score 0.574, so d=64 separates states better than the reference run.
`COV_DEFAULT_DTILDE` is now **64**. Calibrated on IgG ONLY -- confirm on Ribosembly, where GT can
count recovered states (d=128 gave 14/16 there).

**CAVEAT added 2026-08-06 (audit).** This curve was scored on `spread` ALONE, and
[[flex-pca-ab-testing-traps]] trap 4 says spread must never be read without `cmin` because a
degenerate blob decorrelates from everything and wins. The curve is monotone in d
(250:0.628 -> 64:0.542), which is what a blur metric looks like; and it ranks the run the user
considers GOOD (0.574) BELOW the setting it recommends. Every independent criterion asked for more
rank, including cross-halfset generalization at 118 -- the one statistic measured to track signal
(+0.60, [[flex-pca-basis-rank-over-specified]]). **The d=64 default is not established; re-score
the sweep with cmin before trusting it.** The timings above are unaffected. Note also that this
default has NEVER EXECUTED -- it lives in an uncommitted diff that was never copied into
`build/bin/` (see [[flex-pca-ab-testing-traps]] trap 1).

## No analytic criterion finds it -- all four missed, in both directions

    energy floor (1e-6 relative)          256   4x too high
    memory budget                     250/177   machine-dependent
    spectrum noise-bulk (cov_signal_rank)  86   closest, still high
    cross-halfset generalization v1         6   10x too low  (tested capture vs ABSOLUTE chance;
                                                cap(q) decays with eval(q), so low-ENERGY
                                                directions were rejected for being low-energy)
    cross-halfset generalization v2       118   118 of 128 column directions reproduce

v2 is the important negative: **reproducible is not conformational.** Alignment residual, CTF and
ice reproduce across halfsets at cos ~ 1 and contribute nothing to state separation. The optimum is
fixed by a DOWNSTREAM objective (state separation) that column-space statistics cannot see -- the
same wall as [[flex-pca-signal-ranking-is-blocked]], one level lower. Both estimators remain in the
tree REPORT-ONLY (`SIMPLE_COV_DSIGNAL=1` to bind the first).

**Therefore a constant, biased LOW.** The loss is asymmetric: below the optimum costs a little
quality and no time; above costs quality AND time as d^4. Bias low is the safe automatic policy.

## Krylov is dead as a priority

The solve was 426 s / 48% of the run at d=250. At d=64 it is ~19 s. Even an infinitely fast solver
now buys <1.2x end-to-end, and it would change the numerics. The handoff ranked randomized block
Krylov top of the not-started list *because* the solve looked like 48% -- that was an artifact of
the bad default, not a property of the algorithm.

## Where the time is now, and what scales

At 100k the profile is flat (reconstruction 68, trajectory 40, states_combined_eo 44, solve ~19).
Stages capped by the subsample are CONSTANT in N; the embedding and all reconstruction stages are
O(N). Extrapolating to 3M: reconstruction ~83% of runtime, solve ~8%. So at scale the insertion
path is the target, not the solve and not polar (polar accelerates PROJECTION; the 3M cost is
INSERTION into a Cartesian grid). Measured internals: insertion is 94% of the reconstruction stage,
`mean_live_states` = 5.7-6.2 of 20, so tightening the GMM responsibility floor is a direct linear
lever. See [[flex-pca-ab-testing-traps]] before running any comparison.

## lp is likely over-set too (open, sweeps running)

Per-shell column FSC on IgG dies between shell 14 and 15 (~27 A), but the default lp=15 sets
kstop=25 -- about eleven shells of pure noise in the covariance band. That FSC is computed from the
RAW halfset accumulators before Wiener regularization, so it is not circular. Unlike d_tilde this
one should be automatable: "does this shell carry reproducible signal" is what FSC answers, and the
code already logs it every run. lp also drives box_crop through auto_box_crop (safety 1.25), so
lp=22 -> bc=44 and lp=30 -> bc=32, i.e. 3.1x and 8x fewer voxels.
