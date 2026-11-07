---
name: flex-pca-em-iterations-overfit-cross-half
description: "2026-09-14: on cnga1 the probe EM is NOT converged at 4 iterations (loglik still falling steeply) but cross-half update agreement DECLINES monotonically from iteration 2 at every rank -- likelihood and reproducibility pull in opposite directions, so more iterations is not automatically better"
metadata:
  type: project
---

**The trajectories are in every run log already and I had read those logs without reading them.**
Per-iteration lines: `PROBE ITER n -2logL/N`, `update agreement (even|odd, prev basis deflated)`,
`mean principal-angle cosine vs previous basis`.

**cnga1 `sweep_n10`** (box 128, lp 8, rank 10, n_probe_iters=4, maxits_pcg=8):

| iter | -2logL/N | delta | principal-angle cos vs prev | agreement | n>=0.9 |
| 1 | 17.34 | -- | -- | -1.000 | 0 |
| 2 | 9.85 | -7.49 | 0.9616 | 8.679 | 2 |
| 3 | 5.65 | -4.20 | 0.9766 | 8.212 | 1 |
| 4 | 3.79 | -1.86 | 0.9867 | **7.542** | 1 |

**Two real effects in tension.** (a) UNDER-ITERATED: the likelihood is still dropping steeply at
iteration 4 (deltas -7.5, -4.2, -1.9, roughly halving, so several more iterations of real movement
remain) and the basis is still rotating (cos 0.987, not ~0.999). (b) OVERFITTING: cross-half
`update agreement` falls monotonically after iteration 2 at EVERY rank -- n10 8.679/8.212/7.542,
n16 13.637/13.019/11.785, n24 20.086/19.161/17.438 -- and `n>=0.9` peaks at 2 on iteration 2.
So the half-fits track their own noise more closely as EM proceeds.

**Consequence: "use 15-25 iterations because the winning recipes did" is NOT safe to assume here.**
Those recipes (ribosome 25 iters -> visible 40S rotation; PfCRT injection 15 iters + Wiener off ->
recovered a 16 deg motion end to end) were higher SNR. cnga1 may sit below that threshold, in which
case the optimum is FEW iterations with early stopping on the agreement measure -- which is exactly
what `SIMPLE_COV_PROBE_CONV` / `EO_PATIENCE` were for and what the recipes deliberately disable
(see [[flex-pca-em-speedup-plan]]).

**Open, needs the endpoint measure:** per-iteration bases are NOT saved (only `pcit000_*`) and
`flex_pca_paired.txt` holds only the final merge, so reproducibility-vs-iteration cannot be
recovered from finished runs -- it needs an explicit iteration ladder (4 / 8 / 20) scored on the
endpoint matched |cos|. Baseline to beat: 8 of 10 at mean |cos| 0.798, K=3.
Related: [[cnga1-is-a-continuum-not-discrete-states]], [[flex-pca-probe-iters-and-whitening]].

---

**RESOLVED 2026-09-14 by a 6-arm ladder on cnga1** (2/3/4/8/14/20 probe iterations, everything else
identical: maxits_pcg=8, rtol=0, Wiener ON, rank 10, box_crop 128, lp 8, refine3D envelope mask,
opt12 binary which reproduces the opt8 anchor digit-for-digit end to end).

| n_probe_iters | >=0.7 both | **>=0.85 both** | mean min\|cos\| | loglik/particle | K |
| 2 | 7 | **5** | 0.778 | -34.011 | 4 |
| 3 | 6 | 3 | 0.754 | -32.094 | 3 |
| **4** | **8** | 3 | **0.798** | -30.487 | 3 |
| 8 | 6 | 1 | 0.716 | -26.011 | 4 |
| 14 | 5 | 1 | 0.700 | -21.818 | 7 |
| 20 | (pending) | | | -19.296 | 9 |

**1. n_probe_iters=4 is the OPTIMUM on the primary metric.** The prediction that 2-3 would beat it
(from the declining in-run agreement) was WRONG. The default was already right.
**2. n_probe_iters=2 gives the most STRONGLY reproducible components (5 at >=0.85 vs 3).** Early
stopping keeps the dominant modes crisp (0.95-0.97) while never resolving the weak ones (comp 3 at
0.41/0.46); 4 iterations spreads quality evenly instead. Prefer 2 when you want a few crisp modes.
**3. MORE ITERATIONS ACTIVELY HARM low-SNR data.** 8/14 give 6 then 5 components, ONE above 0.85.
The "use 15-25 because the ribosome and PfCRT recipes did" advice was WRONG here -- those were higher
SNR. Do not generalise iteration count across datasets.
**4. LIKELIHOOD IS NOT A USABLE MODEL-SELECTION CRITERION.** It improves MONOTONICALLY by ~15 nats
from 2 to 20 iterations while reproducibility strictly degrades, and it prefers the WORST basis on
the ladder. A convergence check on loglik would happily run to 20 and return 5 of 10 components.
**5. K tracks fit flexibility, not structure.** On this continuum (dip 0.00000 at d=10 and d=16) K
goes 3 -> 4 -> 7 -> 9 as the fit converges: each extra cluster is a boundary drawn through
continuous density. Same pattern wherever the fit degraded (box 96, lp 6, cyl45 mask all gave K=4).

**Caveat on the in-run diagnostic:** `update agreement (even|odd, prev basis deflated)` declines
monotonically (8.679 / 8.212 / 7.542 / 6.829 at iters 2-5, 5.318 by iter 11) but it is referenced to
the PREVIOUS DEFLATED BASIS, so it partly measures how much the basis is still moving. Its decline
DECELERATES (~0.7/iter to iter 5, ~0.25/iter to iter 11), so a linear extrapolation from the first
few points overstates the collapse -- I made that error. Rank arms on the ENDPOINT cross-half matched
\|cos\| from flex_pca_paired.txt, never on the in-run agreement.
Related: [[cnga1-is-the-apo-closed-dataset]], [[flex-pca-probe-iters-and-whitening]].
