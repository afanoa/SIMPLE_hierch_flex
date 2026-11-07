# Stage-B reference shrinkage toward relatives (prior-treeshrink): verdict

2026-11-06/07, side session of the "Tree as full PPCA mixture" orchestrator. Merged build
`/mnt/beegfs/elmlund/afan/hetero_tree_20261028/merge_master2/install`; gate fixtures under `math_reworks/gates/`.

## Verdict

WORSE than production at every dose tried (half, tree_wiener's full Wiener, strong one-part): each loses the RNase P closed
split on at least one of the two real trees, and none changes the toys or the null. Not worth wiring into the stage loop in Fortran.

Shrinking each state's deviation from the node mean per shell (tree_wiener's emp2 cross-half Wiener, applied between three
single-iteration stage-B calls) does not change the toys' truth readouts. It removes the closed-state split on both real RNase P
trees: the closed corner is still in the stage-B maps, but the node's split merges it, so the tree would deliver the closed state
inside one leaf (hard gate fails). The cause is that the shrinkage
suppresses exactly the band the split tests (14.6-10.2 A, shells 16..23 at box 88): stage B stops partitioning on it, the
states' maps differ half as much there, and the split loses its evidence. The loss is monotone in the dose and no dose gains
anything that replicates.

## Numbers (the control reproduces production: .997-1.000 of labels shared, same split, on every node)

Doses: copy = control (no shrinkage), dev2h = half dose, dev2 = tree_wiener's emp2 estimator (full Wiener), dev1 = one part
(no focus mask, the strongest in the split's band). "moved" = labels differing from production.

RNase P closed-state nodes (closed corner by the gate rule: closed arm >= .95, open arm <= .6 on a stage-B map):

| node | arm | moved | closed-row enrichment (best state) | s02 map closed/open | split: test-band separation (threshold) / stage-B states vs EM groups (nats, SE) / effect floor P/N | children |
|---|---|---|---|---|---|---|
| mm2_r131 r2 (seed 2) | production | - | 1.112 | 1.032/.519 corner | .134 (.065) / +4.4 (.7) / 2.74 | 2 |
| | copy | .2 % | 1.112 | 1.032/.520 corner | .133 (.070) / +4.0 (.7) / 2.73 | 2 |
| | dev2h | 6.6 % | 1.159 | 1.040/.518 corner | .090 (.041) / +2.0 (.4) / 2.49 | 2 |
| | dev2 | 13.6 % | 1.217 | 1.047/.524 corner | .072 (.047) / +0.3 (.3) -> EM 17984/3182 / .657 merged | **1** |
| | dev1 | 11.1 % | 1.159 | 1.037/.519 corner | .073 (.049) / -0.2 (.4) -> EM 17473/3693 / .619 merged | **1** |
| mm2_r131_s1 r2 (seed 1) | production | - | 1.209 | .997/.544 corner | (production split 2) | 2 |
| | copy | .3 % | 1.207 | .997/.544 corner | .114 (.072) / +4.0 (.7) / 2.96 | 2 |
| | dev2h | 6.6 % | 1.213 | 1.002/.550 corner | .101 (.069) / +0.4 (.5) -> EM 15881/4553 / .945 merged | **1** |
| | dev2 | 12.4 % | 1.198 | 1.003/.556 corner | .053 (.071) FAIL | **1** |
| | dev1 | 10.2 % | 1.199 | .995/.547 corner | .050 (.084) FAIL | **1** |

Toys and null (best-group F1 for closed / holo; split = state_split children):

| node | arm | moved | stage-B NMI / ARI | closed F1 (rec, pur) | holo F1 | split children |
|---|---|---|---|---|---|---|
| mm2_qtoyflex r | production = copy | .1 % | .168 / .139 | .17 (.73, .10) | .62 | 3 = states |
| | dev2h | 2.7 % | .168 / .139 | .17 (.72, .09) | .62 | 3, EM regrouped: closed .16, holo .64 |
| | dev2 | 5.3 % | .166 / .137 | .17 (.71, .09) | .62 | 3, EM regrouped: closed .16, holo .64 |
| | dev1 | 5.3 % | .164 / .136 | .17 (.71, .09) | .61 | 3, EM regrouped: closed .16, holo .64 |
| mm2_qtoy2 r | production = copy | .0 % | .656 / .536 | .55 (.99, .38) | .94 | 5 = states |
| | dev2h | .8 % | .652 / .532 | .55 (.99, .38) | .94 | 5 = states |
| | dev2 | 1.5 % | .648 / .528 | .56 (.99, .39) | .94 | 5 = states |
| | dev1 | 1.4 % | .652 / .532 | .55 (.99, .38) | .94 | 5 = states |
| mm2_qtoynull r | production = copy | .2 % | - | - | - | 1 |
| | dev2h | 5.4 % | - | - | - | 1 |
| | dev2 | 11.2 % | - | - | - | 1 |
| | dev1 | 11.3 % | - | - | - | 1 |

Every dose leaves the toys within +-.02 F1 / .008 NMI of production and the null one leaf; on qtoyflex any dose flips the
split from the stage-B states to its EM regrouping (same truth readout).

Per-band effect of the shrinkage (shrink_it<k>.log; kept power of the deviation inside / outside the focus mask, dev2):
>20 A .99 / .3-.7, 20-14 A .82-.88 / .02-.13, 14-10 A .13-.46 / ~0, <10 A ~0. Focus mask 1.4-4.8 % of the sphere on the real
and toy nodes, 22 % on qtoynull, whose deviation cross-half FSC above 20 A is .98 (a null split looks like signal to it).

## What ran

- `shrinkB.py SRC DST MODE --counts --mskdiam`: reads the K state maps a refine3D iteration ships (recvol_stateNN{,_even,_odd}.mrc,
  the ML-regularised pair the next iteration aligns against), writes them under the same names to DST (all sidecars copied).
  copy = unchanged (control); dev2 = each half's deviation from the count-weighted node mean, per Fourier shell, times
  W = C (C + A diag(s2_half) A^T)^+ (C = PSD cross-half covariance of the deviations, A = I - 1 w^T), in two parts inside /
  outside tree_wiener's 20 A between-state focus mask, as emp2; dev2h = dev2 at half dose; dev1 = one part (no focus mask).
  The node mean is untouched; merged = (even + odd) / 2.
- `shrinkB_run.sbatch` (MODE, LIST): per node, the node's own stageB/hetero_step.sh refine3D call with maxits=1, run three times
  in one directory with startit=1,2,3 (maxits is run length, so call k is iteration k with production's seed, keyed draws,
  committed sigma2 and fsc files); refs passed through shrinkB.py MODE between calls as vol1..volK; then the node's
  split/hetero_step.sh state_split under SIMPLE_SPLIT_CONTRAST=maps SIMPLE_TREE_KEEPMAPS=yes; readout f2r_ab_readout.py vs
  production (+ copy) and, on r131, make_view.py + arm_readout.py on the stage-B maps.
- `shrinkB_chain.sbatch` runs two such passes in one allocation; `shrinkB_readout.sbatch` re-reads all arms side by side;
  `shrinkB_test.sbatch` validated shrinkB.py on production maps (copy byte-identical, node mean kept to 2e-7).
- Jobs: 55010199 (test), 55010293 copy / 55010294 dev2 (four gate nodes), 55011976 copy on r131_s1 r2 then dev2h,
  55011977 dev2 on r131_s1 r2 then dev1 (five nodes). ~80 cores at a time. All finished (00:17); no jobs open.
- Runs: `/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/shrinkB/runs/<key>_<node>_<mode>/` with stageB/ (step<k>.log,
  cmd<k>.txt, refs_it<k>/ = the maps passed to call k+1, shrink_it<k>.log = kept power and deviation FSC per band, b_ptcl3D.txt),
  split/ (step.log, split_stats.txt, split_children.txt), prod/ (production's dump), readout/ (labels.txt, arm.txt).
  Logs `shrinkB/logs/`.

## How to resume

On batch2 (login node: sbatch only): `cd /mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/shrinkB`, then
`MODE=<copy|dev1|dev2|dev2h> LIST=$PWD/<list> sbatch --export=ALL -J shrinkB_<mode> shrinkB_run.sbatch` (a run directory that
exists is skipped; move it aside to redo it), and `sbatch shrinkB_readout.sbatch` for the side-by-side readout. Lists:
shrinkB_list.txt (the four gate nodes), shrinkB_list_s1.txt (mm2_r131_s1 r2), shrinkB_list_dose.txt (all five). A new dose is
one more `choices` entry in shrinkB.py. Each node takes ~5 min at 40 cores (4 workers x 8 + the 8-core driver).

## Why not pursue it further

- No gain anywhere that replicates: the toys move <= .02 in F1 / .008 in NMI at every dose; r131's closed-row enrichment moves within the
  build/seed spread (~.07, the f2r session's production vs argmin replay on r131 r1): +.105 on seed 2, -.011 on s1 at full dose.
- The harm is structural: the multivariate Wiener removes the deviation where its cross-half SNR is low, which at a node with a
  subtle real difference (the closed corner) is the 14-10 A band the split uses to decide. A dose small enough to keep the split
  (half dose keeps it on seed 2 but not on s1) keeps the references close to production and buys nothing.
- On the null the deviation's low-resolution cross-half FSC is .98: the shrinkage keeps a null split as signal, so it cannot
  act as a noise guard either (the nulls stay one leaf because the split's test rejects them, as in production).
- Already regularised per state: in pcg_focus mode the ML regularisation of each state takes its prior from the FSC of the
  state's deviation from the consensus background (simple_rec3D_pcg_strategy.f90, dev_even/dev_odd), so the shipped references
  are filtered by the deviation's own SNR before this step; the multivariate shrinkage adds the cross-state coupling on top.
