# Stage-B reference shrinkage toward relatives (shrinkB): status

FINAL (00:20, 2026-11-07): WORSE than production at every dose (half / full emp2 / strong one-part): each loses the RNase P
closed split at node r2 on at least one of the two real trees (full and strong on both; half on s1), toys within +-.02 F1 /
.008 NMI, nulls one leaf. Not worth Fortran. Handoff: /mnt/beegfs/elmlund/afan/hetero_tree_20261028/handoff/prior-treeshrink/
(VERDICT.md with all tables, NOTES.md, scripts); repo copy untracked in worktree confident-hodgkin-9ba093 at
scripts/hetero_tree_handoff/prior-treeshrink/ (commit blocked by the session's permission policy). No jobs open.

Owner: the shrinkB side session (worktree confident-hodgkin-9ba093). Scripts: local `scripts/hetero_tree_shrinkB/` (untracked) =
batch2 `/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/shrinkB/` (shrinkB.py, shrinkB_run.sbatch, shrinkB_readout.sbatch,
shrinkB_test.sbatch, shrinkB_list.txt). Runs: `shrinkB/runs/<key>_<node>_<mode>/{stageB,split,prod,readout}`; logs `shrinkB/logs/`.
Production nodes are read only.

## Design
- Stage B = the node's own stageB/hetero_step.sh refine3D call with maxits=1, run three times in ONE directory with startit=1,2,3
  (refine3D's maxits is run length, so call k runs iteration k with production's which_iter keys, seed, committed sigma2 and fsc
  files). Between calls the shipped references (recvol_stateNN{,_even,_odd}.mrc = the ML-regularised pair the next iteration aligns
  against) go through shrinkB.py into refs_it<k>/ and are passed as vol1..volK; every sidecar of the stem is copied.
- Modes: copy (control, byte-identical copies), dev2 (main arm: emp2-style two-part, inside/outside tree_wiener's focus mask), dev1
  (one part). Shrinkage: each state's deviation from the count-weighted node mean, per Fourier shell, multivariate Wiener
  C (C + A diag(s2_half) A^T)^+ with C = PSD cross-half covariance of the deviations, A = I - 1 w^T; each half shrunk against
  half-map noise (tree_wiener's W2); the node mean is untouched; merged = (even + odd) / 2.
- Then the node's state_split call under the gate env (SIMPLE_SPLIT_CONTRAST=maps SIMPLE_TREE_KEEPMAPS=yes); readout
  f2r_ab_readout.py vs production (and copy), r131 per-map closed/open arm occupancy (make_view.py + arm_readout.py).

## Validation (job 55010199, production final stage-B maps)
copy byte-identical; dev1/dev2 keep the node mean to 2e-7. dev2 kept power of the deviation (inside focus) >20 A .99, 20-14 A .84-.88,
14-10 A .22-.46, <10 A ~0; outside the focus .30-.72 / .06-.13 / ~0. Focus mask 4.8 % (qtoyflex r), 1.7 % (r131 r2), 22 % (qtoynull r).
NOTE: on qtoynull the deviation's cross-half FSC above 20 A is .98 inside / .77 outside: the cross-half covariance keeps a null's
low-resolution split as "signal" (the same blindness as the f2r projected-deviation hazard).

## Jobs
55010293 copy, 55010294 dev2 (fixtures in order: qtoyflex r, qtoy2 r, r131 r2, qtoynull r), submitted 23:33, ~80 cores together.

## Results
- Harness: the copy control reproduces production (qtoyflex r: populations per iteration within 4 particles, 0.999 of labels
  shared, NMI 0.989; same split): the restart effect of three startit-continued calls is nil.
- qtoyflex r: production = copy: stage-B NMI .168, closed F1 .17 (rec .73 pur .10), holo .62 (rec .90 pur .47), split = states.
  dev2: NMI .166, closed .17 (rec .71 pur .09), holo .62; split reallocates (32862/17058/10062): children closed .16, holo .64.
  5.3 % of labels move (NMI .769 vs copy); truth readouts within +-.02 (seed noise of the replays). Same as production.
  The split explains the reallocation: control stage-B states beat the split's EM regrouping by 15.6 nats cross-half (SE 2.6),
  dev2's by -1.4 (SE 1.4) -> EM groups; weakest-pair statistic .277 (z 30) -> .202 (z 20): less distinct states, no truth gain.
- qtoy2 r: production = copy (1.000 labels): NMI .656 ARI .536, closed F1 .55 (rec .99 pur .38), holo .94, split = 5 states.
  dev2: NMI .648 ARI .528, closed .56 (pur .39), holo .94, split = 5 states; 1.5 % of labels move. Same (NMI -.008).
- r131 r2 (seed 2; the closed state separates here in mm2): production = copy (.998 labels): states 10051/11115, closed-row
  enrichment of the best state 1.112, stage-B maps s01 .777/.519, s02 1.032/.519 (corner by the gate rule >= .95 / <= .6),
  split 2 children (test band 16..23 weakest pair .134 z 8.6; stage-B states beat EM groups +4.4 SE .7; effect floor P/N 2.74).
  dev2: states 9972/11194 (13.6 % of labels move), enrichment 1.217, maps s01 .755/.507, s02 1.047/.524 (corner kept in the
  stage-B maps) BUT the split REFUSES: weakest pair .072 (z 6.6, passes), stage-B states vs EM groups +0.3 (SE .3) -> EM groups
  17984/3182, effect floor P/N .657 < 2 -> merged, K* 1. The tree would deliver r2 as one leaf: the closed state is lost (hard gate).
  Mechanism: at r2 the shrink keeps .13 of the deviation power at 14-10 A inside the focus and ~0 outside, which is the split's test
  band (shells 16..23 = 14.6-10.2 A at box 88); stage B stops partitioning on it, the states' maps differ half as much there.
- qtoynull r: production = copy (.998 labels) and dev2 (11 % of labels move): one child each. Null passes.
- REPLICATE mm2_r131_s1 r2 (seed 1): production = copy (.997 labels): enrichment 1.209 / 1.207, s02 .997/.544 (corner), split 2
  (weakest pair .114 z 6.1 PASS; stage-B states +4.0 SE .7; P/N 2.96). dev2: 12.4 % of labels move, enrichment 1.198, s02
  1.003/.556 (corner kept in the maps), split REFUSES at the test itself: weakest pair .053 (z 2.6) < threshold .071 -> K* 1.
  The closed split is lost on both real trees; the seed-2 enrichment gain (1.112 -> 1.217) does not replicate (1.209 -> 1.198).
- Second batch (23:52): chain 55011976 = copy on mm2_r131_s1 r2 (closed corner at r2 there too: s02 .997/.544) then dev2h
  (half dose) on r131 r2, r131_s1 r2, qtoyflex, qtoy2, qtoynull; chain 55011977 = dev2 on r131_s1 r2 then dev1 (one part, strong).
