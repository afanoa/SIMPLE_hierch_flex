# Prior: child poses held near the parent's stage-B poses (pose anchor): verdict

Session: pose-anchor side session (worktree loving-wilbur-d4817e, branch claude/loving-wilbur-d4817e), 2026-11-06/07.
Build: merge_master2 (`/mnt/beegfs/elmlund/afan/hetero_tree_20261028/merge_master2/install`). No Fortran written.
Work dir: `/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/anchor/` (A below).

## Verdict: not worth wiring in

The anchor never improves a readout, and its losses grow as it tightens. With a ~2.7-3 deg direction cell (r131, qtoyflex)
it is neutral to slightly worse: the r131 corner is kept, enrichment stays within seed noise, the closed map gains 0.01-0.02
open-arm density, and qtoyflex NMI drops 0.007. Where the data width is below the direction grid and the direction freezes
(qtoy2), it costs the split children 0.01-0.11 in purity and 0.03-0.14 in NMI at both seeds, against a seed spread of
0.002-0.003; r3 and r5 stay leaves. It is 1-2 points more seed-stable but settles 4-18 % away from the default
partition, and that partition is worse.

Why: there is no pose drift to limit, and the poses it holds the child to are noisy. Under prob_beta 0.3 the stored
stage-B poses are single tempered posterior draws: the parent re-draws each iteration (median 2.7 deg on r131, 1.2-1.8 on
qtoy2; 0.25 in the older greedy builds). Holding a child near them locks in that draw noise, which the child's own draws
otherwise average away. This is the same direction as prob_beta 1, which lost the closed corner.

No Fortran prior was built. A Gaussian prior at the data width sits between the default cell and these arms. The
key-only arm already removes two thirds of the drift tail, in-plane flips included (they come with the shift disk and the
direction freedom), and what remains is junk on the toys and non-closed particles on r131.

Suggested row for the verdict table of `doc/implementation_notes/hetero_tree_latent_integration_2026_11_06.md`:
| Pose anchor: child poses held near the parent's stage-B poses (geom cell and shift_sigma from the parent's spread) | r131 corner kept, enrichment within seed noise, closed map +0.01-0.02 open arm; qtoyflex -0.007 NMI; qtoy2 children purity -0.01 to -0.11 (NMI -0.03 to -0.14) where directions freeze, leaves kept; 1-2 points more seed-stable | no |

## Mechanism (production mm2 trees, incoming = parent stage-B poses -> child stage A -> stage B)

- `prob_neigh_mode=geom` already anchors directions. The neighbourhood is the coarse cell (`nspace_sub`, default 500,
  radius ~5 deg on the full sphere) holding the previous direction, and the draw stays inside it. A direction therefore
  never leaves the cell of the node's INCOMING direction. In-plane is a free posterior draw over all rotations: no key
  bounds it in the tables (`prob_inpl` is refine=neigh only). Shifts use a window of +-trs (5 px) per iteration. A
  particle whose last shift posterior is wider than trs/3 searches a disk of 3*shift_sigma (default 30 A).
- Net child drift, incoming -> stage B (`A/drift/*_drift.txt`), on r131 r1/r2: direction median 2.9 deg (one
  5000-grid step), p90 6.0; in-plane median 0.4, p90 1.4; shift median 0.3 A, p90 0.9 A. That equals one parent
  iteration's draw width: nothing accumulates over the six iterations. The toys show the same or less.
- The large moves form a tail (in-plane > 30 deg or shift > 10 A): r131 1.5-1.6 %, qtoyflex r1 0.6 %, qtoy2 0.03-0.3 %.
  On the toys it is junk: 97-100 % of the tail on the qtoy2 nodes, 61 % on qtoyflex r1. On r131 it is not the closed
  particles: closed rows have a 0.07-0.13 % tail vs 1.6-1.7 % for the rest, with in-plane p99 2.7-3.4 deg and shift p99
  1.7-1.9 A (`A/logs/pa_tail_55010190.out`).
- Depth is checked in the only trees that kept depth >= 2 node projects (`ancestral_dsnr/gates/r1211_s1_off`, the real
  RNase P subtree, and `qtoy2_off`; older greedy-pose builds). From depth 1 (16k particles) to depth 3 (3.7k) the drift
  does not grow: direction median 0 and p99 5.7-6.0 deg, in-plane p99 2.0-2.4 deg, shift p99 1.1-1.5 A, tail <= 0.07 %,
  none of it closed (`A/logs/pa_drift2_55010471.out`).

## The key-only test

- Arm `anc`, on stage A and stage B only:
  - `nspace_sub` = 2/(1-cos r), with r the parent's median per-iteration direction change (the 'dist' of the parent's
    stage B over the child's rows), clamped to nspace. This gives r131 r2 1802, r1 1762 and qtoyflex r1 1488. Every
    qtoy2 node clamps, so its direction is frozen.
  - `shift_sigma` = refine3D's own shift-prior estimate from the parent's increments (`shift_prior_sigma2`). The measured
    spread is below its 0.25*lp floor, so it sits at the floor: 3.75 A at stage A, 2.5 A at stage B.
- Arm `def`: the saved calls, unchanged.
- Per node: production (def at its seed), def at the other seed (the replicate) and anc at both seeds. Production seeds
  are 2 for r131 and 1 for qtoy2 and qtoyflex.
- Harness check: def at production's seed reproduces production exactly (r131 r2, 1.000 of 21166 labels).
- On r131 r2 the anchor cuts the drift tail from 1.6 % to 0.5-0.6 % (in-plane p99 62 -> 4.8 deg; closed rows' tail -> 0),
  halves the last-iteration jitter (dist median 2.7 -> 1.35 deg) and relabels 12 % of particles, against 8 % for a seed
  change. Closed rows relabel at 10.6-11.3 % under the anchor vs 6.7 % for a seed change (`A/logs/pa_change_55010884.out`).

r131 r2, the closed corner (production seed 2):

| arm | closed enrichment | closed-state map closed/open | split | labels = production |
| --- | --- | --- | --- | --- |
| production | 1.112 | 1.032 / 0.519 | 2 | 1.000 |
| default s1 | 1.099 | 1.032 / 0.524 | 2 | 0.920 |
| anchor s2 | 1.111 | 1.029 / 0.542 | 2 | 0.882 |
| anchor s1 | 1.088 | 1.028 / 0.544 | 2 | 0.878 |

r131 r1 (the closed state separates at r2 in this tree, not r1; production seed 2):

| arm | closed enrichment | best closed-arm map closed/open | split | labels = production |
| --- | --- | --- | --- | --- |
| production | 1.247 | 0.873 / 0.625 | 2 | 1.000 |
| default s1 | 1.245 | 0.873 / 0.620 | 2 | 0.926 |
| anchor s2 | 1.243 | 0.869 / 0.631 | 2 | 0.877 |
| anchor s1 | 1.237 | 0.867 / 0.633 | 2 | 0.873 |

qtoyflex r1, the closed cut (production seed 1). The arm continuum is split near chance in every arm (s1 AUROC ~0.57):

| arm | stage-B NMI | children NMI | children closed F1 (recall/purity) | arm s1 AUROC between children | labels = production |
| --- | --- | --- | --- | --- | --- |
| production | 0.158 | 0.179 | 0.53 (0.83/0.39) | 0.569 | 1.000 |
| default s2 | 0.157 | 0.177 | 0.53 (0.84/0.39) | 0.572 | 0.957 |
| anchor s1 | 0.149 | 0.172 | 0.52 (0.83/0.38) | 0.563 | 0.920 |
| anchor s2 | 0.150 | 0.172 | 0.52 (0.83/0.38) | 0.565 | 0.918 |

qtoy2 depth-1 nodes (production seed 1). Each cell is the split children's mean purity over genuine particles / their NMI
with the structure combo over genuine particles; the children equal the stage-B states (`A/final_readout.txt`):

| node | production | default s2 | anchor s1 | anchor s2 | labels anchor = production |
| --- | --- | --- | --- | --- | --- |
| r1 (conformer 1 vs 2) | 0.744 / 0.222 | 0.745 / 0.223 | 0.639 / 0.084 | 0.639 / 0.084 | 0.816 |
| r2 (49 vs 43) | 0.928 / 0.678 | 0.928 / 0.675 | 0.920 / 0.644 | 0.920 / 0.643 | 0.959 / 0.958 |
| r4 (22 vs closed 18) | 0.914 / 0.694 | 0.914 / 0.696 | 0.905 / 0.652 | 0.905 / 0.653 | 0.964 |
| r3 (combo 28) | leaf (1 child) | leaf | leaf | leaf | 0.797 / 0.795 |
| r5 (holo) | leaf (1 child) | leaf | leaf | leaf | 0.801 |

r3 and r5 stay leaves in every arm and at every seed. The r1 loss is not drift removal: the anchor's net displacement is
nearly the default's (incoming -> stage B, direction p90 3.0 vs 3.5 deg, in-plane p90 0.74 vs 0.73 deg). What the
frozen direction removes is the child's own re-drawing of its poses at each iteration.

Seed reproducibility, as stage-B labels matched over permutations (`A/logs/pa_agree_55012024.out`):

| node | anchor s1 vs s2 | default s1 vs s2 |
| --- | --- | --- |
| r131 r2 | 0.940 | 0.920 |
| r131 r1 | 0.938 | 0.926 |
| qtoyflex r1 | 0.963 | 0.957 |

## What ran where (all under A)

- Production drift: `pa_drift.sbatch` (job 55009949) writes the node dumps `drift/<tree>_<node>_{n,a,b}.txt`, the parent
  dump `drift/<tree>_r_b.txt` and the tables `drift/*_drift.txt`. `pa_tail.sbatch` (55010190) splits the tail by truth
  class. `pa_drift2.sbatch` (55010471) measures drift by depth.
- Replays: `pa_replay.sbatch` array 55010141 over `pa_list1.txt` (25 replays, `%2`). Each replay is in
  `replays/<tree>_<node>_<arm>_s<seed>/` with `stageA/ recB/ stageB/ split/` (each step's `cmd.txt` and `step.log`) and
  `readout/{labels,drift,arm}.txt`. Logs are `logs/pa_replay_55010141_<task>.out`.
- Readouts:
  - `pa_post.sbatch`: the full drift and tail per replay.
  - `pa_agree.py`: labels, truth readouts and pairwise agreement.
  - `pa_change.py`: which particles an arm relabels.
  - `pa_qflex.py`: the qtoyflex arm cut.
  - `pa_qtoy.py`: qtoy2 children purity.
  - `pa_collect.sh`: one block per node.
  - `pa_final.sbatch`: all of the above into `A/final_readout.txt`.

## Open jobs and how to resume

- No open jobs. All 25 replays finished (array 55010141). Every table above is in `A/final_readout.txt`, written by
  job 55013462 (`pa_final.sbatch`); rerun that script to regenerate it.
- More seeds or nodes: append lines "KEY NODE ARM SEED TRUTH XA XB" to a list (XA/XB are the stage A/B extra keys joined
  by '+', or '-' for none), then run `LIST=<list> sbatch --array=1-N%2 --export=ALL,LIST=<list> pa_replay.sbatch`. A
  replay directory that already exists is skipped; move it aside to rerun. The script exports
  `SIMPLE_SPLIT_CONTRAST=maps` (the gates' split env).
