# Pose anchor (child poses held near the parent's stage-B poses): status

VERDICT (00:15, 2026-11-07): NOT worth wiring in. Handoff: /mnt/beegfs/elmlund/afan/hetero_tree_20261028/handoff/prior-parentpose/
VERDICT.md (+ scripts); all tables in math_reworks/anchor/final_readout.txt; no open jobs. Repo copy untracked in worktree
loving-wilbur-d4817e (scripts/hetero_tree_handoff/prior-parentpose/); the commit was blocked by the session's permission
classifier. The anchor never improves a readout. r131 and qtoyflex (~2.7-3 deg cell) are neutral to slightly worse. qtoy2
(frozen direction) costs the children 0.01-0.11 purity and 0.03-0.14 NMI; r3 and r5 stay leaves. Why: the stored poses are
single prob_beta 0.3 draws, and holding a child to them removes its own re-drawing.

Owner: the pose-anchor side session (worktree loving-wilbur-d4817e). Cluster dir `math_reworks/anchor/` (scripts `pa_drift.py`,
`pa_tail.py`, `pa_replay.sbatch`, `pa_list1.txt`, `pa_collect.sh`; replays under `replays/<key>_<arm>_s<seed>/`, drift dumps
under `drift/`, `drift2/`). No Fortran written.

## Mechanism facts (production mm2 trees, incoming = parent stage-B poses -> child stage A -> stage B)

- prob_neigh_mode=geom confines a particle's direction to the coarse cell (nspace_sub=500, radius ~5 deg, full sphere) that
  holds its INCOMING direction: the neighbourhood is the cell of prev_proj, the draw stays in it, so the next iteration's
  cell is the same. Directions can never random-walk inside a node. In-plane is a free posterior draw over all 256
  rotations (no key restricts it in the tables; prob_inpl is refine=neigh only). Shifts: window +-trs (5 px) per
  iteration; a particle whose last shift posterior is wider than trs/3 searches a disk of 3*shift_sigma (default 30 A).
- mm2's per-iteration direction change (parent stage-B 'dist', median 2.7 deg r131, 1.2-1.8 qtoy2, 3.0 qtoyflex) is the
  prob_beta=0.3 posterior draw width; the older greedy builds have 0.25 deg.
- Net child drift in -> stage B (r131 r1/r2): direction median 2.9 deg (one 5000-grid step), p90 6.0; in-plane median 0.4,
  p90 1.4; shift median 0.3 A, p90 0.9 A. Equal to ONE parent iteration's jitter: no accumulation over the 6 iterations.
- The only large moves are a tail (in-plane > 30 deg or shift > 10 A): r131 1.5-1.6 %, qtoyflex r1 0.6 %, qtoy2 0.03-0.3 %.
  It is junk on the toys (97-100 % of the tail on qtoy2 nodes; 61 % on qtoyflex r1, junk tail rate 4.5 %) and NOT the
  closed particles on r131 (closed reference rows tail rate 0.07-0.13 % vs 1.6-1.7 % for the rest; closed in-plane p99
  2.7-3.4 deg, shift p99 1.7-1.9 A). Same at depth 2-3 in the older qtoy2_off tree (direction median 0, in-plane p99 ~1 deg).
- Depth (older greedy-pose builds, ancestral_dsnr/gates, the only trees with kept depth >= 2 node projects): the real RNase P
  subtree r1211_s1_off from depth 1 (16k) to depth 3 (3.7k): direction median 0 / p99 5.7-6.0 deg, in-plane p99 1.95 -> 2.0-2.4,
  shift p99 1.1 -> 1.2-1.5 A, tail <= 0.07 %, none on closed rows. No drift growth with depth or with weaker references.

## Key-only anchor replays (merge_master2 build; production seeds; arm keys on stage A and stage B only)

Width from the data: geom cell radius = the parent's median per-iteration direction change over the child's rows
(nspace_sub = 2/(1-cos r): r131 r2 1802, r1 1762, qtoyflex r1 1488, qtoy2 nodes clamp to nspace = frozen direction);
shift_sigma = refine3D's own shift-prior estimate from the parent's increments (shift_prior_sigma2), which sits at its
0.25*lp floor because the measured spread is smaller (3.75 A stage A, 2.5 A stage B; default 10 A).
Harness: default keys at production's seed reproduce production exactly (r131 r2: 1.000 of 21166 labels).

r131 r2 (closed corner; production seed 2):
  arm          closed enrich | closed-state map closed/open | split | labels = production
  production   1.112         | 1.032 / 0.519                | 2     | 1.000
  default s1   1.099         | 1.032 / 0.524                | 2     | 0.920
  anchor s2    1.111         | 1.029 / 0.542                | 2     | 0.882
  anchor s1    1.088         | 1.028 / 0.544                | 2     | 0.878
  The anchor removes the drift tail (in->B geodesic p99 62 -> 9 deg) and halves the last-iteration jitter (dist median
  2.7 -> 1.35 deg), relabels more than a seed change does (12 % vs 8 %), and leaves the corner where it was (open arm
  +0.02 in the closed map at both seeds). It removes 2/3 of the tail (in-plane > 30 deg 1.3 -> 0.4-0.5 %, shift > 10 A
  1.4 -> 0.2-0.3 %, closed rows' tail -> 0): the in-plane flips go with the shift disk and direction freedom, so a Fortran
  in-plane prior would only act on the remaining 0.5 %. Relabelled: seed change 8.0 % (closed rows 6.7 %), anchor
  11.8-12.2 % (closed rows 10.6-11.3 %).

qtoyflex r1 (closed cut; production seed 1):
  arm          stage-B NMI | children NMI | children closed F1 (rec/pur) | arm s1 AUROC children | labels = prod
  production   0.158       | 0.179        | 0.53 (0.83/0.39)              | 0.569                 | 1.000
  default s2   0.157       | 0.177        | 0.53 (0.84/0.39)              | 0.572                 | 0.957
  anchor s1    0.149       | 0.172        | 0.52 (0.83/0.38)              | 0.563                 | 0.920
  anchor s2    0.150       | 0.172        | 0.52 (0.83/0.38)              | 0.565                 | 0.918
  Slightly worse at both seeds, beyond the seed spread (+-0.002 NMI).

r131 r1 (closed separates at r2 in this tree, not r1; production seed 2):
  arm          closed enrich | best closed-arm map closed/open | labels = prod
  production   1.247         | 0.873 / 0.625                   | 1.000
  default s1   1.245         | 0.873 / 0.620                   | 0.926
  anchor s2    1.243         | 0.869 / 0.631                   | 0.877
  anchor s1    1.237         | 0.867 / 0.633                   | 0.873
  Same pattern as r2: the closed map loses a little closed density and gains open density under the anchor at both seeds.

Seed reproducibility (stage-B labels, permutation-matched; logs/pa_agree_55012024.out): anchor s1 vs s2 / default s1 vs s2 /
anchor vs production: r131 r2 0.940 / 0.920 / 0.88; r131 r1 0.938 / 0.926 / 0.875; qtoyflex r1 0.963 / 0.957 / 0.92.
The anchor is 1-2 points more seed-stable but settles 8-12 % away from the default basin, and that basin is not better.

qtoy2 (production seed 1; anchor = frozen direction, the data width is below the grid). Each cell is the children's mean
purity / NMI over genuine particles:
  node  production     default s2     anchor s1      anchor s2
  r1    0.744 / 0.222  0.745 / 0.223  0.639 / 0.084  0.639 / 0.084
  r2    0.928 / 0.678  0.928 / 0.675  0.920 / 0.644  0.920 / 0.643
  r4    0.914 / 0.694  0.914 / 0.696  0.905 / 0.652  0.905 / 0.653
  r3, r5: one child (leaf) in every arm and seed.
  On r1 the net displacement equals the default's (direction p90 3.0 vs 3.5 deg): the loss is the child's own re-drawing.
