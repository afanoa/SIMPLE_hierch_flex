# f2r-global handoff: flex -> refinement couplings and the global coordinate (2026-11-06)

Session: the "flex -> refinement" side session of the tree-as-PPCA-mixture work. Everything below ran on batch2; scripts are
copied beside this file (`f2r_*.py`, `f2r_*.sbatch`, `f2r_*_list*.txt`) and in `scripts/hetero_tree_proto_latent/` of the Mac
checkout (untracked). The detailed record with every table is `f2r_flex_to_refinement_status.md` (copied here too).
Runs live under `/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/latent/` (`L` below): `runs/<key>/f2r_*`,
`ab_sm/`, `ab_seed/`, `global/`, `shared/`, `dumps/`, `cache/`, `logs/`.

## Verdicts (final)

1. Flex state weights as a prior on the refinement's state posterior (frozen poses, Python stand-in): real information
   (beats a stickiness control) but the held-out calibration switches it off as references improve (qtoy2 it-1 NMI gain
   +0.040 at 15k -> +0.018 at 60k; r131 full node alpha 0 from iteration 3); frozen, it overrides held-out data
   (-0.004..-0.025 nats/particle) and anchors null splits. Alpha = 0 on eight null runs. NOT worth building.
2. Two-distribution (learned per-state Gaussians over latent + objective indicators): worse than the calibrated product
   (qtoy2 -0.044 vs +0.077 nats), overconfident and negative on nulls. Closed.
3. Held-out latent gain as a discrete-vs-continuum node test (inode fixtures): f50 0, cont +0.0006, null +0.0005, f15 +0.028:
   responds to a noisy minority reference, not to continuity. Closed.
4. Latent-span regularised state references (two FSC weightings): stops the drift on r131 but costs signal on qtoy2 and
   makes null deviations look reproducible; the raw-FSC weighting restores calibration and returns the drift. Closed.
5. prob_state_marg (RELION-style state evidence summed over poses) in production (merge_sm build = claude/state-eval tip +
   switch; fast gate 15/15): r131 enrichment 1.450 -> 1.466, 95 % labels unchanged, nulls one child. No gain. The seven-file
   switch and the three-file `write_states` program are parked on branch `claude/prob-state-marg` by the parent session.
6. Posterior temperature by held-out prediction over a 65-pose neighbourhood x states: pose-only beta* 0.10 qtoy2, 0.15-0.20
   r131, 0.05-0.10 nulls; joint 0.2-0.3. The parent's gate trees at 0.15 lost the closed corner and merged a qtoy2 pair, so
   0.3 stands; no per-iteration estimator.
7. Latent seeding of stage A (flex + site-factor mixture labels replacing site_em's seeds; node replays with the tree's own
   calls, `f2r_replay*.sbatch`): better seeds (qtoy2 NMI 0.571 vs 0.309) but the finished stage B ties on both toys
   (qtoy2 0.650 vs 0.656; qtoyflex 0.172 vs 0.168) and loses the closed-row enrichment on r131 r1 (1.02 vs 1.25 at both
   seeds); nulls one child. Harness facts: the replay reproduces production bit-for-bit at production's seed; the split step
   needs `SIMPLE_SPLIT_CONTRAST=maps` exported (gates/<arm>/split_env.txt); fixture/child projects carry the parent's stacks
   and `prg=selection` refuses them (`write_states` writes states tree-style). No port.
8. Composed global chart coordinate (child -> parent affine maps on fflex charts, `f2r_global.py`, `L/global/`):
   closed AUROC global / root chart r131 .67-.69 / .62, qtoy2 .994 / .982, 692k m33 .75 / .77, base92 .77 / .77, base80
   .82 / .77; d' <= 1.3 on real data; closed-enriched leaves no more isolated than other sibling pairs (r131 d' 1.35-1.51);
   the 692k closed pair r1111221_s1/_s2 at d' 1.12 (overlap .74) is in the lower quartile of 112 sibling pairs, so a latent
   merge rule would fuse it. Reproducibility: root chart .998-1.000 across seeds and recipes, composed coordinate .89 across
   r131 seeds (leaf NMI .19) and .45-.96 across 692k recipes (NMI .24-.45). DECISIVE: valleys of .8-1.0 on split axes (r131,
   the qtoyflex continuum, all three 692k roots) vanish in the root chart (one affine frame): they are SEAMS of the
   piecewise-affine composition (children reach the root through different maps, R^2 .1-.6); cross-node d' of 20-54 on 692k
   are seams too. Only qtoy2's discrete conformers show a real root-chart valley (.6). Nothing to port.

## In flight: the seam-free coordinate (`f2r_shared.py`, `f2r_shared.sbatch`)

One tree-wide basis (all charted nodes' flex modes, band-limited to 15 A, orthonormalised, top q = 30 by SVD; `--basis root` =
the root's modes as the control), every particle embedded at its final leaf stage-B pose (`f2r_leafposes.py` assembles
`dumps/root_<key>/final_ptcl3D.txt` from the leaf nodes' stage-B dumps) or at the input poses, with its own b, G; z_ML with
covariance (G + ridge)^-1 and z_MAP under one Gaussian. Tests: closed AUROC / d', per-leaf isolation, the named closed pair,
truth equivalence, valleys on the root-split and closed-pair axes, Bhattacharyya merge criterion; cross-tree CCA with
`f2r_shared.py compare`.
Jobs: gates array 55008359 (`f2r_shared_list.txt`: mm2_r131, mm2_r131_s1, pb15_r131, mm2_qtoy2, pb15_qtoy2, mm2_qtoyflex);
692k 55008360 (full692_m33_s1, makes `cache/f2r_shared_full692.npy`) then 55008361 (base92, base80), `f2r_shared_list692.txt`.
Outputs: `L/shared/<key>_<all|root>_<final|input>/report.txt`, `summary.json`, `shared.npz`; logs `L/logs/f2rshared_<job>_<task>.out`.
Resume: `cd L && sbatch --array=1-6%3 --export=ALL,LIST=L/f2r_shared_list.txt f2r_shared.sbatch` (gates) and
`sbatch --array=1-3%1 -c 16 --mem=160G -t 08:00:00 --export=ALL,LIST=L/f2r_shared_list692.txt f2r_shared.sbatch` (692k;
dumps, poses and caches are reused when present). Cross-tree comparisons: `python3 f2r_shared.py compare L/shared/A L/shared/B A B`
with `MEST_DIR=.../math_reworks/mest MEST_LATENT=L PPINFO=.../math_reworks/ppinfo` exported (python3 = /mnt/nasapps/production/python/3.12/bin/python3).
Decision rule agreed with the parent: the seam-free coordinate is worth porting only if it isolates the closed state better
than the root chart (AUROC / d' and the closed pair's rank among sibling pairs), reproduces across seeds and recipes at least as
well as the root chart, and supports merges that keep the closed pair distinct.

## Partial numbers (gate trees done 2026-11-06 23:55; 692k running)

Seam-free coordinate, tree-wide basis q = 30 (root-modes control q = 10 in the same pipeline), z_MAP unless stated:
- Closed AUROC (d'): r131 s2 .694 (.70) final poses / .699 (.71) input poses / root-modes control .577 (.28); r131 s1 .705 /
  .714 / .570; pb15 r131 .705 / .720 / .573. For reference flex's own root chart .616, the composed (seamed) coordinate .67-.69.
  qtoy2 .988 (3.2) / .984 / .970 (2.7) (composed 4.5); qtoyflex .843 (1.40) / - / .753 (.97) (composed 1.30).
- Isolation: r131 leaf pairs d' .68-1.09, Bhattacharyya overlap .83-.92 for every pair, closed-enriched leaves no more isolated
  than the others (the composed frame's cross-node d' 5.6-6.7 were seams); qtoy2's closed leaf is the most isolated sibling
  (d' 1.57; holo leaf 2.05) but BC > .5 for 15 of 28 distinct pairs; qtoyflex pairs d' .93-1.20, BC .80-.83.
- Valleys on the root split axis (per half): r131 (all three trees) none, K* 1; qtoyflex continuum none; qtoy2 .35 @ -.25 sd
  (root chart .6): the seams are gone and so are the valleys on real data and on the continuum.
- Truth equivalence (pb15_qtoy2): equivalent pair d' .65 vs distinct >= 1.01, AUROC 1.0 (root-modes .32 vs .73).
- Reproducibility (canonical correlations 1-3 of z_MAP; variate-1 extrema overlap, chance .05 | leaf NMI):
    final poses: r131 s2 vs s1 .89 .89 .87 | .69 (NMI .19); s2 vs pb15 .88 .87 .86 | .68 (NMI .34); qtoy2 .97 .96 .95 | .74-.82
    input poses: r131 s2 vs s1 .996 .995 .994 | .94; s2 vs pb15 1.000 | .98; qtoy2 .999 | .97
    root-modes control at final poses: .86 / .84 / .97.
  => the loss of reproducibility is the final poses (each tree's leaf refinements imprint its partition), not the basis;
  at the input poses the coordinate is reproducible and as informative, but on real data it is weak (d' <= .8).
- Merge criterion (Bhattacharyya > .5): would merge every r131 pair and half of qtoy2's distinct pairs: unusable as is.
Pending: 692k m33 / base92 / base80 (closed pair r1111221_s1/_s2 rank, valleys, CCA across recipes); the RECOVAR map route is
not warranted unless the coordinate isolates the closed state (it does not on r131).
