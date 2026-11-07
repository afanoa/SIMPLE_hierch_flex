# hetero_tree ports: status and how to finish (2026-11-07, porting session of "Tree as full PPCA mixture")

Branch `claude/tree-ports` (worktree mest-formulations on the Mac; the coordinator pushes it to fork/hetero-tree as WIP), rebased
onto fork/hetero-tree e1d5fdd4c. Never rebase or amend at or below e1d5fdd4c. No Co-Authored-By trailers. Nothing here pushes.

Cluster root: M=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks. Builds: `sbatch --export=ALL,S=$M/src_portsN,V=$M/build_portsN
$M/build_ports3.sbatch` (cmake Release, BUILD_TESTS=ON, 32 cores, ~5 min, runs the fast gate; result in $V/fast_gate.exit,
$V/fast_gate.log, $V/make.log; binaries in $V/install). Never rebuild a build dir a running job uses: make a new one.
Gate fixtures: `bash $M/launch_gates.sh $M/build_portsN/install ARM qtoy2 r131 qtoynull olenull micnull` and
`SEEDX=1 RUNSFX=_s1 bash $M/launch_gates.sh ... ARM r131` (-> gates/ARM_<ds>; ~20-50 min each, workers via slurm); readouts
`sbatch --export=ALL,ARM=ARM $M/readout.sbatch` (-> gates/ARM_readout.txt) and the seed-1 r131 readout (copy of
bc_s1_readout.sbatch with bc -> ARM). References: gates/mm2_* (production), gates/ports_* (piece 1 + parked piece 2).

## Piece 1: depth-tied tree_wiener (commit d142ea8af)
`depth_wiener=yes` (default): tree_wiener fits one covariance variance per level of the tree's topology, shared by two units down to
their lowest common ancestor, delivers `depth2`; `no` restores the clustering-tree rule. Files: simple_tree_wiener.f90,
simple_hetero_tree.f90 (final_maps), parameters/parse/UI.
- Own-tree fast gate 15/15 (build_ports1 from src_ports1). Gate trees ran on build_ports (piece 1 + parked piece 2):
  gates/ports_{qtoy2,r131,r131_s1,qtoynull,olenull,micnull}, all COMPLETED. Readouts: jobs 55011290 (ports_readout.txt) and
  55011291 (ports_s1_readout.txt); compare with gates/mm2_readout.txt / mm2_r131_s1 readout (leaf counts, label chains, qtoy2
  leaves vs truth, r131 closed corner). Expected: identical leaves/labels (the final maps only change their Wiener filter).
- GATED 2026-11-07 (readouts gates/ports_readout.txt, ports_s1_readout.txt vs mm2_readout.txt): nulls one leaf each (identical);
  qtoy2 identical leaves/counts, NMI .780, same F1s; r131 seed 2 four leaves, closed corner retained (0.96/0.54 at r1_s2; mm2
  1.03/0.52 at r2_s2, root groups 22484/21850 vs 23168/21166 = run-to-run variation); seed 1 four leaves, closed corner 1.00/0.53
  at r1_s1, SEPARATES(map). Final maps deliver depth2 (R^2 .53 vs tree .60). Reported to the coordinator for the push.

## Piece 2: global chart coordinates (commit 480bc9f2f on claude/tree-ports-charts, PARKED)
Built and fast-gated earlier; on the old base (ee1f7749b over ee26d9bbf). Not to be pushed as a port unless the user decides.

## Piece 3: site_split (commits 043b6bcfc + c4aa7b562 + 378254f33; squash before the final push if wanted)
`site_split=no` default. Worker `simple_site_split.f90` (program site_split), the deferred leaf in `decide` (S_SITESPLIT, xsplit/),
discover_islands, state_mergeup's site path (site_focus, site_split_em, site_refs, cross_em(l_site)), tester suite `site split`,
doc/implementation_notes/hetero_tree_site_split_2026_11.md. Proposal and evidence: side-session notes port_proposal_site_split.md,
discrete_vs_continuous.md (sections 7/7b = the Python per-site table, copied to $M/mest/xs_gate1_pytable.txt).

Gate 1 (the worker alone on the saved gate nodes): `cd $M/mest && sbatch --array=1-17%6 --export=ALL,LIST=$M/mest/xs_gate1.list,BUILD=$M/build_portsN/install xs_gate1.sbatch`
(one node per task; results $M/mest/xs_gate1/<key>/summary.txt and <key>/sN/labels.txt; the dir must not exist: move the old one aside),
then `python3 xs_gate1_compare.py xs_gate1.list xs_gate1 /mnt/beegfs/elmlund/afan/igm_cryolo/tree/fixture_igm xs_gate1_pytable.txt`
(python3 = /mnt/nasapps/production/python/3.12/bin/python3). Pass: calls agree with the Python table (0 DISCRETE except igmfix_root
s6, GAPw ~.09), GAPw within .05, group agreement >= .9 at the strong factors (lf_r211 s1).
- Run 1 (build_ports3): every site collapsed to one group in <= 3 iterations. Cause: the port scored particles against the masked group
  map with the mixing prior and soft weights. Fixed: references = the half's consensus + site-masked deviation of each group from the
  groups' weighted mean (xsbg.py's residual-vs-deviation score), flat prior, hard labels, full band.
- Run 2 (build_ports4, $M/mest/xs_gate1_v2): two groups everywhere, but the cross-half factor is far stronger than the Python's
  (z1 10-35 with lobes 2-9 A apart on RNase P nodes and the nulls; Python z1 <= 3.4 there). Cause (diagnosed, fix in progress):
  cross_em docks each group map onto the other half's consensus by its own 3-D shift, so the two group references differ by a
  consensus-gradient dipole and the EM converges to a translation split; the Python never docks. Fix: no docking and no FSC shrinkage
  of the deviations on the site path (commit c4aa7b562; build_ports5 job 55012207, run 3 chained: gate 1 into $M/mest/xs_gate1 with run 2
  moved to xs_gate1_v2, and gates/ports5_qtoyflex_ss).
- Gate 0 PASSED 2026-11-07 (site_split=no on build_ports3, gates/ports3_* vs gates/ports_*): hetero_tree_labels.txt byte-identical
  for r131 (both seeds), qtoynull, olenull, micnull; qtoy2 identical up to leaf numbering (nodes r2/r3 finished in the other order
  under nodepar=4; the contingency table is a permutation matrix), readouts identical (NMI .780, same F1s, closed corner 0.96/0.54
  and 1.00/0.53). The non-site path is unchanged by the later site-path fixes.
- Run 3 (build_ports5, $M/mest/xs_gate1, 57 of 79 sites in): the RNase P leaf nodes now read like the Python (z1 mostly <= 3,
  r1 <= .5, all no_call), but the root nodes keep a strong reproducible factor at one site (r131s2_r s1 z1 16.5 r1 .95, lobes
  9 A; r131s1_r s1 13.9/.94; ctun_root s1 37.9; olenull_r s1 15.7/.87 lobes 5 A) where the Python read z1 <= 5.8; all stopped by
  G2/B (no DISCRETE). The Python's own EM never converged on these nodes (labels re-randomised each step, score gap median
  .1-.9 nats), so the remaining question is the per-particle evidence magnitude: build_ports6 logs the score-gap P10/median/P90
  per E-step (xs_gate1_gap/ on lf_r1312, r131s2_r, olenull_r) for a direct comparison with xsbg.py's "score gap" lines in
  igm_cryolo/tree/fixture_igm/logs/igmfx_gxsb.<job>.out.
- RUNAWAY BUG found by the first qtoyflex tree (ports4_qtoyflex_ss, run-2 code, failed with 147 leaf lines): a node whose
  site_split worker returned kept its step S_SITESPLIT, so decide re-ran every poll (the same children added 122 times). Fixed in the
  tree (the step returns to S_SPLIT once the summary is read; a worker that left no summary files the leaf; junk-state particles
  stay out of the site children): commit 378254f33, in build_ports6 with the gap diagnostic. The ports5_qtoyflex_ss run was cancelled (same bug).
- r131 seed 2 (piece 1 question): mm2 vs ports partitions agree only .62 (seed 2) / .50 (seed 1) after best leaf matching, while
  ports vs ports3 are byte-identical; merge_master2's source is ee26d9bbf (diff -rq empty) but it was built on fsitgl-xfer03p and
  links libopenblaso (OpenMP OpenBLAS) where the ports builds (x6342 nodes, BLA_VENDOR=OpenBLAS) link libopenblas; same AVX2
  instruction mix. Running: gates/mm2b_r131 (same mm2 binary, seed 2 again: run-to-run reproducibility) and gates/b0_r131
  (build_base = ee26d9bbf with the ports recipe, no piece 1). b0 == ports byte-identical clears piece 1; mm2b == mm2 shows the
  difference is binary-bound (build host / BLAS), not run-to-run noise.
- qtoyflex (the coordinator's positive control; production fails there): `launch_gates.sh` has a qtoyflex case (maxdepth 2, qtoy2's
  keys); `ARMKEYS=site_split=yes RUNSFX=_ss bash launch_gates.sh $M/build_portsN/install ARM qtoyflex` and the control without
  ARMKEYS; score with `python3 $M/latent/qflex_leaves.py gates/<run>/1_hetero_tree` against gates/mm2_qtoyflex. ports4_qtoyflex (control, site_split=no, valid for any build) is running; the layer-on tree reruns as ports6_qtoyflex_ss.
- Then gates 2-6 of the proposal (IgM fixture and null twin, nulls, r131, ctun/ctpk/mm3q) only after the coordinator has seen 0-1.

Scripts: this folder's scripts/hetero_tree_mest/ (the Python M-est/charts prototypes, xs_gate1.sbatch, xs_gate1_compare.py,
xs_gate1_pytable.txt); on the cluster also $M/mest/ (same files, plus logs/ and runs/).
