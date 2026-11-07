# M-estimation / ML formulations for the hetero_tree: status and numbers

Owner: the mest-formulations session (worktree `mest-formulations-e5b42a`). Scripts: `scripts/hetero_tree_mest/` (local, untracked) =
`/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/mest/` on batch2 (`mest_post.py`, `mest_estep.py`, `mest_prior.py`, sbatch files);
runs under `math_reworks/mest/runs/`. Nothing of the f2r session's `f2r_*` files was touched; its dumps and plane caches were read only.

## Mechanism (why soft weights add nothing): stage-B posteriors of the mu gate trees (`runs/post/post.txt`)

node | n | K | mean max | >=.99 | gap p50 (nats) | calibration (accuracy of the argmax per max-posterior bin <.9 / .9-.99 / .99-.999 / >.999)
---|---|---|---|---|---|---
qtoy2/r | 60000 | 5 | .958 | .768 | 11.2 | .435 / .571 / .632 / .717; wrong 32.5 % of non-junk, their mass on a right state median 0.000 (mean .016); purity hard .675 soft .671
qtoy2/r1 | 21937 | 2 | .919 | .449 | 4.1 | .599 / .856 / .938 / .976; wrong 16 %, mass on a right state median .067; purity hard .840 soft .818
qtoy2/r2 | 11472 | 2 | .957 | .698 | 7.6 | .634 / .860 / .936 / .980; purity hard .932 soft .922
qtoy2/r4 | 9461 | 2 | .962 | .761 | 10.2 | .467 / .730 / .863 / .969; purity hard .913 soft .908
qtoy2/r3, r5 (one conformer each, merged leaves) | | 2 | .81-.84 | .09-.15 | 1.7-2.0 | accuracy .97-.99 everywhere
r131/r1 | 25397 | 2 | .826 | .096 | 1.8 | closed refs s1 x1.37 / s2 x0.63; closed rows' mean max .839 vs others .825
r131/r, r2 | | 2 | .82-.86 | .10-.18 | 1.8-2.4 | enrichment x1.03 / x1.22
qtoynull/r, olenull/r | | 2 | .78-.81 | .02-.03 | 1.4-1.6 | the noise-lock gap scale of two identical states

Two regimes. (1) Toy nodes with real conformational differences: gaps 10-30 nats, posteriors one-hot, M-est = hard labels (purity .675 vs .671),
and the mistakes are confident (28 % of the >.999 assignments wrong at the root) with zero posterior mass on the right state: no reweighting by
the stated posterior can recover them, tempering flattens everyone equally. (2) Real nodes (r131) and homogeneous nodes: the evidence gap
distribution sits at the null's (median 1.8 vs 1.6 nats), so the posteriors are soft but uninformative; soft weights then blend states with
noise-level evidence, which is the closed-arm dilution seen on 692k (anc arm r2112 .91 -> .81, r1312 1.07 -> .99; FSC 8.2 -> 7.4 A vs TW 4.5 A).

Algebra: any weighting that depends only on (own leaf, target leaf) makes the M-estimated map a frequency-blind linear combination of the leaf
maps (numerators and CTF^2 denominators add), i.e. a special case of tree_wiener's per-shell linear combination without the per-shell
covariance; dist_0.4 on qtoy2 (FSC 6.8 -> 5.0 A, dev 20-12 .133 -> .121) is exactly that. Per-particle information beyond it exists only where
posteriors vary within a leaf, which is regime (2), where they carry no signal.

## Label-level formulations at frozen stage-B poses (`runs/estep_*`, 15k subsamples, tree_estep_v3 machinery, 4 iterations)

- B0 = tree_estep_v3's D0b1 reproduced exactly (qtoy2 draw NMI .542 -> .489 vs production .658; r131 enrichment 1.34 -> 1.13; nulls 50/50).
- BV (predictive E-step, reference half-map noise inflating the per-state noise + log-det): inflation 0.08-0.12 % of sigma2 on every fixture;
  identical to B0 (qtoy2 .578 vs .579 amax NMI, r131 1.341 vs 1.339). Map uncertainty cannot calibrate the responsibilities.
- BT (Student-t, nu by ML, IRLS weights): nu -> Gaussian (1e4) on qtoy2 and both nulls, weights within +-1 %, junk AUROC of the residual .52
  (junk share 10 %); on r131 nu = 30, weights .58-1.78 (p1-p99) but enrichment unchanged (1.341 vs 1.339 at it 1, 1.151 vs 1.131 at it 4).
  At cryo-EM SNR the whole-image residual cannot see junk; the per-particle amplitude already absorbs scale.

## Final-map formulations (`runs/prior_*`, mest_prior.py = tree_wiener in numpy, validated: emp2 reproduces TW_emp2 to 3 decimals)

Estimators: emp2 (delivered; PSD empirical covariance, two parts), emp1 (no focus mask), clu2 (= tree2, complete-linkage tree),
topo2 (NNLS tree covariance on the hetero_tree's own topology), depth2 (one innovation variance per depth per shell), anc2 (sequential
hierarchical MAP: child = parent + Wiener-shrunk innovation). All are per-shell linear combinations of the leaf maps; a hierarchical
Gaussian tree prior with MAP reconstruction is exactly the tree-covariance Wiener (the sufficient statistics are the leaf maps).

tree | covariance R^2 (40-8 A): clu2 / topo2 / depth2 | half-map FSC median | state readout
---|---|---|---
mm2 qtoy2 (8 leaves) | .293 / .136 / .003 | all 4.73 A (raw 6.79) | dev 20-12 A vs truth: raw .135, emp1 .163, emp2 .153, depth2 .147, anc2 .141, clu2 .126, topo2 .118; 12-8 A all <= 0
mm2 r131 (4 leaves) | .650 / .650 / .564 | all 6.69 A (raw 9.19) | every estimator keeps the r2_s2 corner 1.05/0.56, others identical to +-.01
692k f692_lf_s1 (66 leaves) | - / .157 / .053 | emp2 4.73, topo2 4.42, depth2 4.26, anc2 4.26 A | top-8 closed leaves unchanged to +-.02 (r2112 .91/.68-.70, r1312 1.06-1.07/.65), 0 corners for all

Reading: the leaves' map covariance does not follow the tree's depth or topology (R^2 .003-.16 on the toy and the 692k tree, .56-.65 only
on the 4-leaf r131 where clustering = topology). Structured priors add no state signal (qtoy2 dev <= emp2); on the 66-leaf tree they lift the
nominal FSC by pooling the high-resolution content whose between-leaf covariance is unmeasurable, with the 10 A state readouts untouched.

## Verdict
No formulation beats hard labels + tree_wiener on the fixture readouts; none dilutes the closed state or splits a null (the final-map variants
are label-neutral). Not worth a Fortran port, except optionally a depth-tied design in simple_tree_wiener (~40 lines, leaf names -> depths)
if a deliverable wants the higher nominal leaf resolution on large trees (4.7 -> 4.3 A on 692k) at no readout change.
f2r session (reply 2026-11-06): latent-span regularised references hold r131 (1.377 -> 1.367) but cost qtoy2 (NMI .574-.588 vs .615) and
harden the nulls (.63 -> .84): the same two regimes.

# Port status (branch claude/tree-ports, rebased onto fork/hetero-tree e1d5fdd4c; worktree mest-formulations; builds under
# math_reworks/build_ports*, gates gates/ports_*)

- d142ea8af (was ee1f7749b) piece 1: depth-tied tree_wiener (`depth_wiener=yes` default; the switch restores the clustering-tree rule). Own-tree
  fast gate 15/15 (build_ports1). GATED on the piece 1+2 build (gates/ports_*, readouts ports_readout.txt / ports_s1_readout.txt
  vs mm2): nulls one leaf (identical), qtoy2 identical leaves and counts (NMI .780), r131 closed corner retained in both seeds
  (0.96/0.54 r1_s2; seed 1 1.00/0.53 r1_s1 SEPARATES(map)); final maps deliver depth2 (R^2 .53 vs tree .60).
- 480bc9f2f piece 2 (global chart coordinates, charts/global_coordinates.txt + split_gaps.txt): PARKED on claude/tree-ports-charts
  at the coordinator's request (the f2r session found the composed coordinate reproduces worse than the root chart and seams the
  split valleys). Built and fast-gated; not on the push branch.
- 043b6bcfc (was 27c42dcf6) piece 3: site_split=no default; simple_site_split.f90 (worker: per-site cross-half two-state EM on
  state_mergeup with site-masked references and half-B renumbering, GAPw with a 3-D convex hull), discover_islands in the tree
  (solvent-mean-corrected half variance), the deferred leaf in decide (S_SITESPLIT), the site_split program (commander, router,
  UI), the site split unit suite (hull, GAPw on synthetic two-blob vs rod maps), doc/implementation_notes/hetero_tree_site_split_2026_11.md.
  Built on batch2 (build_ports3 = piece 1 + piece 3, no piece 2), fast gate 15/15 with the new suite. Gates of the proposal running:
  gate 0 PASSED: the six fixtures with site_split off on build_ports3 (gates/ports3_*) give the ports_* partitions, five label files
  byte-identical and qtoy2 identical up to leaf numbering (node completion order under nodepar=4), readouts identical;
  gate 1 = the worker alone on the 17 saved gate nodes (mest/xs_gate1/<key>/summary.txt, array 55008841; every site of each node,
  Python's mskdiam 232 / 240 c2 / 600), compared with the direct readout of discrete_vs_continuous.md sections 7 and 7b
  (mest/xs_gate1_pytable.txt, mest/xs_gate1_compare.py: call, GAPw within .05, label agreement). Expected: 0 DISCRETE except
  igmfix_root site 6 (GAPw .09).
  Gate-1 run 1 (build_ports3): every site collapsed to one group within three EM iterations (the Python kept two everywhere): the
  port scored particles against the masked group map with the mixing prior and soft weights. Fixed in 043b6bcfc: the site path's
  references are the half's consensus plus the site-masked deviation of each group from the groups' weighted mean (xsbg.py's
  residual-vs-deviation score), flat prior, hard labels, full band; GAPw -1 when a lobe lies on solvent. Rebuilt as build_ports4,
  gate 1 rerunning there; gate 0 stays on build_ports3 (its non-site path is unchanged by the fix).
  Run 3 (c4aa7b562: no docking, no FSC shrinkage on the site path; build_ports5): the RNase P leaf nodes read like the Python
  (z1 <= 3, no_call everywhere), the root nodes and olenull keep one strong reproducible factor each (z1 14-38, lobes 5-9 A,
  stopped by G2/B) where the Python's EM never converged; build_ports6 logs the per-particle score gap to compare evidence.
  The first qtoyflex tree with the layer on (run-2 code) ran away: a node whose site_split worker returned kept step S_SITESPLIT
  and was re-decided every poll (122 child additions, 147 leaf lines, crash). Fixed (step returns to S_SPLIT once the summary is
  read; no summary files the leaf; junk particles kept out of the children) and rerunning as gates/ports6_qtoyflex_ss.
  r131 seed 2 (coordinator's question): mm2 vs ports partitions agree .62 / .50 (seeds 2 / 1) after best matching while ports
  vs ports3 are byte-identical; merge_master2's source is ee26d9bbf but built on xfer03p against libopenblaso. gates/mm2b_r131
  (same binary again) and gates/b0_r131 (ee26d9bbf with the ports recipe) settle run-to-run vs binary vs piece 1.

# Second task: the tree as a hierarchical mixture of local PPCA charts (mest_charts.py, runs/charts_<tree>/charts.txt)

Charts = every node's fflex fit (flex_pca_coordinates.txt = delivered MAP z after the gauge rotation; flex_pca_embedding.bin = raw MAP z,
per-particle posterior precision, trailing noise scale + deconvolved z). Frame facts: delivered z = raw z rotated (canonical correlations
.96-.997, the LS map is not orthogonal, |R^T R - I| .7-.9, so the precisions are carried by the LS map); the budget that balances
(observed = flex's deconvolved + measurement) is the MAP z un-shrunk by the map prior with the noise-scaled covariance ns * P^-1
(qtoy2 root axis 1: 177 = 152 + 0.46 * 52). Child -> parent affine maps by precision-weighted LS (R^2 per parent axis .1-.6), composed
to the root; readouts on the arm particles of the toys.

tree | coordinates | PC1 rho(s1) / conformer means / supervised ceiling | closed AUROC (d') | holo AUROC (d')
---|---|---|---|---
qtoy2 (6 discrete arm conformers) | composed | .643 / .600 / .847 | .995 (4.7) | .997 (6.6)
 | root chart alone | .674 / .771 / .757 | .982 (3.0) | .994 (3.7)
qtoyflex (30-conformer continuum) | composed | .571 / .779 / .590 | .819 (1.3) | .919 (2.3)
 | root chart alone | .524 / .760 / .531 | .784 (1.1) | .865 (1.6)

Within nodes (s1 ceiling local / composed / root): qtoy2 r1 .50/.49/.40, r2 .58/.59/.42 (local charts resolve the discrete pairs the root
chart blurs, and the composition keeps it); qtoyflex r1 .22/.21/.26, r2 .25/.25/.16, r3 .23/.23/.35 (on the continuum no chart tracks
s1 within a node beyond .35: the per-particle noise swamps adjacent conformers). The tree's two big arm leaves overlap in s1
(r2 -9.6..-0.1, r1 -8.1..3.9), and the global axis orders their conformer means at .78.

Split criterion (1-D XD along the children's mean-difference direction in the node's chart, tied-variance mixture, K by cross-half
held-out likelihood, gap depth D between the children's modes, p against a uniform-continuum null with the same deconvolved mean/sd
and noise, 100 replicates; the uniform null itself gives D .3-.5 at K* >= 4, the ripple floor):
- qtoy2 root: 8 of 10 conformer-group pairs gaps (D .50-.86, p .01; chart d' 2.8-5.3); the two misses involve the mixed group 22+18.
  r4 (22 vs closed 18, d' 3.5): D .42 (halves .56/.31), p .33 (borderline); r2 (49 vs 43, d' 2.6) and r1 (1 vs 2, d' 1.7): no gap.
  Same-conformer nodes r3, r5 (merged by the tree): D 0.
- qtoyflex: every split D 0: the arm cuts (root 1|2, r2) correctly, but also closed (r1 1|2, d' 2.1) and holo (root 1|3, d' 2.3).
- nulls: D 0, K* 1 (no false gap). r131: D 0 at all three splits (closed split r1: d' 2.2, enrichment 1.25/0.71).

Verdict: coordination gives a per-particle continuous coordinate the tree lacks, cheaply (affine maps on existing outputs), with a modest
gain over the root chart on a continuum (+.06 s1 ceiling, +.035/.05 AUROC) and a real one where nodes resolve discrete sub-states
(closed d' 3.0 -> 4.7, holo 3.7 -> 6.6). The gap criterion has no false gaps and separates discrete from continuum cuts when the chart's
d' >= 3, but is blind at d' <= 2.5, which is where closed/holo sit on qtoyflex and the closed split on r131: it cannot replace the split
test, only annotate strong discrete splits. Not worth a Fortran modelling layer; the only port with value is writing the composed global
coordinates (and the per-split D, p, d') as tree outputs.
