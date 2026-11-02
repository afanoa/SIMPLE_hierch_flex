# hetero_tree: mask core, site factors and a variance stop (audit items 6.1 / 6.2 / stop, November 2026)

Branch `claude/tree-mask-discovery-stop` off `claude/tree-dist` (9005a8c77). Three gaps of
`hetero_tree_algorithmic_audit_2026_11_01.md` worked through in the order asked: (1) the node mask's core by a relative
share rule with a data-derived falloff, (2) site discovery by a weighted SVD with a share threshold, (3) a stopping rule
from the variance map. Each item was prototyped in Python on the saved fixture maps first (the base arm of the math-reworks
gate replays, `math_reworks/gates/base_*` on batch2, maps kept; the r131 seed-1 fixture `route9/guardfix/runs/r131_md1_s1`
and the arms root `subcls_mask/arms/tree/a` from the 2D sub-class session), then, where the prototype left a case, ported to
Fortran and replayed on the gate fixtures with the same seeds and keys as the base arm (`mask_disc/gates/m23_*`).

Prototype: `scripts/cls_expansion_proto/tree_mask_discovery_stop/vmask_proto.py` (port of `half_variance`, `falloff_mask`
and `discover_sites`, validated against the trees' own `mask.mrc` at Dice 0.98-1.00 and against their `sites.txt`), run per
node by `proto.sbatch`; work dir `/mnt/beegfs/elmlund/afan/hetero_tree_20261028/mask_disc/` (proto/<node>/report.txt,
builds/, gates/). The three map sets of a node are the ones the tree uses: the node-mask maps (flex states at a project
root, the ten cheap children at a subtree root, the parent's recB0 state maps below the root under `oneflex`), the recB0
state maps (discovery) and the stage-B state maps (the split).

One statistic recurs below: the share of a node's positive cross-half variance held by the top x % of the molecule's
support voxels, and its inverse, the molecule fraction covered by the level set that holds a share S of that variance.

## 1. Mask core by a relative rule (audit 6.1): no-go at the prototype

The proposal: replace `core = v > 0.2 max(v)` by the smallest voxel set inside the molecule support holding a fixed share
S of the cross-half variance (S from the fixtures in 0.5-0.8), and take the Gaussian falloff length from the variance map's
own correlation length instead of 6 px; keep the 3 px dilation and the falloff shape.

What the maps show: the cross-half variance is diffuse on every node, real heterogeneity included. The top 5 % of the
support voxels hold 0.11-0.14 of the positive variance on the nulls and 0.16-0.35 on the real nodes; a level set holding
half of it already covers 15-25 % of the molecule on the real nodes (14-16 pieces at the r131 root), and the 3 px dilation
plus the falloff then fill the rest. The 0.2-of-peak core is smaller than every share core at S >= 0.5 on every node, so
the share rule enlarges every mask instead of localising it.

| node (mask maps) | top-5 % share | old core / mask, % of molecule | S 0.6 core / mask | S 0.7 core / mask | arm coverage old -> S 0.7 | acf half-lag (px) |
|---|---|---|---|---|---|---|
| r131 root, seed 2 (10 cheap children) | 0.25 | 13.2 / 94.0 | 22.0 / 98.3 | 30.1 / 99.9 | closed arm 0.97 -> 1.00 | 5.6 |
| r131 root, seed 1 (10 cheap children) | 0.25 | 12.8 / 93.5 | 21.9 / 97.8 | 30.1 / 99.9 | closed arm 0.98 -> 1.00 | 5.6 |
| arms root (10 cheap children) | 0.16 | 12.2 / 99.8 | 31.8 / 100 | 40.7 / 100 | four arms 1.00 -> 1.00 | 4.6 |
| qtoy2 root (8 flex states) | 0.28 | 11.7 / 51.7 | 17.7 / 56.0 | 25.1 / 65.5 | - | 8.0 |
| qtoy2 r1 (parent recB0; its split is refused) | 0.19 | 16.2 / 90.0 | 29.3 / 96.4 | 39.1 / 98.9 | - | 6.6 |
| qtoy2 r2 / r3 / r4 (parent recB0) | 0.35 / 0.34 / 0.35 | 6.0 / 29.4, 7.3 / 25.6, 8.4 / 45.2 | 15 / 74, 16 / 67, 14 / 64 | 23 / 93, 24 / 87, 21 / 86 | - | 6.8 / 8.9 / 11.0 |
| r131 r1 / r2 (parent recB0) | 0.22 / 0.23 | 8.6 / 61.2, 8.1 / 60.5 | 27 / 95, 27 / 94 | 37 / 98, 36 / 98 | closed arm 0.98 / 0.95 -> 1.00 | 7.3 / 7.7 |
| ole root (8 flex states, box 360) | 0.27 | 14.9 / 45.1 | 19.0 / 53.5 | 26.2 / 66.3 | - | 18.0 |
| qtoynull / olenull / micnull root (8 flex states) | 0.14 / 0.14 / 0.11 | 11 / 99, 32 / 100, 63 / 100 | 34 / 100, 38 / 100, 42 / 100 | 44 / 100, 49 / 100, 52 / 100 | - | 3.5 / 5.5 / 4.3 |

The 93 % mask at the r131 root is therefore not the 0.2-of-peak rule's doing: that core is 13 % of the molecule, in
fourteen pieces spread over it (the variance peaks inside the closed arm, as the 2D sub-class session found, but 75 % of
it lies elsewhere), and the dilation and falloff (mask above 0.5 out to about 10 px from every core voxel) join the
pieces. A relative rule cannot localise a variance that is not localised. The falloff length from the autocorrelation
half-lag of the variance map would be 3.5-5.6 px on the nulls and the cheap-children maps and 6.6-11 px (18 px at OLE's
0.95 A/px) on the real nodes' parent-recB0 maps, so it would widen the masks of real depth-1 nodes further; it is not
adopted either. Verdict: **no-go**, no Fortran default change, no gate run (the gate inputs would be a 95-100 % mask
against today's 93 % at r131 and larger masks everywhere else). `MFRAC`, `MDIL`, `MFALL` stay.

What does localise: the recB0 state maps after the focused stA0 pass. On the arms root their variance puts 0.48 in the top
5 % of voxels and the 0.2-of-peak mask covers 22.8 % of the molecule with all four arms inside it (the cheap-children mask
is 99.8 %); on the r131 root 0.23 and 58.6 % with the closed arm covered 0.96. That is the map the children already use
(`oneflex_inputs`), and it says the root's cheap-children mask is the weak one, not the rule.

## 2. Site discovery by a weighted SVD with a share threshold (audit 6.2)

The change (`discover_sites`): for a core component with at least 2 x DMINVOX voxels, the voxel-by-state deviation
matrices of the two halves (columns weighted by the square root of the state weight) replace the 2-means on normalised
profiles. The even half's K x K Gram matrix gives the state directions and the singular values (`jacobi`, as `top_dirs`
did); each factor's spatial map is the even and the odd deviations projected on its direction; a factor is kept when the
two maps correlate at a Fisher z of at least `DFACT_Z` = 3 over `n_eff` = half the component's voxel count times the
15 A band's share of Fourier space (the `mask_neff` form), and when its singular value holds at least `DFACT_SHARE` = 0.2
of the component's variance. With two or more factors kept, every voxel goes to the kept factor whose half-averaged map is
largest there; the pieces of at least max(20, DMINVOX / 4) voxels become the factor's core, as before. One kept factor
leaves the component one candidate (`fac = 0`), as the 2-means did. `DFACT_COS`, `DFACT_ITS` and `top_dirs` go;
`DFACT_SHARE` changes meaning (from a voxel group's share of the smoothed variance to a singular value's share).

Prototype on the recB0 maps (old core; the share cores of section 1 as the discovery core were also tried and lose the
closed arm, coverage 0.71-0.73, so the discovery core stays at 0.2 of the peak):

| node | 2-means sites (tree) | SVD, share 0.10 / 0.15 / 0.20 / 0.25 | what changes at 0.20 |
|---|---|---|---|
| r131 root, seed 2 | 5 | 5 / 4 / 4 / 3 | component 10 (682 voxels, n_eff 1) no longer splits into 10f1 / 10f2; the strongest site stays 5f1 (seed state 6, share 0.49 -> 0.51, mask 14 % of the molecule), site 2 keeps the closed arm at 0.83; at 0.25 5f1 / 5f2 (0.15) fold too |
| r131 r1 | 7 | 6 / 5 / 5 / 5 | best closed-arm site coverage 0.88 -> 0.96 |
| r131 r2 | 3 | 3 / 3 / 3 / 3 | identical, closed arm 0.86 |
| arms root | 5 | 5 / 4 / 4 / 4 | 1f2 (2-means share 0.13; SVD 0.04) folds; arms 1, 2, 4 recalled >= 0.94 as before, arm 3 missed as before |
| qtoy2 root, r2, r3 | 6, 4, 2 | identical at every share | - |
| qtoy2 r1 (refused node) | 14 | 17 / 16 / 14 / 13 | - |
| qtoy2 r4 | 6 | 4 / 4 / 4 / 3 | - |
| ole root | 5 | 5 / 4 / 4 / 4 | - |
| qtoynull / olenull (c1 port) / micnull | 20 / 7 / 0 | 23, 8, 0 / 20, 7, 0 / 20, 7, 0 / 20, 7, 0 | share 0.10 adds sites on the nulls |

Share 0.10 over-counts (the audit's reproducibility-only finding again), 0.25 starts folding a real seed-2 factor at the
r131 root; 0.15 and 0.20 agree everywhere but qtoy2 r1 (16 vs 14). `DFACT_SHARE` = 0.2. Every component's factor
shares, cross-half r and n_eff are now one log line per component.

Gate replays: section 4 (m23 under the old ranking, m23d under variant D, a qtoy2 seed-2 pair, the corner nodes' split
tests over five draws); verdict and commit in section 5.

## 3. A stopping rule from the variance map

The proposal: a node whose variance map after refinement shows no localised structure, by the share rule of section 1
(the level set holding the fixed share covers more than a set fraction of the molecule), becomes a leaf regardless of
population and depth.

Which map: the stage-B states' cross-half variance separates the nulls from every real split; the earlier maps do not.
Molecule fraction covered by the level set holding 0.7 of the positive variance (top-5 % share in brackets):

| stage | nulls (qtoynull / olenull / micnull) | refused splits (qtoy2 r1, r131 r2) | real splits |
|---|---|---|---|
| node-mask maps | 44 / 49 / 52 (0.14 / 0.14 / 0.11) | 39, 36 (0.19, 0.23) | r131 root 30, arms 41, qtoy2 root 25, qtoy2 r2-r4 21-24, r131 r1 37, ole root 26 (0.16-0.35) |
| recB0 states | 47 / 48 / 50 (0.13 / 0.13 / 0.12) | 39, 20 (0.16, 0.39) | r131 root 37, arms 10, qtoy2 root 27, r2-r4 16-28, r131 r1 24 (0.23-0.48) |
| stage-B states | 41 / 38 / 33 (0.19 / 0.20 / 0.18) | 31, 11 (0.22, 0.51) | r131 root 6 (seed 1: 4), r1 13, arms 3, qtoy2 root 16, r2-r4 7-10, ole root 14, ole r1-r5 6-23 (0.26-0.87) |

At the node-mask and recB0 stages a real root (r131 at 30 / 37 %, arms at 41 %) sits where the nulls and the refused
qtoy2 r1 sit, so a stop there would cost real splits. At stage B the largest real value is ole r2 at 22.5 % (an accepted
split, top-5 % share 0.26) and the smallest diffuse value is qtoy2 r1 at 31.3 %, then micnull at 33.3 %. The rule is
coded at the split step (`advance`, case `S_STAGEB`, `diffuse_variance`): the stage-B half maps' variance through
`half_variance` and `share_core`, `MSHARE` = 0.7, leaf when the core covers more than `MDIFF` = 0.3 of the molecule
(one group, `decide` makes the multi-state leaf the split test would have made); otherwise `state_split` runs as before.
Both constants are set by three nulls against fifteen real nodes with a 7-point margin on the real side.

What it saves: the split-test job of the stopped node, which in the base arm costs 0.15 (qtoynull), 0.72 (olenull) and
0.22 (micnull) core-hours, 0.05 at qtoy2 r1; a node's other steps cost 6-25 core-hours. Nothing else changes on the
fixtures, because the nulls were already one leaf by the split test and qtoy2 r1 already a refused leaf. On the 692k tree
the rule cannot shorten the deep lineages: the deep splits of memory `null-guards-do-not-stop-692k-deep-splits` have
localised differences (the locality guard kept 11-12 of 12), which is the opposite of what this rule stops. The stage-B
diffuseness is the quantity the locality guard already reads per state pair (`L / L_ref`), so the stop confirms the prior
negative result rather than beating it.

Gate replay (m23 and m23d arms, section 4): the stop fired at the root of every null (m23 40.7 / 37.8 / 34.6 %, m23d
40.7 / 46.1 / 34.6 %) and at no node of the three qtoy2 trees, the three r131 trees or the OLE tree; each null ended as the
one leaf the split test makes of it. Saved: three split jobs, about 1.1 of the 36-48 core-hours the three null trees cost.
By the feature-freeze and science-before-speed rules a leaf rule that buys 1 % on a 7-point margin is not adopted: the
commit keeps the measurement as a log line at every split step (`variance_spread`, `MSHARE` = 0.7, no `MDIFF`), so the
692k tree reports the fraction at every node for free; the stop is the one-line `frac > 0.3` branch of the m23 build.

## 4. Gate replays

Items 2 and 3 were built together (`builds/m23`, 9005a8c77 + both changes, fast gate 15 / 15) and replayed with the base
arm's keys and seeds. The two changes are separable in the logs: item 3 acts only where a node logs
`stage-B variance diffuse` (it then skips `state_split`), item 2 everywhere a node has two or more recB0 states. The
coordinator then asked for the same replays under the rework session's site ranking (variant D, commit 6b4697c38:
one subspace test on the endpoints' cross-half Grams, pole strength = the weakest endpoint's hold-out z), because the
old ranking's margins come from four quarter replicates and flip on small mask changes; that commit was cherry-picked
(dd7fa0f45) and built as `builds/m23d`.

### 4.1 m23: items 2 + 3 under the old ranking

| fixture | base (9005a8c77) | m23 | verdict |
|---|---|---|---|
| qtoynull md0 | one leaf (split test refused) | one leaf: root stopped by item 3 at 40.7 % (limit 30), no split test | pass |
| olenull md0 | one leaf | one leaf: stopped at 37.8 % | pass |
| micnull md0 | one leaf | one leaf: stopped at 34.6 % | pass |
| r131 md1 seed 2 | root seeds from site 1 (margins 12.2 vs site 3 7.0); 25401 / 18933; r1 splits (T 0.154, z 4.31) into r1_s2 13332 = closed corner 0.99 / 0.51, r1_s1; r2 refused (T 0.125, z 3.43) | 4 sites (10f1 / 10f2 fold into 10); root seeds from site 3 (site 1 4.98, site 3 7.55); 24988 / 19346; r1 splits (T 0.195, z 5.29) into r1_s2 13205 = closed corner 0.99 / 0.52, r1_s1 11783; r2 refused (T 0.089, z 2.28) | corner kept; site order flipped |
| qtoy2 md1 seed 1 | 7 leaves, NMI 0.712; combo 18 (closed) F1 0.87, 101 0.71, 43 0.32 | 6 sites as before but the factor masks of components 1 and 2 differ; root 23101 / 16940 / 13252 / 6707 (base 22813 / 17574 / 12893 / 6720); 9 leaves, NMI 0.667; combo 18 F1 0.59 (its leaf 0.42 closed + 0.38 combo 22), 101 0.71, 43 0.54, 49 0.75 | worse at this seed |
| ole md1 | 10 leaves, largest 11152 (0.138); root K* 5; sites per node 5 / 2 / 6 / 2 / 3 / 4 | 10 leaves, largest 10725 (0.133); root K* 5 from site 3's seeds (margins site 1 8.5, site 2 13.5, site 3 37.1); sites per node 4 / 2 / 2 / 4 / 2 / 2; landscape against base (the rework's scorer, base root flex embedding): W1 0.002, extrema 0.001, coverage 0.48 / 0.80, against base's own null 0.002 and 0.65 / 0.51 (variant D against base: 0.003, 0.20 / 0.68) | span unchanged, 6 fewer sites |

The r131 root flip is the fragility of memory `site-strength-ranking-decides-r131-corner`: site 1 keeps the same pole in
both arms (2949 vs 2942 seeds, the closed population, mask Dice 0.953 between the arms) but its pole-rule margin fell from
(17.3, 12.2) to (5.0, 5.1), while site 3's pole (now 3151 seeds, seed state 8; was 9076, seed state 2; mask Dice 0.92)
rose to (30.9, 7.5). The corner was still delivered because the root partition hardly moved (24988 / 19346 against
25401 / 18933) and r1's re-test passed with more margin. The qtoy2 change runs through `site_states=yes`: the minority
poles of sites 3-6 seed extra root states, the SVD factor masks of components 1 and 2 differ from the 2-means', the root
K-way partition shifts, and r3 then separates combo 43 instead of the closed combo 18 from 22. The same fixture gives
NMI 0.707-0.764 and combo-18 F1 0.51-0.89 across the rework session's variants, so one seed cannot tell an effect from
the seeding lottery; a seed-2 pair (base and m23) is in section 4.3.

Single-draw caveat: every r131 verdict above rests on one equal-N subsample draw of the split test. The draw is seeded by
the constant `RND_SEED` in `state_mergeup`, a file of the rework session; that session's replay instrument (`build_df`:
`SIMPLE_SPLIT_SEED` added to the constant, the production test logged, then its draw-free single pass) was run on the
saved split projects of the corner nodes in section 4.4, where m23's r1 statistic above (0.195) has the standing of the
base's 0.154: one draw of a test whose spread at these nodes is as large as its margin.

### 4.2 m23d: items 2 + 3 under the variant D ranking

Reference: the rework session's variant D arm (`sub4_*`: 6b4697c38 on 9005a8c77, same keys and seeds, `math_reworks/gates`).
Site discovery runs before the ranking on the same maps, so m23d's sites are m23's.

| fixture | variant D alone (sub4) | m23d (D + items 2 and 3) | verdict |
|---|---|---|---|
| qtoynull / olenull / micnull md0 | one leaf each (split test, locality guard) | one leaf each: root stopped by item 3 at 40.7 / 46.1 / 34.6 % | pass |
| r131 md1 seed 2 | root 5 sites (5f1, 2, 5f2, 10f1, 10f2), poles kept at sites 1 and 3, site 1 first (5.73 vs 4.13); 25401 / 18933; r1 splits (T 0.154, z 4.31): closed corner 0.99 / 0.51; r2 refused (T 0.125, z 3.43) | root 4 sites (10f1 / 10f2 one factor), the one kept pole at site 1 (5.50); 25397 / 18937; r1 splits (T 0.147, z 4.17) into r1_s1 12536 = closed corner 1.00 / 0.49 and r1_s2 12861; r2 refused (T 0.064, z 1.13) | corner kept, site 1 first |
| r131 md1 seed 1 | root 3 sites (3f1, 1, 3f2), poles at sites 1 and 3, site 1 first (5.13 vs 3.78); 23319 / 21015; r1 refused (T 0.116, z 1.93); r2 splits (T 0.197, z 5.61) into 11601 / 9414 = closed corner (base: at r2, 0.97 / 0.54) | root 2 sites (3, 1: component 3's SVD has one factor over z 3, shares 0.75 / 0.19), the one kept pole at site 1 (5.11); 23449 / 20885; both children seed from their fflex states as in every arm ("no site gives a split"); r1 refused (T 0.181, z 3.27; its map s02 0.92 / 0.62); r2 refused (T 0.107, z 2.30, threshold 0.165) although its stage-B map holds the corner (s02 0.97 / 0.53) | corner refused by one draw; 4.4 |
| qtoy2 md1 seed 1 | root site 1 first (22.1, 15.0, 13.3, 11.5, 13.0, 12.7); K* 5: 22101 / 11471 / 10522 / 9289 / 6617; r2 and r4 seed from site 1 and split (7150 / 4321, 4955 / 4334); r1 seeds from site 7 (margin 3.92) -> 12700 / 9401; r3 (6 sites) and r5 (4) seed from their fflex states and are refused (T 0.165, 0.113); 8 leaves, NMI 0.764; F1 combo 18 0.89, 43 0.85, 49 0.83, 22 0.72, 28 0.94, 1 0.73, 2 0.59, 101 0.94 | root site 1 first (22.1, 15.0, 13.4, 13.3, 11.3, 11.4), the six sites with masks of Dice 0.92-1.00 to the 2-means'; K* 5: 22118 / 11473 / 10516 / 9272 / 6621; r2 and r4 as in sub4 (7155 / 4318, 4926 / 4346); r1 keeps no pole and seeds from its fflex states -> 12053 (combo 1 0.85) / 10065 (combo 2 0.73); r3 (4 sites) refused (T 0.140); r5 split (T 0.200, z 4.34) into 3618 / 3003, both combo 101; 9 leaves, NMI 0.796; F1 combo 18 0.89 (leaf 4346, 0.82 closed), 43 0.85, 49 0.83, 22 0.72, 28 0.94, 1 0.85, 2 0.80, 101 0.71 | same root; NMI up through r1's seeding, r5 over-split |

The r131 trees under D differ from the reference only where item 2 merges a factor pair (10f1 / 10f2 at seed 2, 3f1 / 3f2 at
seed 1): the root site order is site 1 first at both seeds as in the reference, the root partitions agree to within 130
particles, and the depth-1 nodes seed identically (fflex states). The qtoy2 root partition agrees to within 20 particles per
group, and six of the eight combo F1 values are identical; the two that differ are r1 (no kept pole under item 2, so the
fflex states seed instead of site 7's 3.9-margin pole: combos 1 and 2 go from 0.73 / 0.59 to 0.85 / 0.80) and r5 (T 0.200
against 0.113 at a node the rework session measured at 1 pass in 4 draws: a split-test draw that over-splits combo 101,
0.94 to 0.71). The NMI gain is the r1 gain net of the r5 loss.

### 4.3 qtoy2 at seed 2, base against m23

Same keys, seed 2, old ranking, `builds/base` against `builds/m23`.

| | base seed 2 | m23 seed 2 |
|---|---|---|
| root sites / kept poles | 6 sites; poles at sites 1, 2, 6 (104.6, 88.3, 72.8), site 1 first | 6 sites (2f1 / 1f2 order swapped, masks shifted); poles at 1, 2, 4, 6 (104.6, 88.3, 63.5, 30.9), site 1 first |
| root partition | 22954 / 17015 / 13317 / 6714 | 22981 / 17181 / 13129 / 6709 |
| leaves, NMI | 9, 0.715 | 7, 0.733 |
| combo F1 | 1 0.84, 2 0.80, 18 0.59, 22 0.49, 28 0.95, 43 0.53, 49 0.75, 101 0.65 (r4 over-split 3320 / 3394) | 1 0.84, 2 0.80, 18 0.59, 22 0.49, 28 0.94, 43 0.52, 49 0.73, 101 0.94 (r4 refused, T 0.132 at z 3.23) |

At seed 2 the arms agree on every closed-state number (combo 18 at 0.59 in both: r3 separates 18 from 22 in neither) and
differ at r2 (3 against 2 leaves) and r4 (an over-split of combo 101 in the base). With seed 1 (0.667 against 0.712) the
item's qtoy2 effect under the old ranking is within the seeding lottery; under variant D (4.2) it is +0.03 with the same root.

### 4.4 The split test over several draws at the corner nodes

The rework session's instrument (`math_reworks/build_df`, its `draws.sbatch` recipe) replayed on the saved split projects
(`mask_disc/draws`, 16 runs, 3 core-hours): the production test at `SIMPLE_SPLIT_SEED` 1-4 (the tree's own draw is seed 0)
and the draw-free statistic (`T_df`, section 8 of the rework note) on the same labels. Thresholds are each draw's own null at
z_crit 3.71 (m 460). The base seed-2 r1 and sub4 seed-1 r1 rows are the rework note's.

| node | own draw | T at seeds 1-4 (threshold) | passes | mean T / mean z | T_df (z over the seeds) |
|---|---|---|---|---|---|
| m23d r131 seed 2, r1 (the corner) | 0.147 P | 0.139 (0.139) P, 0.134 (0.160) F, 0.121 (0.128) F, 0.148 (0.153) F | 2 / 5 | 0.138 / 3.6 | 0.126 (7.4-11.0) |
| base r131 seed 2, r1 (the corner) | 0.154 P | 0.074 F, 0.137 F, 0.166 P, 0.097 F | 2 / 5 | 0.126 | 0.117 (9.2-10.9) |
| m23d r131 seed 1, r2 (the corner in the map) | 0.107 F | 0.151 (0.190) F, 0.105 (0.215) F, 0.096 (0.172) F, 0.164 (0.177) F | 0 / 5 | 0.125 / 2.3 | 0.122 (6.8-10.4) |
| sub4 r131 seed 1, r2 (the reference's corner) | 0.197 P | 0.160 (0.163) F, 0.108 (0.154) F, 0.088 (0.139) F, 0.118 (0.183) F | 1 / 5 | 0.134 / 3.1 | 0.121 (6.1-10.1) |
| m23d r131 seed 1, r1 | 0.181 F | 0.088 (0.207) F, 0.177 (0.166) P, 0.131 (0.239) F, 0.162 (0.195) F | 1 / 5 | 0.148 / 2.6 | 0.122 (9.0-10.7) |
| sub4 r131 seed 1, r1 | 0.116 F | 0.132, 0.168, 0.180, 0.117, all F | 0 / 5 | 0.143 | 0.119 (4.6-8.6) |

At every corner node the production verdict is one or two draws in five, in the reference trees as in m23d: the reference's
seed-1 corner (sub4 r2, T 0.197) is its one passing draw in five, the base's seed-2 corner (0.154) one of two. Item 2 moves
neither the mean single-draw statistic (0.138 against 0.126 at seed 2, 0.125 against 0.134 at seed 1) nor the draw-free
one (0.126 against 0.117, 0.122 against 0.121), and under the draw-free test every one of these nodes splits at z >= 6, so
both seeds deliver the corner in both arms. Read as the coordinator asked (average over draws before calling the corner lost
or kept): kept at seed 2; at seed 1 exactly as well supported as in the reference, and which single tree shows it is the
split test's draw.

## 5. Verdicts, commits and open questions

| item | fixture evidence | before -> after | gate verdict | cost (core-h) | go / no-go |
|---|---|---|---|---|---|
| 1 mask core by a share rule | 19 nodes' variance maps (prototype) | masks 93 -> 98-100 % (r131 root), 51 -> 56-66 (qtoy2 root), 25-60 -> 65-99 (depth 1) | not run: every mask grows | 4.5 (prototypes, all items) | no-go, no code |
| 2 SVD site factors | prototype: identical sites where the 2-means agreed, one fragment fewer at the r131 root, arms 5 -> 4 at the same recall, nulls unchanged. Gates: nulls 3 / 3 one leaf; r131 corner kept at seed 2 under both rankings and at seed 1 under D as well supported as the reference over five draws (4.4), root site 1 first at both seeds under D; qtoy2 NMI 0.796 vs 0.764 (D, seed 1), 0.733 vs 0.715 (old, seed 2), 0.667 vs 0.712 (old, seed 1); OLE 10 leaves as base, landscape within the base's null, 16 sites over six nodes instead of 22 | 2-means factor split (cos < 0.8, shares >= 0.1) -> SVD factors kept at share >= 0.2 and cross-half Fisher z >= 3 | pass: not worse on the corner, nulls and OLE; better on qtoy2 under D | m23 151, m23d 84, seed-2 pair 48, draws 3, builds 6 | go: committed as the default (413c0ce39), one commit to drop if the user disagrees |
| 3 diffuse-variance stop | stage-B share-0.7 core covers 33-46 % of the molecule on the nulls, <= 23 % on every real split; the stop fired at 3 / 3 null roots and at no real node in seven trees | split job skipped at a diffuse root (0.15-0.72 core-h each) | pass as harmless; saves about 1 % of a node, nothing on real trees | in the arms above | log line only: the fraction is reported at every split step, no leaf rule |

Commits on this branch (off `claude/tree-dist` 9005a8c77 plus the cherry-picked dd7fa0f45):
- 413c0ce39 hetero_tree: item 2 (`discover_sites`).
- item 3 as the log line (`variance_spread`, `share_core`, `MSHARE`), built as `builds/final` (fast gate 15 / 15; qtoynull md0
  smoke tree `gates/final_qtoynull`: the line printed at its root (40.7 %), the split test ran as before and the locality
  guard merged it to one leaf; the r131 and qtoy2 paths are the m23d build's, unchanged by the edit).
- this note and `scripts/cls_expansion_proto/tree_mask_discovery_stop/` (the prototype and its array driver).

Open questions for the user:
1. Item 2 lands as the default. Its corner evidence at r131 seed 1 is "as well supported as the reference over five draws",
   not "the tree shows the corner": under the production single-draw split test no arm shows it reliably at that seed
   (reference 1 / 5, m23d 0 / 5; the base's seed-2 corner 2 / 5). The rework session's draw-free statistic (their commit
   9ee5d5c72, available, not the default) makes every one of these nodes split at z >= 6 in every arm; without it the corner
   at ~20k-particle nodes stays a draw whatever the site rule. Keep 413c0ce39, or hold it until the split test is settled?
2. Item 3: keep the log line only (committed), or also the leaf rule (`frac > 0.3`, section 3) that saved three split jobs on
   the nulls and nothing on the real trees?
3. The qtoy2 gain under D is one node: r1's seeding (item 2 keeps no pole there, the fflex states seed; combos 1 and 2 F1
   0.73 / 0.59 -> 0.85 / 0.80), net of an r5 split-test draw that over-splits combo 101 (0.94 -> 0.71); a second qtoy2 seed
   under D was not run (the seed-2 pair ran under the old ranking).
