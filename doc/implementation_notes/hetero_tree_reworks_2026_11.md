# hetero_tree: the mathematical reworks (November 2026)

Branch `claude/tree-math-reworks` off `claude/tree-dist` (9005a8c77). Spec: sections 4, 6 and 7 of
`hetero_tree_algorithmic_audit_2026_11_01.md`, amended on 2026-11-01 (the junk screen stays, with two soft uses). Every
rework is a same-seed replay on the gate fixtures before any default changes: qtoy2 (maxdepth 1), the r131 closed
subtree (maxdepth 1, the closed corner must be retained), the nulls qtoynull / olenull / micnull (maxdepth 0, the largest
leaf must hold >= 80 % of the particles) and, for reworks 1 and 2, the OLE span. Replay runs live under
`/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/` on batch2; builds there follow the dist-build recipe (an
x6342 node, Release, ARCHOPT, fast gate), so the numerics match the baseline builds.

Gate keys are the recipe18d gate runs' (the latest accepted arms): qtoy2 `nparts=4 mskdiam=232 minpop=900 site_layer=yes
site_states=yes kway_maxdepth=2` at maxdepth 1; r131 `nparts=12 mskdiam=232 minpop=3000 site_layer=yes site_states=yes
kway_maxdepth=-1 split_ntests=460` at maxdepth 1, seed 2, from `ancestral/fixtures/r131`; the nulls at maxdepth 0 with their
recorded keys; all with `oneflex=yes sta0_npeaks=2 junk_screen=no` and `SIMPLE_SPLIT_CONTRAST=maps`, maps kept.

## 1. Endpoint reproducibility by one subspace test (audit 6.3)

(in progress)

### 1.1 What changed

`site_rank` keeps the four endpoint EM fits (the lo / hi tails of fflex directions 1 and 3) and drops every halves and
quarters refit. Each fit already yields the endpoint's cross-half deviation maps at N* (even and odd, inside the site
mask). The hold-out-band Gram matrices of these maps (even.even, odd.odd, even.odd; `dev_gram`) carry every decision:

- an endpoint is accepted when the cross-half correlation of its deviation is `HALVES_KSIG` (3) Fisher SEs above zero,
  with `n_eff` from `mask_neff` (the same n_eff the hold-out SE in `site_starts.txt` used);
- two accepted endpoints are one pole when they correlate positively and neither's residual after the other reproduces
  across halves (the principal-angle test, `gram_resid`: project the other's direction out in each half, correlate the
  residuals, Fisher z below 3); groups merge closest first (highest disattenuated correlation) when every cross pair is
  such a pair; the two ends of one motion are anti-correlated and stay distinct poles by construction;
- a pole whose deviation lies in the span of two or more other positively correlated poles is redundant and dropped;
- the rest (best pole posterior below 0.5) is ranked by one minus that posterior and scored once (no EM); it is a state
  when its deviation reproduces and is distinct from every pole it correlates with positively;
- the pole rule keeps a pole when its share is below 0.5 and its seed set, scored once by the pole's posterior, gives a
  deviation that reproduces and is distinct from the other poles' directions. `site_pole_rule.txt` keeps its columns:
  the two margins are now the Fisher z of reproducibility and of distinctness, so `strongest_pole` reads them unchanged.

Per site: 4 EM fits + 1 scoring pass for the rest + 1 per minority pole, against 4 x 7 + 6 + 6 x (kept poles) fits
before. `site_halves.txt`, `site_poles.txt`, `site_seeds.txt`, `site_starts.txt`, `site_neff.txt` keep their formats
(the controller's parsers are unchanged); `site_gram.txt` is new. The `rowmap` key is still accepted and now unused.

Prototype (`subspace_proto.py` on the Gram dumps of 32 saved sites, qtoy2 / qtoynull / olenull roots): on qtoy2 the
grouping reproduces the halves verdict where that verdict was decisive (sites 1, 2, 5) and keeps apart endpoints the
halves test merged at a disattenuated correlation of 0.85 with a residual that reproduces at z 8 (sites 3, 4); on the
nulls the Fisher-z test rejects every endpoint at n_eff below about 100 and accepts the same quality-driven deviations
the halves test accepted at high n_eff (those are left to the pole rule and the split test, as before). The eigenvalue
rank of the symmetrised cross-half correlation matrix has no power against near-collinear ends and is not used.

### 1.2 Same-input comparison on 32 saved sites

The halves build (9005a8c77 + the Gram dump, `build_dump`) and the subspace build (`build_sub`) re-ranked the same 32
saved sites (same project, mask, latent and seed): the qtoy2_keep_s1 root (6 sites), the qtoynull_s2 root (20) and the
olenull_s2 root (6), 16 threads each.

| sites | same `site_poles.txt` | kept poles old / new | wall old / new (s per site) | EM fits old / new |
|---|---|---|---|---|
| olenull (6) | 6 / 6 | 0 / 0 | 292-382 / 83-119 | 28 / 4 |
| qtoynull (20) | 16 / 20 (the 4 differ only in whether an empty rest line is written; 0 poles both ways) | 0 / 0 | 40-56 / 18-24 | 28 / 4 |
| qtoy2 (6) | 2 / 6 | 3 / 6 | 49-110 / 20-44 | 28 / 4 |

On both nulls no pole is kept by either test. On qtoy2 the two sites the halves test decided (sites 1, 2) come out
identical, seeds included. Sites 3 and 4, where the halves test merged everything into one pole and then rejected it on
share, now give a two-ended site (z1lo against z3hi, anti-correlated; the merged z1hi+z3lo group is rejected on share
and its rows join z3hi at posterior 0.5). Sites 5 and 6 keep the same main pole; site 5 additionally keeps a rest
state of 88 seeds (its deviation reproduces at z 3.3 and is distinct at z 10) that the halves accept-margin refused. A
rest or pole of a few dozen seeds is legal (the floor is PROB_MULTISTATE_MIN_POP + 1 = 6) but is a state refine3D cannot
refine; whether the tree is harmed by such states is what the gate replay measures.

## 3. Null-calibrated mask core (audit 6.1): prototype verdict

`mask_disc_proto.py` recomputes `falloff_mask`'s cross-half covariance from the root flex half maps (exact on the qtoy2
and qtoynull roots: correlation 1.000 / 0.999 with `mask_variance.mrc`) and calibrates it against noise surrogates built
from the half differences through the same low-pass, scaling and smoothing (a derangement pairing, per-map weights
restored). The result on three roots (qtoy2, qtoynull, the r131 fixture):

| root | peak v / sd_null | sphere with z > 3 | old core (0.2 x peak) covers | z-core covers |
|---|---|---|---|---|
| qtoy2 | 3500 | 80 % | 52 % of the molecule | 100 % |
| qtoynull | 640 | 99.7 % | 99 % | 100 % |
| r131 | 173 | 75 % | 95 % | 100 % |

The between-state variance of the flex states is above the half-map noise ceiling almost everywhere: the state maps
differ reproducibly by small amounts over the whole molecule (scale, B-factor and pose residuals that the least-squares
scaling does not remove), so significance does not localise anything. The 0.2-of-peak rule asks a different question,
where the between-state variance concentrates, and it is that question the mask needs. **No-go** for 6.1 as proposed;
the null calibration has no binding use here (as a floor it never binds). The dilation and falloff are unchanged.

The r131 fixture root takes the cheap-children path (10 maps, equal weights), not the flex maps; recomputed from those
the covariance matches `mask_variance.mrc` exactly (correlation 1.000) and the picture is the same: z > 3 on 34 % of
the sphere, the z-core one 692k-voxel blob covering the whole molecule against the 0.2-of-peak core's 94 %. What this
node adds is in section 4: none of its core components carries a reproducible factor (cross-half r 0.2-0.3 at n_eff
30-40), so the 93 % mask the sibling session flagged is the 0.2-of-peak rule applied to a covariance map that is mostly
noise. A node whose components reproduce nothing should get no focus (the whole molecule) rather than an arbitrary
blob; that veto is the one knob-free use of the calibration found here, and it belongs to discovery (section 4), not
to the mask threshold.

## 4. Cross-validated SVD for site factors (audit 6.2): prototype verdict

The same script factorises each core component's voxel-by-state deviation matrix (columns weighted by sqrt of the state
weight) by SVD of the even half and projects the odd half on the right singular vectors. Per factor it reports the
cross-half correlation of the factor's spatial map (Fisher z with n_eff from surrogate noise pairs on the component's
voxels, n_eff 7-37 at 15 A) and the singular value against a noise floor (the largest singular value of pure-noise
matrices from the half differences on the same voxels).

- Reproducibility over-counts: on the qtoy2 root the leading 2-4 factors of every component reproduce (r 0.9-0.999) and
  2-5 singular values sit above the noise floor, including factors holding 1-4 % of the component's variance. Small
  reproducible differences are everywhere, as in section 3.
- The SVD disagrees with the 2-means in both directions: component 1 (variance shares 0.93 / 0.04) is split by the
  2-means into "two factors" at 0.65 / 0.35 (equal-weight normalised voxel profiles), while components 3 and 4 with a
  clear second factor (0.77 / 0.20 and 0.76 / 0.21) are called one factor by the 2-means.
- Verdict: a **partial go**. The weighted SVD is the right factorisation (deterministic, gives each factor's state
  direction, spatial map and share) and should replace the 2-means, but the number of site factors cannot be set by
  reproducibility or by the noise floor alone; one variance-share threshold on the singular values (the DFACT_SHARE
  analogue, one constant instead of three) stays. Not ported to Fortran in this round; the factor masks of the SVD
  against the 2-means masks on qtoy2's known motions are the next prototype step.

## 2. The single-pass split test (audit 6.7)

### 2.1 What changed

`state_mergeup%test_once` replaces `test` in `state_split` and in the merge-up (`run`); `decide_once` replaces `decide`;
the old routines stay in the file until the rebase onto the cleanup branch that deletes the switched-off variants inside
them (deleting them here would conflict with that cleanup).

- One cross-half EM, started from the stage-B labels of every state present in both halves (no board decision in the
  coarse band first).
- The equal-N maps of every EM group are built once per half at the global `neq` (a group the EM left below `neq` is
  drawn whole, with its scale noted), docked onto the other half's consensus, and every pair statistic is read off them
  in the hold-out band (`kt`..`khi`). After a merge the merged group is redrawn once at `neq`; the old path rebuilt four
  maps per pair at the pair's own `min(|p|, |q|)` and re-tested merged groups at a larger N than the null's.
- One stratified shuffle null at the same `neq`, in the same band, one `z_crit` (`split_zcrit(ntests)`; the merge-up
  uses the largest `ntests` the tree's nodes used instead of a max-of-20 threshold); the null is still extended while
  a statistic sits within reach of the threshold.
- The node writes `split/split_stats.txt`: the weakest kept pair statistic, the weakest pair at the start, the null's
  mean and sd, `z_crit`, `neq`, K* and `ntests`. `equiv_tree_inputs` reads these numbers for the merge-up's margin
  (`split_effect`: the same `r_hat = 2y / ((1 - y) n)` conversion, on the weakest kept pair) and the equivalence test
  uses the null's sd directly instead of back-estimating it from the median and the maximum of a parsed log.
- The junk screen's SSNR unit screen is unchanged; `junk_report`'s reconstructions (a log line) are not run; the
  switched-off `route_members`, `em_null` and the floor are not part of the new path.
- The quality-split flag (the audit amendment's second soft use): after the test, `fit_contrast` gives each particle's
  contrast and zero-contrast posterior; `quality_flag` appends to `split_stats.txt` each delivered group's mean contrast
  and p_empty with its standardised offset `d` from the node mean (the within-group contrast sd as the unit) and a flag
  when any |d| >= 1. It gates nothing. The first soft use, the `1 - p_empty` particle weight in the site EM and the node
  reconstructions, is not in this commit: it changes every reconstruction and needs its own replay.

## 5. Map-space same-variable test and score-at-fixed-maps climbs (audit 6.4 / 6.5): prototype notes

6.5 cannot be prototyped on saved outputs: scoring a parent's particles against a child's two fixed pole maps is an
E-step over the particle planes. It is, however, already available in the code path `site_rank` uses: `site_em` with a
given score skips the EM and `nn_score` scores particles against fixed references, so a climb becomes one `site_rank`
call with the child's pole maps as references (one scoring pass per climb instead of a 4-fit EM plus the halves). The
memory record says the pooled sites drove no split, so the first question for the user is whether the climbs are wanted
at all before this is built. Not implemented in this round.

6.4 is prototyped below on a finished tree's factorial layer (pole maps, posteriors and refined masks of the kept sites).

`pole_angle_proto.py` on the closed_r11 fixture tree (r11_new9, three kept root sites with pole half maps; the refined
masks were pruned, so each site's footprint is the top 10 % of its own deviation magnitude): the pole deviation maps
(pole half map minus the particle-weighted leaf mean, 15 A, least-squares scaled) are reproducible at r 0.99 across
halves and pairwise at 63-88 degrees on the union footprints (55-77 degrees on the whole support), residual z 33-59 at
n_eff 134-503. The posterior Spearman of the same pairs is -0.03 to 0.05. Both tests call the three sites distinct, so
this tree has no fold to contrast. The decisive case is a tree where `FGATE` folded sites (the OLE deliverable trees);
the base_ole / sub_ole gate trees carry a factorial layer and are the next input for this script. Verdict so far: the
map-space test is cheap (the pole half maps already exist) and gives the angle and its noise ceiling directly; whether it
disagrees with the posterior Spearman where it matters is untested.

### 2.2 Null controls (maxdepth 0, same seeds and keys as the baseline)

| arm | qtoynull | olenull | micnull |
|---|---|---|---|
| baseline (9005a8c77) | 1 leaf (root K* 2, guard-merged) | 1 leaf | 1 leaf |
| rework 1 (subspace sites) | see 1.3 | | 1 leaf |
| rework 2 (single pass) | 1 leaf (K* 2 before the guard, weakest pair 0.083 vs threshold 0.059 from 30 shuffles) | 1 leaf (0.186 vs 0.072) | 1 leaf (0.795 vs 0.077) |

The single pass reproduces the baseline's root verdicts on the nulls: the re-test passes a quality / view split at the
root in every arm (as recorded for this fixture set since 2026-10-29) and the locality guard merges it, so the largest
leaf holds 100 % of the particles. The null thresholds agree with the baseline's to the second decimal (0.059 vs
0.061, 0.072 vs 0.067, 0.077 vs 0.077); the extension of the null while a pair sits near the threshold still runs
(30 of 50 shuffles on qtoynull). A bookkeeping bug in the first rework 2 build (the weakest kept pair statistic written
to `split_stats.txt` was overwritten by the threshold) affects only the merge-up's margin, not the splits; it is fixed
in the committed code and the rework 2 trees' merge-ups are re-run with the fixed binary.

### 1.3 Gate replay (same seeds and keys; baseline = 9005a8c77, `build_dump`; rework 1 = `build_sub`)

| gate | baseline | rework 1 |
|---|---|---|
| qtoynull md0 | 1 leaf (100 %) | 1 leaf (100 %) |
| olenull md0 | 1 leaf | 1 leaf |
| micnull md0 | 1 leaf | 1 leaf |
| qtoy2 md1 | 7 leaves, NMI(leaves, structure) 0.712; root 4-way | 8 leaves, NMI 0.764; root 5-way |
| r131 md1 (closed corner >= 0.95 / <= 0.6) | r1_s2 0.99 / 0.51 (corner) | variant A: lost at seed 2, delivered at seed 1; variant B: delivered at both seeds (below) |
| OLE md1 | 10 leaves (root 5-way) | variant A: 8 leaves (root 4-way); leaf landscapes as close as the null (extrema differ by 0.002 span units, W1 0.003 against the null's 0.002) |

qtoy2 best-leaf F1 per structure combo, baseline -> rework 1: combo 1 0.67 -> 0.73, 2 0.52 -> 0.59, 18 (closed) 0.87 -> 0.89,
22 0.54 -> 0.72, 28 0.94 -> 0.94, 43 0.32 -> 0.85, 49 0.73 -> 0.83, 101 0.71 -> 0.94. The root's seeding site is the same
(site 1, z1lo); the extra states come from the sites whose poles the subspace test kept where the halves test kept none
(site 3 seeds a state instead of site 5), and the root splits five ways instead of four. Quality NMI stays 0.00.

**r131, variant A (pole strength = the seed set's deviation z).** Both arms find the same two poles at the root (site 1:
z1lo+z3hi, 2981 seeds, the size of the closed population; site 3: z1lo+z3hi, 9076 seeds), with the same seed counts.
`strongest_pole` ranks the sites by the pole rule's margin: the halves margin put site 1 first (12.2 vs 7.0), the
seed-set Fisher z puts site 3 first (4.66 vs 3.29, three times the seeds give a more reproducible seed deviation). The
root therefore seeds from site 3, its partition shifts (23910 / 20424 against 25401 / 18933), node r1 still carries the
closed-corner state at stage B (1.02 / 0.49, baseline 0.99 / 0.51) but its re-test falls from z 4.31 to 2.45 against
3.71 and r1 becomes one leaf: the corner is not delivered at maxdepth 1. **Hard gate failed** for variant A.

The pole tests themselves did not change the outcome; the ranking statistic did. The audit names the hold-out peak as
the principled selection in `site_rank`, so variant B ranks a pole by its weakest endpoint's hold-out reproducibility
and distinctness (`zself`, `zdist` of the endpoint fits; site 1 5.8 vs site 3 5.1 here) and keeps the seed-set tests as
the keep criterion only; `site_pole_rule.txt` carries the strength in columns 4-5 and the seed-set z's in two new
columns before the verdict. Variant B is replayed on r131 (seeds 2 and 1) and qtoy2 (`sub2` arm); the nulls keep no
pole in either variant and are unchanged. Seed-1 replays of the baseline and variant A test whether the corner's
delivery at maxdepth 1 is itself seed-robust.

### 2.3 A variant without the choose routing (build_r2 = rework 1 variant A + the single pass, membership = the raw EM groups)

This build predates the port of the production choose rule into `test_once` and could not be cancelled, so it stands as
an ablation of the routing: qtoy2 root groups are mixtures of stage-B states (group 1: 160 / 21728 / 81 / 1059 / 845
particles of states 1-5) where the choose rule delivers pure states; r131 root 24393 / 19941 and r1 splits into r1_s1
(1.01 / 0.49, the closed corner) and r1_s2 (0.76 / 0.57), so the corner is delivered here although the root seeded from
site 3 as in variant A. Nulls one leaf each. Core-hours: nulls 6-14, qtoy2 22, r131 15, OLE 24.

**A defect found by the replay, and its fix.** At node r1 of the routing-less qtoy2 tree stage A refused to start:
"state 3, population 4, required minimum 6". The pole rule's re-seeding (unchanged from production) puts rows no kept
pole claims at posterior 0.5 into a rest bucket; with the new rest and pole tests accepting small states more readily
(a kept pole of 248 seeds, a rest bucket of 4) the tree now meets refine3D's floor. Fix (variant C, `build_sub3` /
`build_r2e`): a state out of `site_rank` needs at least SITE_GRID(1) = 500 seeds, the smallest N at which the hold-out
curve measures a deviation, below which `site_rank` has no evidence about it; `seed_poles` uses that floor for poles and
the rest, and `fold_small_states` folds any smaller final state (pole or rest bucket) into the remaining state with the
highest mean posterior for its rows, marking the pole "folded" in `site_pole_rule.txt` (the controller's parser ignores
it). refine3D's own floor (6) is unchanged.

### 2.4 Rework 2 against the subspace arm (same upstream, variant A sites; `build_r2c` = single pass + choose routing)

| gate | subspace arm, old split test | subspace arm, single pass |
|---|---|---|
| nulls (3) | 1 leaf each | 1 leaf each (root K* 2 before the guard in every arm) |
| qtoy2 md1 | 8 leaves, NMI 0.764 | 10 leaves, NMI 0.728: r3 and r5 split inside combos 28 and 101 (the halves differ in junk share, 0.06 vs 0.10 and 0.05 vs 0.12, not in the quality quintiles: NMI with quality 0.001) |
| r131 md1 | r1 refused (z 2.45), corner not delivered | r1 split (0.197 vs threshold 0.144), corner delivered as r1_s2 1.02 / 0.49 |

The three changed verdicts have one cause. After the cross-half EM one group sits far below the stage-B equal N (r131
r1: 2274 against 5877 per half; qtoy2 r3: 83 %, r5: 80 %). The old test drew both pair maps at the smaller EM group
(so a pair of 2274-particle maps) but kept its null at the stage-B N of 5877: the mismatch the audit names, biased towards
refusing. The first single pass drew the large group at the full N and the small one whole, so its pair was cleaner than
the null's pairs: biased the other way, towards passing (the qtoy2 within-combo splits). Neither is one N. The committed
single pass sets the equal N after the EM to the smallest EM group in either half and draws both the pair maps and the
null at it, with the map regularisation scaled to that N (`build_r2f`); its replay is in 2.5.

**Seed 1 (both arms, same keys).** The baseline's halves margin itself ranks site 3 first at seed 1 (5.14 against
3.64; at seed 2 it ranked site 1 first, 12.2 against 7.0) and delivers the corner through r2 (r2_s1, 0.97 / 0.54).
Variant A at seed 1 ranks site 1 first (5.02 against 4.47; the reverse of its seed-2 order) and delivers the corner
through r2 as well (r2_s1, 0.97 / 0.51). So the two root sites are a near-tie under both statistics and which one seeds
the root changes with the EM's random subsample; the corner's delivery at maxdepth 1 then depends on how the root
partition falls. Tally so far: baseline 2 of 2 seeds, variant A 1 of 2. The variant B and C arms add two seeds each.

**Variant B with the unequal-N single pass (`build_r2d`), seed 2.** r131: the root seeds from site 1 again and the tree
is the baseline's to the particle (root 25401 / 18933, r1 leaves 13332 / 12069, corner r1_s2 0.99 / 0.51, r2 one
leaf). qtoy2: identical to the variant A + single pass tree (10 leaves, NMI 0.728, r3 and r5 split by quality), so the
ranking change touches r131 only and the qtoy2 over-splits belong to the split test's N handling (2.4).

**Variant B alone (`build_sub2`).** r131 seed 2: site 1 strongest (5.73 against 4.13), the baseline's tree to the
particle, corner r1_s2 0.99 / 0.51. r131 seed 1: site 1 strongest again (5.13 against 3.78), corner r2_s1 0.97 / 0.51.
Two of two seeds, as the baseline; the endpoint hold-out z ranks the two sites the same way at both seeds where the
halves margin flipped. qtoy2 seed 1: identical to variant A (8 leaves, NMI 0.764; site 1 is strongest at the root
under either statistic). qtoy2 at seed 2 (no baseline counterpart): 6 leaves, NMI 0.710, a measure of the tree's own
seed variance.

`pole_angle_proto.py` on the base_ole factorial layer (three kept root sites, refined masks pruned as above): r_s2 against
r_s3 and r_s5 are distinct under both tests (Spearman -0.04 / 0.03; angles 90 / 82 degrees, residual z 23-27); r_s3
against r_s5 is the fold case: posterior Spearman -0.60 (folded by `FGATE`) and in map space a disattenuated cosine of
-0.81 (144 degrees) with a residual that reproduces at z 18, that is, the two ends of one motion axis with a 34 %
component the ends do not share. Both tests agree on the six pairs seen (two trees); the map test adds the angle and
the residual and would need its own threshold on the angle, so it is not knob-free either. Verdict: **partial**, untested
where it matters (a deliverable OLE tree with many folds); not pursued further in this round.

**Variant C (`build_sub3`: B + the 500-seed floor and fold).** r131: seed 2 the baseline's tree (corner r1_s2
0.99 / 0.51), seed 1 corner r2_s1 0.97 / 0.51; two of two. qtoy2 seed 1: 10 leaves, NMI 0.707 (baseline 0.712,
variant B 0.764): the root partition is variant B's to within two particles, but at the depth-1 nodes the floor folds
poles of a few hundred seeds that variant B kept as states, the children's stage A / B then settle on two states, and
r2, r3 and r5 split where variant B refused (r3 and r5 inside combos 28 and 101, by junk share rather than quality
quintile, F1 0.94 -> 0.65 and 0.94 -> 0.64; combo 43 0.85 -> 0.75). The crash in 2.3 needs only refine3D's floor of 6 applied to every state
`site_rank` writes (the rest bucket of the pole rule had none), not 500: variant D keeps B's ranking and folds below 6
(`build_sub4`); its replay is below. With the unequal-N single pass on top of C (`build_r2e`) the trees are identical to
C's to the particle on both qtoy2 and r131 (the choose rule delivered the stage-B states in every node).

### 2.5 The one-N single pass (`build_r2f`, on variant C's sites, seed 2)

| gate | variant C, old test | variant C, one-N single pass |
|---|---|---|
| nulls (3) | 1 leaf each | 1 leaf each |
| qtoy2 md1 | 10 leaves, NMI 0.707 | 9 leaves, NMI 0.728 (r3 refused, combo 28 whole again; r5 still split) |
| r131 md1 | corner r1_s2 0.99 / 0.51 | r1 refused (0.109 against 0.137), corner not delivered |

At r131's r1 both runs share the EM (MAP changed 0.344 / 0.325) and the equal N after it (1923 per half, from a
stage-B N of 5966). The null threshold is the same whether the null is drawn at 5966 (0.136) or at 1923 (0.137): the
null statistic hardly depends on N. What differs is the pair statistic of one random equal-N subsample of the same two
groups: 0.154 in the old run, 0.109 in the new one, a spread as large as the margin over the threshold. The verdict at
this node, in either test, is set by the draw, the "variance of one seeded random subsample" the audit lists; the N
mismatch it also lists moves the threshold by nothing here. The remedy is to average the pair statistic (and the null's)
over several equal-N draws, which neither test does yet.

**Rework 2 verdict: no-go as a default.** The single pass is implemented (one EM from the stage-B labels, maps built once,
one N for pairs and null, one z_crit, numeric `split_stats.txt` read by the merge-up, the choose routing retained, the
quality-split flag) and passes the nulls, but it loses the r131 corner at seed 2 on a node the old test passes by one
subsample draw, and on qtoy2 it trades one within-combo over-split for another. The code is committed so the user can judge
the diff; the production path should not switch to it until the pair statistic is averaged over draws and re-gated.

## 6. Overnight items from the coordinator (2026-11-02)

### 6.1 Item 7: one ranking quantity for every site

The quantity asked for, the held-out evidence of two states over one inside the site mask, is what `site_rank`
already measures per endpoint: the hold-out curve's peak is the cross-half correlation, in the hold-out band inside the
mask, of the deviation of the top-N* particles of the endpoint's class from the node mean, and its Fisher z over
`mask_neff`'s n_eff is that evidence as a z. Variant B / D ranks a pole by the smallest such z over its endpoints
(`site_pole_rule.txt` columns 4-5, read unchanged by `strongest_pole`), so one quantity ranks every site and every
pole, the same one that accepts an endpoint; the seed set's own deviation z is the keep criterion only. Variant A
(the seed set's deviation z as the strength) and the halves margin (self-match over an SE from four quarter replicates)
are the two it replaces. Root site order beside the corner verdict, r131 fixture:

| arm | seed | root site order (strength) | corner at maxdepth 1 |
|---|---|---|---|
| baseline, halves margin | 2 | site 1 (12.2) > site 3 (7.0) | r1_s2 0.99 / 0.51 |
| baseline, halves margin | 1 | site 3 (5.14) > site 1 (3.64) | r2_s1 0.97 / 0.54 |
| variant A, seed-set z | 2 | site 3 (4.66) > site 1 (3.29) | lost (r1 refused, z 2.45) |
| variant A, seed-set z | 1 | site 1 (5.02) > site 3 (4.47) | r2_s1 0.97 / 0.51 |
| variant B / C, endpoint hold-out z | 2 | site 1 (5.73) > site 3 (4.13) | r1_s2 0.99 / 0.51 |
| variant B / C, endpoint hold-out z | 1 | site 1 (5.13) > site 3 (3.78) | r2_s1 0.97 / 0.51 |

Site 1's pole holds 2981 seeds, the size of the closed population; the endpoint hold-out z ranks it first at both
seeds. Variant D (B with the fold below refine3D's floor of 6 instead of 500) is the chosen form; its replay table is
in 1.3 and 6.3.

### 6.2 Item 13: quality as a nuisance

(a) The site EM already weights every particle's planes by `con` from `fit_contrast`, the posterior-mean contrast
(1 - p_empty) x E[amplitude] (`site_em`: `cn = self%con`), so the soft weight in the site EM is the existing behaviour;
nothing is removed. For the node reconstructions `state_split` now writes `split/ptcl_weights.txt` (project row,
1 - p_empty, contrast) for the controller to apply at the merge, which this session does not edit. On the fixtures
p_empty is 0.0000 for every delivered group, so the weight reduces to the amplitude term.

(b) `quality_flag` computes, for the split the test proposed (before the locality guard), each group's mean contrast and
mean p_empty, the standardised offset d of its contrast from the node mean (unit: the pooled within-group contrast sd)
and a flag at |d| >= 1, appended to `split/split_stats.txt` as numbers; nothing reads it as a gate. The first version
ran after the guard and so never saw the one-conformer null's root split (merged to K* 1 by the guard); the standalone
split replays below run it before the guard.

Standalone `state_split` replays with the flag computed before the guard (`build_r2h`, the base trees' split projects):

| node | split proposed (K* before / after the guard) | group contrast means | d | flag |
|---|---|---|---|---|
| qtoynull root (the one-conformer null's quality split; guard merges at L / L_ref 0.225) | 2 / 1 | 0.9897, 0.9893 | +0.002, -0.002 | no |
| olenull root | 2 / 1 | 0.9952, 0.9895 | +0.05, -0.08 | no |
| r131 root (structural) | 2 / 2 | 0.9887, 0.9825 | +0.03, -0.04 | no |
| qtoy2 root (structural, 4 groups) | 4 / 4 | 1.0015 .. 0.9813 | +0.20 .. -0.19 | no |

p_empty is 0.0000 for every particle on every fixture. The contrast of `fit_contrast` is each particle's amplitude
against the other half's mean map at its own noise level, so it is normalised to the particle's SNR and does not
respond to the quality (SNR / B-factor) spread that the null's split follows; the flag built on it cannot light up on
the fixture it must light up on. The split test already holds two per-particle terms that do see quality: the noise
variance it estimates for every particle when it loads the planes and the stage-B correlation. The next replay below
adds both to the flag.

With the three terms (`build_r2i`, standalone replays before the guard) no node flags: qtoynull root d = 0.00 / 0.01 /
0.03 (contrast / log noise variance / stage-B corr), olenull 0.05-0.08 / 0.01 / 0.05-0.07, micnull 0.05-0.07 / 0.01 /
0.03-0.04, r131 root and qtoy2 root and r1 all below 0.25. The simulation truth explains why: the qtoynull root's
pre-guard groups (31139 / 28861) have snrq 1.074 and 1.050 (d = +0.03 / -0.03 at a within-group sd of 0.38) and B-factors
100.4 and 98.9 (d = +0.01 / -0.01), so that split is not a quality split; it is a false positive of another kind (view or
featureless) and the locality guard merges it at L / L_ref 0.225. The stage-B correlation is a sound quality proxy
(r = 0.66-0.70 with snrq on both toys; contrast 0.22-0.55; the log noise variance 0.00, it is a whitening scale), so
the flag would respond to a genuine quality partition. qtoy2's root groups differ in junk share (d +0.34 for group 3)
and not in quality (|d| <= 0.02).

**Item 13 verdict.** (a) The soft weight is already in the site EM through `con`; the node-reconstruction weights are
written per particle (`split/ptcl_weights.txt`) for the controller's merge; p_empty is 0 on every fixture, so there is
no RNase P gate to run yet (junk_screen=no behaviour unchanged). (b) The flag is implemented as numbers in
`split/split_stats.txt` (per group: n, mean contrast, mean log noise variance, mean stage-B corr, their d's, mean
p_empty; a yes / no at |d| >= 1), computed on the split the test proposed and gating nothing. It does not light up on
the one-conformer null because that null's split is not a quality split; a fixture with a genuine quality partition is
needed to see it fire. The memory note that called the qtoynull root split a quality split is corrected.

### 6.3 Variant D (`build_sub4`): the committed rework 1

| gate | baseline (9005a8c77) | variant D |
|---|---|---|
| qtoynull / olenull / micnull md0 | 1 leaf each | 1 leaf each (no pole is kept on any null, so D = A there) |
| qtoy2 md1, seed 1 | 7 leaves, NMI 0.712 | 8 leaves, NMI 0.764; best-leaf F1 better on 7 of 8 combos (22: 0.54 -> 0.72, 43: 0.32 -> 0.85, 101: 0.71 -> 0.94), equal on one |
| r131 md1, seed 2 | site 1 > site 3; corner r1_s2 0.99 / 0.51 | site 1 (5.73) > site 3 (4.13); the baseline's tree to the particle, corner r1_s2 0.99 / 0.51 |
| r131 md1, seed 1 | site 3 > site 1; corner r2_s1 0.97 / 0.54 | site 1 (5.13) > site 3 (3.78); corner r2_s1 0.97 / 0.51 |
| OLE md1 | 10 leaves | variant A: 8 leaves, landscapes as close as the null; variant D's tree in 6.4 |
| site_rank cost | 28 EM fits per site; 49-382 s | 4 EM fits + 1-3 scoring passes; 20-119 s (2.2-3.3x less) |

**Rework 1 verdict: go**, as variant D (the subspace test, the endpoint hold-out z as the pole strength, refine3D's floor
applied to every state `site_rank` writes). The one default it changes is the site test itself, which the gates
support; the user decides whether the ranking statistic's seed stability on r131 (two of two seeds, where the halves
margin flipped between seeds) is evidence enough.

With the one-N single pass on top of D (`build_r2g`): nulls one leaf each, qtoy2 9 leaves at NMI 0.742 (combo 28 split
again, 0.94 -> 0.69), r131 r1 refused at seed 2 (corner lost): the same verdicts as on variant C's upstream, so rework
2's no-go does not depend on which site test feeds it.
