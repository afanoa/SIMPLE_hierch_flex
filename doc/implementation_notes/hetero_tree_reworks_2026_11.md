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

At seed 1 (variant C's sites, `build_r2f`) the one-N single pass refuses r2 as well, the node whose stage-B state is
the corner (0.97 / 0.51), so it delivers the corner at neither seed where the production test delivers it at both.

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
| OLE md1 | 10 leaves | 8 leaves (identical for variants A, C, D), landscapes as close as the null (6.4) |
| site_rank cost | 28 EM fits per site; 49-382 s | 4 EM fits + 1-3 scoring passes; 20-119 s (2.2-3.3x less) |

**Rework 1 verdict: go**, as variant D (the subspace test, the endpoint hold-out z as the pole strength, refine3D's floor
applied to every state `site_rank` writes). The one default it changes is the site test itself, which the gates
support; the user decides whether the ranking statistic's seed stability on r131 (two of two seeds, where the halves
margin flipped between seeds) is evidence enough.

With the one-N single pass on top of D (`build_r2g`): nulls one leaf each, qtoy2 9 leaves at NMI 0.742 (combo 28 split
again, 0.94 -> 0.69), r131 r1 refused at seed 2 (corner lost): the same verdicts as on variant C's upstream, so rework
2's no-go does not depend on which site test feeds it.

## 7. Summary (2026-11-02)

| rework | state | verdict |
|---|---|---|
| 1. subspace site test (6.3) | committed 6b4697c38 (variant D) | **go**: nulls one leaf each, qtoy2 NMI 0.764 against 0.712, r131 corner at both seeds as the baseline, OLE within the null; site_rank 2.2-3.3x faster per site |
| 2. single-pass split test (6.7) | committed ee764dda6, not the default | **no-go as default**: nulls pass, qtoy2 within-combo splits move, r131 corner lost at both seeds (two upstreams at seed 2); the verdict is set by the one equal-N subsample draw, so average the pair statistic over draws and re-gate |
| 3. null-calibrated mask core (6.1) | prototype only | **no-go**: between-state variance is above the half-map noise on 34-100 % of the sphere on every root, null included; the 0.2-of-peak rule stays |
| 4. cross-validated SVD factors (6.2) | prototype only | **partial**: the weighted SVD should replace the 2-means, but reproducibility alone over-counts factors; one share threshold stays |
| 5. map-space same-variable test, scored climbs (6.4 / 6.5) | prototype / design | **partial / design**: agrees with the posterior Spearman on 6 pairs (one fold case), needs an angle threshold; climbs as one scoring pass are available through site_em but untested |
| item 7 (one ranking quantity) | in 6b4697c38 | the weakest endpoint's hold-out z, the quantity that accepts an endpoint, ranks every pole and site; r131 order and corner per arm in 6.1 |
| item 13 (quality as a nuisance) | in ee764dda6 | weights: in the site EM already (con); per-particle weights written for the node reconstructions; flag: numbers in split_stats.txt, gating nothing; it stays off on the one-conformer null because that split is not a quality split (truth d 0.03) |

Cost of the replay campaign (sacct, 2026-11-01 18:00 to 2026-11-02 00:30): about 1350 core-h over 13 arms, of which
the baseline 178, the rework 1 arms (sub, sub2, sub3, sub4) 430, the rework 2 arms (r2, r2c, r2d, r2e, r2f, r2g) 590,
site and split replays, prototypes and builds 150. The routing-less rework 2 arm (128 core-h) ran because a cancel
was refused by the session's permission layer.

Not done: the p_empty weight in the node reconstructions needs the controller to read `split/ptcl_weights.txt` (the
controller is another session's tonight); the averaged pair statistic for the split test; the SVD factor masks against
the 2-means masks; the map-angle test on a deliverable OLE tree; the deletion of the old `test` / `decide` /
`eqn_stats` / `node_effect` path, which waits for the rebase onto the cleanup branch.

### 6.4 OLE at maxdepth 1 (80,599 particles, seed 1, `landscape_score.py --leaves` on the baseline root's flex embedding)

| arm | leaves | root | against the baseline: extrema max / mean, W1, hull coverage | against its own null: W1 |
|---|---|---|---|---|
| baseline | 10 | 5-way | - | 0.002 |
| variant A | 8 | 4-way | 0.002 / 0.001, 0.003, 0.20 / 0.68 | 0.002 |
| variant C | 8 | 4-way (identical to A) | 0.002 / 0.001, 0.003, 0.20 / 0.68 | 0.002 |
| variant C + one-N single pass | 8 | 4-way (identical to A) | 0.002 / 0.001, 0.003, 0.20 / 0.68 | 0.002 |

| variant D (committed rework 1) | 8 | 4-way (identical to A) | 0.002 / 0.001, 0.003, 0.20 / 0.68 | 0.002 |

The four rework arms give the same OLE tree to the particle, and all differ from the baseline by no more than the
baseline differs from its own shuffled null, so the OLE span is retained. At maxdepth 1 the leaf means barely spread
along the leading flex axes (ends at 0.48-0.56 of the span), so this gate has little power here; a deliverable OLE tree
at maxdepth 4 is where span differences would show.

## 8. The draw-free pair statistic (overnight item, 2026-11-02)

### 8.1 Definition (commit 9ee5d5c72, `pair_stat_df`, `noise_power`, `null_upto_df`)

Every group's maps are built from all its particles in each half (merged groups are the particle-weighted means of
their members, no redraw). Per shell k, the between-group signal is the cross-half product of the full-N difference maps,
C_k = <A_p - A_q, B_p - B_q>_k (the halves' noise is independent and cancels), and each group's per-particle noise power
is its half difference over the half sizes, sigma_k^2 = |A_g - B_g|_k^2 / (1/nA + 1/nB). The statistic is the expected
cross-half FSC of a pair of n-particle difference maps, FSC_k(n) = C_k / (max(C_k, 0) + (sigma_p,k^2 + sigma_q,k^2) / n),
3-shell smoothed and peaked over the hold-out band as `pair_stat` does, with n the smallest EM group in either half (the
same n the one-N pass used). The null is the same quantity from within-strata label shuffles at full N, stored whole and
evaluated at the same n. Nothing in the statistic or the null is drawn, so its only randomness is the shuffles'.

### 8.2 Python on the saved stage-B half maps (`drawfree_proto.py`, box 88 cropped to the test's 48, node mask, per-half state counts from the split log)

| node | states per half (A / B) | n | draw-free T at n/4, n/2, n, 2n | actual full-N cross-half FSC peak | old one-draw T (log) |
|---|---|---|---|---|---|
| r131 r1, seed 2 | 6103+6712 / 5966+6620 | 1923 | 0.121 0.215 0.350 0.513 | 0.671 | 0.154 |
| r131 r1, seed 1 (baseline) | 5231+5931 / 5147+6050 | 2878 | 0.234 0.372 0.534 0.689 | 0.718 | 0.210 |
| r131 r1, seed 1 (variant D) | 5533+6148 / 5497+6141 | 2718 | 0.244 0.387 0.551 0.705 | 0.761 | 0.116 |
| qtoy2 r3 (seed 1, D) | 2434+2885 / 2358+2845 | 1967 | 0.164 0.280 0.434 0.602 | 0.515 | 0.165 |
| qtoy2 r5 (seed 1, D) | 1671+1703 / 1588+1655 | 1268 | 0.089 0.164 0.281 0.439 | 0.382 | 0.113 |
| qtoynull root | 15643+14357 / 15496+14504 | 6828 | 0.079 0.146 0.252 0.398 | 0.420 | 0.089 |
| olenull root | 16544+23756 / 16524+23775 | 7992 | 0.215 0.350 0.512 0.672 | 0.722 | 0.088 |
| micnull root | 12139+16265 / 12184+16220 | 12580 | 0.655 0.786 0.878 0.934 | 0.888 | 0.801 |

The construction behaves as intended: deterministic given the labels, monotone in n, and at full N it reproduces the
actual cross-half FSC of the stage-B states. Its values on the stage-B states are not the Fortran test's values, which
are read on the cross-half EM's groups (a third of the particles move, and at r131 r1 the EM leaves one group at a
third of the stage-B size), so the comparison with the old statistic's draw spread is made in Fortran below, on the
same EM groups, by an instrumented replay build that logs the production test and then runs the single pass.

### 8.3 Draw variance: the production statistic over four subsample seeds against the draw-free one (same node, same labels)

Standalone `state_split` replays on the saved split projects with an instrumented build (`build_df`: the production test
runs and is logged, then the single pass with the draw-free statistic decides; the subsample seed comes from
`SIMPLE_SPLIT_SEED`, a replay-only instrument that is not committed). "old T" is the production re-test's weakest pair,
"T_df" the draw-free weakest pair at the start; thresholds are each test's own null at z_crit.

| node | old T over seeds 1-4 (tree's own draw) | old verdicts | T_df (every seed) | T_df threshold | z_df | df verdicts |
|---|---|---|---|---|---|---|
| r131 r1, seed-2 tree (the corner node) | 0.074 0.137 0.166 0.097 (0.154) | F F P F (P) | 0.117 | 0.045-0.054 | 9.2-10.9 | P P P P |
| r131 r1, seed-1 baseline tree | 0.203 0.205 0.225 0.221 (0.210) | P P P P | 0.189 | 0.080-0.115 | 6.6-9.4 | P P P P |
| r131 r1, seed-1 variant D tree | 0.132 0.168 0.180 0.117 (0.116) | F F F F | 0.119 | 0.059-0.101 | 4.6-8.6 | P P P P |
| qtoy2 r3 (seed 1, D) | 0.189 0.202 0.143 0.203 (0.165) | F P F P (F) | 0.182 | 0.114-0.135 | 4.9-6.5 | P P P P |
| qtoy2 r5 (seed 1, D) | 0.120 0.137 0.153 0.174 (0.113) | F F F P (F) | 0.133 | 0.068-0.102 | 4.8-7.3 | P P P P |
| qtoynull root | 0.062 0.059 0.065 0.055 (0.089) | P P P F (P) | 0.054 | 0.023-0.028 | 7.8-9.6 | P P P P |

Three facts. The old statistic's spread over the draw at these nodes (sd 0.02-0.04) is as large as its margin, so its
verdict at the corner node is one draw in four; the draw-free statistic has no spread, by construction, and lies at the
centre of the old one's distribution (0.117 against a mean of 0.119 at the corner node), as the expectation of a
single-draw FSC should. Its null is about three times tighter (sd 0.01 against 0.03): the signal estimate comes from the
full-N maps and only the analytic noise term is at the common n, so the statistic's own noise is that of the full maps,
in the null as in the real pairs. Every marginal pair the old test passed or failed by the draw is therefore clearly
significant under the draw-free test, including the one-conformer null's root split, which both tests pass and the
locality guard merges. The gate trees in 8.4 show what the extra power does to the leaves.

### 8.4 Gate trees on the draw-free single pass (`build_df`, same seeds and keys as variant D, which is its upstream)

| fixture | variant D (production test) | draw-free single pass | gate |
|---|---|---|---|
| qtoynull md0 | 1 leaf | 1 leaf: root T 0.054, z 7.8 PASS, guard merges (L / L_ref diffuse) | pass |
| olenull md0 | 1 leaf | 1 leaf | pass |
| micnull md0 | 1 leaf | 1 leaf: root T 0.80, z 82 PASS, guard merges (L / L_ref 0.24) | pass |
| qtoy2 md1, seed 1 | 8 leaves, 5-way root, NMI 0.764 | 9 leaves, the same root and r1 / r2 / r4 splits, plus r3 (combo 28, 4792 / 5730); NMI 0.742 | baseline 0.712; below D |
| r131 md1, seed 2 | 3 leaves; corner at r1 (0.99 / 0.51) | 4 leaves; corner at r1 (0.99 / 0.51, z 10.9) plus r2 (9818 / 9115; maps 0.77 / 0.64 and 0.93 / 0.58) | corner retained |
| r131 md1, seed 1 | 3 leaves; corner at r2 (0.97 / 0.51) | 4 leaves; corner at r2 (0.97 / 0.51, z 7.1) plus r1 (12289 / 11030; maps 0.93 / 0.63 and 0.79 / 0.60) | corner retained |
| OLE md1, seed 1 | 8 leaves, 4-way root, within the baseline's null | variant D's tree to the particle: 8 leaves, 4-way root, W1 0.000 against D (every split at z 8.5-28); against the baseline W1 0.003, extrema 0.002, the baseline's own null 0.002 | pass |

Cost: 114 core-h for the arm against 110 for the baseline arm, with the instrumented build running both tests at
every node, so neither the core-hours nor the wall times (r131 59 against 71 min, qtoy2 49 against 44) measure the
single pass alone.

The draw-free trees agree with variant D's at every node the production test passed and add exactly one split per
tree, each at a node the production test rejected by its draw: qtoy2 r3 (old T 0.165 at z 3.25 against z_crit 3.38;
draw-free 0.182 at z 5.6 against a threshold of 0.126), r131 r2 at seed 2 (old 0.125 at z 3.43 against 3.71; draw-free
0.056 at z 4.6) and r131 r1 at seed 1 (old 0.116 at z 1.93; draw-free 0.119 at z 5.6). None of the three is marginal
under the new test, whose null is a third as wide. What they split differs. On qtoy2 the extra split lies inside one
structure combo (28, both leaves 0.87-0.92 pure), is localized (L / L_ref 1.15) and is not a quality split by the flag
(|d| <= 0.1 on contrast, noise variance and stage-B correlation): a reproducible difference the truth table does not
name, and the reason the NMI falls from 0.764 to 0.742, still above the baseline's 0.712. On r131 both extra splits
produce a state at 0.93 closed-arm correlation, 0.02 under the closed-corner threshold, beside the retained corner at
0.97-0.99; the fixture cannot say whether that is a second closed-like substate or an alignment family.

Verdict. The statistic does what it was built for: the verdict no longer depends on a subsample draw (8.3 shows the
production test's retention of the r131 corner at seed 2 to be one draw in four), the nulls stay one leaf, the closed
corner is retained at both seeds with z 7-11, and the OLE tree is variant D's to the particle (no node there is marginal under either test). Its one cost is power: with a null a
third as wide it accepts one more marginal split per tree, and on qtoy2 that split is within a conformer. The gates
therefore support the statistic but do not clearly support it as the default on the NMI alone, so the default stays on
the production test (commit 680bda2d5 puts `exec_state_split` back on `split%test`; `test_once` with the draw-free
statistic stays in `simple_state_mergeup`, one call away) and the switch is the user's decision. If it is taken, the
equivalence merge-up is where the within-conformer split would have to be caught: at full N the two combo-28 leaves
differ reproducibly, so a merge-up that fuses equivalent maps leaves them apart (the df tree's merge-up kept 9 units).

## 9. The rare-state lever and the seed anomaly (2026-11-02, session tree-lever-floor, branch `claude/tree-lever-floor`)

Spec: the rare-state note `tree_rare_state_limit_2026_11.md` (sections 2.3 and 5) and audit 6.5. Builds and replays under
`math_reworks/lever/` on batch2 (`build_lever` = climb + complement variant, `build_bc` = the committed form + the effect
floor of section 10, `build_head` = c7ac3b0bd as the baseline arm); saved-site replays in `lever/sites/`, scored against the
fixture truth by `lever/site_score.py`.

### 9.1 The seed anomaly: what happens at the 8 % root

`site_starts.txt` of `detlim2/runs/new_s1/P_L08/node_r/site` (the arm site, n_eff 217): the two rare-ward endpoints have
a small, pure class at posterior >= 0.5 (z1lo 1,668 particles at rare purity 0.78, z3hi 2,088 at 0.74; their mean posterior
over the node, the "site class share", is 0.125 / 0.138) and N* 40,000; z3lo's EM grew its class to 0.834 of the node
(59,756 particles at posterior >= 0.5, rare purity 0.079 = the node's) and its hold-out peak sits at N* 10,000 (the
10,000 most typical particles, r 0.273). The subspace test accepts all three (z 4.0-4.3). `seed_poles` then gives every
particle to the pole with the highest posterior at >= 0.5: z3lo's flat 0.9 wins against the rare pole's 0.6 almost
everywhere, so pole 1 = {z1lo, z3hi} draws 717 seeds (the few where its mean posterior exceeds z3lo's) and pole 2 = {z3lo}
59,145. The pole rule then rejects both (seed z 1.36; share 0.99): no pole, the fflex fallback seeds. Seed 2 is the same
with the roles of z3hi and z1hi exchanged (545 / 59,368). The halves test on the same fits had merged z3lo and z3hi
into one pole (mean posterior of a 0.83 class and a 0.14 class >= 0.5 for 18,600 particles) by accident of the averaging.

A class that grows past half the node has not found a minority state; the pole rule says so (share < 0.5) but only after
the class has claimed the real poles' particles. The nulls show the same endpoints as the quality axis (qtoynull z1hi /
z3lo at shares 0.96 / 0.97, olenull z1lo / z3lo at 0.91 / 0.93, all accepted at z 10-27) and are clean today only because
the share rule rejects them at the end.

Three forms were replayed on the saved sites (`lever/sites/c0_*`: the complement; `c1-c3`: the complement plus 1-3 climb
passes; `bc_*`: the committed form):

- **the complement** (the pole of a took-off endpoint is its minority class, 1 - posterior, before the hold-out curve):
  at 8 % z3lo's complement reproduces (r 0.42, z 6.6; rare AUC 0.57) and pole 1 grows to 1,845 seeds (seed z 2.27, still
  rejected); at the r131 seed-2 root site the complement pole {z1hi', z3lo'} (4,566 seeds) is kept and the closed-population
  pole {z1lo, z3hi} falls to 349 seeds and is rejected: the hard-gate site changes hands. Rejected.
- **the committed form**: an endpoint whose class share exceeds 0.5 is no pole (`merge_poles_subspace` takes it as not
  accepted; `exec_site_rank` logs "its site class grew past half the node"); and when the rest is rejected as a state it
  joins its best pole only while it is a minority, a majority rest stays the complement of the poles (before, the took-off
  pole was the rest's backstop: rejected by share at the end, its rows fell into the rest). The pole rule no longer treats a
  rest without endpoints as a pole. Every gate site keeps its decision (the took-off endpoints were rejected by share
  anyway); the 8 % pole becomes {z1lo, z3hi} with its own seeds.

### 9.2 The climb (audit 6.5) does not lift the endpoint: no-go

`site_em` was given a climb: after the EM, one hard E-step from the top-N* of the posterior (class maps from that set and
its complement, the per-shell scale and noise refit, the prior at N* / n), then the hold-out curve again; 1-3 passes
(`SIMPLE_SITE_CLIMB` on `build_lever`, removed). The bar at site 5 is r 0.275 (3 Fisher SE at n_eff 116.5 / 112.4).

| site (truth) | endpoint | r at 0 / 1 / 2 / 3 passes (z at 0 / 1) | rare AUC at 0 / 1 | note |
|---|---|---|---|---|
| 6.5 % s1, site 5 | z1hi | .141 / .144 / .138 / .117 (1.5 / 1.6) | .754 / .761 | class at >= 0.5: 0 -> 21,209 particles (share .42) |
| | z3lo | .130 / .120 / .175 / .148 (1.4 / 1.3) | .780 / .789 | |
| 6.5 % s2, site 5 | z1lo | .158 / .206 / .193 / .152 (1.7 / 2.2) | .268 / .733 | the pass turns the anti-rare tail into the best rare endpoint, z 2.2 |
| | z3lo | .133 / .174 / .149 / .161 (1.4 / 1.8) | .787 / .821 | class 609 at purity .81 after one pass |
| 6.5 % s1, site 3 (nuisance) | z3lo | .697 / .502 / .459 / .428 | .48 / .52 | poles 50 / 49 seeds rejected -> one pass: a kept pole of 22,937 seeds (rare purity .069) |
| 8 % s1, site 1 | z1lo | .285 / .212 / .210 / .273 (4.3 / 3.2) | .889 / .160 | one pass: kept pole {z1lo, z3lo} 18,990 seeds at rare purity .028 (the wrong end); 3 passes: 21,980 at .165 |
| | z3hi | .268 / .199 / .276 / .214 (4.0 / 3.0) | .893 / .152 | |
| 8 % s2, site 1 | z1hi | .297 / .214 / .250 / .201 (4.4 / 3.2) | .885 / .165 | kept poles of 22,173 (purity .09) / 25,564 (.14) / 23,364 (.03) at 1 / 2 / 3 passes |
| 10 % s1, site 1 | z1hi | .284 / .236 / .222 / .221 (4.2 / 3.5) | .898 / .139 | kept pole 2,978 at purity .74 -> 28,325 at .19 / 16,345 at .30 / 2,346 at .86 + a 26,097-seed nuisance pole |
| r131 s2 root site 1 | all four | .39-.43 -> .38-.43 | - | the closed pole {z1lo, z3hi} 2,981 -> two poles 5,667 + 19,620 |
| qtoynull root site 1 | all four | .45-.60 -> .44-.46 | - | no pole -> a kept pole of 26,983 seeds (.45 of the null) |
| olenull root site 1 | all four | .74-.89 -> .74-.86 | - | no pole -> a kept pole of 25,676 seeds (.32 of the null) |

No pass lifts any 6.5 % endpoint to the bar (best z 2.2 against 3), while at 8-10 % the pass destroys the pure rare
endpoints (AUC 0.89 -> 0.15: the hard set at N* 25,000-40,000 holds the rare state at purity 0.12 against the node's 0.08,
less contrast than the EM's own class at 0.74-0.78, and the prior at N* / n = 0.4-0.67 pushes the class past half the
node, where it reads as the complement) and hands kept poles to the null roots. The premise of 6.5 does not hold here: N*
is where the deviation map's cross-half SNR peaks (3-10x the class), not the class size, and the EM posterior at
AUC 0.75-0.78 is already at the per-particle limit of the site's band; rescoring against a lower-contrast reference can
only lose ranking. **No-go**; nothing of it is committed. The lever at 6.5 % would have to raise the per-particle
information at the site (band, mask, pose), not the number of scoring passes.

The detection limit is therefore unchanged by this round: the site layer at 6.5 % is the same as before (no endpoint at
3 SE under any form), and at 8 % the (B) pole {z1lo, z3hi} is tested on its own seeds (section 9.3).

### 9.3 Gate replay of the committed (B) form

Arms `head` (c7ac3b0bd, `build_head`) and `bc` ((B) + the floor of section 10, `build_bc`; the floor's actions are
listed in 10.3 and are the only other difference), same seeds and keys; the replays of the 14 saved sites on `build_bc`
(`lever/sites/bc_*`) say where (B) acts:

| fixture | head | bc | (B)'s action |
|---|---|---|---|
| qtoynull / olenull / micnull md0 | 1 leaf each | 1 leaf each | the took-off endpoints (shares 0.91-0.97) are rejected up front instead of by share; the remaining poles (8-15 seeds) rejected as before |
| qtoy2 md1, seed 1 | 10 leaves, 5-way root 22,118 / 11,473 / 10,516 / 9,272 / 6,621 | 5-way root 21,937 / 11,472 / 10,523 / 9,461 / 6,607 (the leaf count is the floor's, 10.3) | root sites 1-2 identical (17,987 / 9,707 seeds); at sites 3-5 the took-off endpoints no longer back the second poles, so site 3 keeps one pole (8,856 equivalent pure members, not 15,096 + 8,856) and site 4's pole seeds the fifth state instead; the partition moves by <= 200 particles per group |
| r131 md1, seed 2 | 4 leaves, root 25,397 / 18,937, corner at r1 (SEPARATES) | the same tree to the particle; corner at r1 | root site 1 identical (2,981 seeds {z1lo, z3hi}, z1hi / z3lo at 0.75 / 0.74 rejected up front); site 3 loses its 9,076-seed pole (it existed through the took-off z1hi at 0.51) but never seeded a state |
| r131 md1, seed 1 | 4 leaves, root 23,449 / 20,885, corner at r2 (SEPARATES) | the same tree to the particle; corner at r2 | root site 1 identical (16,470 + 6,166 seeds) |
| detlim2 roots 6.5 % (seeds 1, 2) | lost: one leaf (new arm) | lost: one leaf (guard merges at L / L_ref 0.33 / 0.42) | nothing took off at the arm site; the 8 % pole {z1lo, z3hi} now has 1,863 / 1,780 seeds at rare purity ~0.75 (was 717 / 545) and is rejected by its own seed test at z 2.23 / 2.31 |
| detlim2 roots 8 % (seeds 1, 2) | found through the fallback, F1 0.457 | found through the fallback, F1 0.457 / 0.459 (guard 0.75, floor P / N 9.8 / 9.9) | as above: no kept pole, the same fallback |
| OLE md1 | 8 leaves, 4-way root 22,700 / 20,864 / 19,639 / 17,396 | the same root and the same leaves at r1, r2, r4 to the particle; r3 moves by 95 particles (10,456 / 9,183 -> 10,551 / 9,088) | the root sites' took-off endpoints (olenull-like quality axis) are rejected up front; the span is the baseline's |

Verdict for (B): the cause is fixed where it was found (the took-off class no longer takes the seeds of the real pole), the
hard gates and the nulls are unchanged to the particle, and the detection limit does not move (16-17, 1.9-2.1 % at 692k):
at 8 % the real pole's own seed deviation is 2.2-2.3 SE, under the pole rule's 3, so the state is still delivered by the
fallback rather than seeded by the site. Cost: head arm 115 core-h, bc arm 112, four roots 32, site replays about 3.

## 10. The effect floor after the locality guard (the rare-state note's P / N)

### 10.1 What changed

`exec_state_split` calls `effect_floor` after `locality_guard`. For every pair of delivered children it builds each child's
half maps as the count-weighted means of its stage-B states' half maps (the particles each state put into the child), scores
the pair exactly as the guard does (`pair_maps`, the body of `pair_locality`: soft sphere of `mskdiam`, 40-10 A band, the
second map scaled onto the first per shell by the cross-half regression, noise-corrected squared difference smoothed 3 px)
and takes P = the reproducible difference power inside the sphere; N = the same smoothed squared power of the node's own
half difference (the count-weighted mean of all stage-B states, even minus odd, over two). The pair with the smallest P / N
merges while it is under `EFF_MIN` = 2, the children are renumbered by size and K* updated; `split_stats.txt` records K*
after the guard and the floor. Nothing else moves: the split test, the guard and its 0.5 are as they were.

### 10.2 The Fortran against the Python port (`detlim2/scripts/split_effect.py`) on the same saved nodes

The df gate trees' split nodes were replayed with `state_split` on `build_floor` (`lever/splits/`, the node's own command and
stage-B maps) and scored by the port (`lever/splits_py_df.tsv`):

| node (children) | P / N Fortran | P / N port | guard L / L_ref (Fortran = port) | floor |
|---|---|---|---|---|
| qtoy2 root, 10 pairs of 5 | 41.7-81.6 | 43.9-88.2 | 0.98-1.10 | kept |
| qtoy2 r1 (12,700 / 9,401) | 5.23 | 4.75 | 0.98 | kept |
| qtoy2 r3 (6,585 / 3,937; inside combo 28) | 0.19 | 0.87 | 1.15 | **merged** (both under 2; P is at the noise level there) |
| qtoy2 r4 (4,955 / 4,334) | 4.63 | 4.71 | 1.25 | kept |
| r131 seed 2 root (25,401 / 18,933) | 7.77 | 7.65 | 0.92 | kept |
| r131 seed 2 r1, the closed corner (13,332 / 12,069) | 3.18 | 3.05 | 1.07 | kept, margin 1.6x |
| r131 seed 1 r2, the closed corner (11,601 / 9,414) | 2.67 | 2.66 | 1.23 | kept, margin 1.3x |
| r131 seed 2 r2 (9,818 / 9,115; the 0.93 closed-arm extra split) | - | 2.09 | 1.14 | at the floor (the bc arm decides) |
| OLE r1 (12,256 / 11,248), r2 (11,643 / 8,733) | 6.86, 5.82 | 6.74, 5.79 | 1.47, 0.98 | kept |
| qtoynull, olenull, micnull roots | guard merges first | 6.0, 9.9, 44.0 | 0.21-0.24 | - |

Agreement is 1-10 % wherever P is well above the noise (the shell assignment of SIMPLE's `fsc` / `spectrum` against the
port's `rint` radius is the residual); where P sits at the noise level (qtoy2 r3: P 3,192 against N 3,677 in the port) the
positive part of the noise dominates and the two differ by 4x, on the same side of the floor. The floor's one action on
the saved df trees is the within-conformer split the draw-free test added at qtoy2 r3 (both leaves 0.87-0.92 pure in combo
28), which is the split that took the qtoy2 NMI from 0.764 to 0.742 in section 8.4.

### 10.3 Gate replay

| fixture | floor's action in the bc arm | nodes saved | states lost |
|---|---|---|---|
| qtoynull / olenull / micnull md0 | none (the guard merges the root pair first, L / L_ref 0.21-0.24) | 0 | none |
| qtoy2 md1, seed 1 | merges r3 (10,523 particles, both children 0.90 pure in combo 28; P / N 0.218) and r5 (6,607, both in combo 101; P / N 0.220): 10 leaves -> 8, NMI 0.773 -> 0.812, best-leaf F1 of combo 28 0.64 -> 0.94 and of combo 101 0.71 -> 0.95, every other combo within 0.01 | 2 | none (both merged pairs lie inside one conformer) |
| r131 md1, seed 2 | none: root 7.8, r1 (the closed corner) 2.6-2.7, r2 (the 0.93 closed-arm substate) 2.08 against the floor of 2 | 0 | none; the corner keeps a 1.3x margin, the r2 substate 1.04x |
| r131 md1, seed 1 | none: corner at r2 2.7 | 0 | none |
| detlim2 roots 8 % (seeds 1, 2) | none: P / N 9.8 / 9.9 (the 6.5 % roots never reach the floor, the guard merges them) | 0 | none |
| OLE md1 | none: root pairs 31-50, depth-1 splits 5.3-7.9 (K* 4 -> 4 at the root, 2 -> 2 at r1-r4) | 0 | none; the 8 leaves are the head arm's |

The kept pair nearest the floor in the whole arm is r131 seed 2's r2 at 2.08, the extra 0.93-closed-arm split the draw-free
test added in section 8.4; the closed corners sit at 2.6-2.7. The floor's two actions are the two within-conformer splits
of qtoy2 that the draw-free test added over variant D, so with the floor the draw-free default no longer costs NMI (0.812
against variant D's 0.764 / 0.796 under the SVD discovery) while keeping its gain at the corner.

## 11. The merge-up on the numeric node record (D): what is done and what is left

Done in this round: `equiv_tree_inputs` reads a node's margin from `split/split_stats.txt` only (weakest kept pair statistic,
null mean, common n -> `split_effect`); a node without numbers or with K* < 2 contributes nothing, and the log-parsing fallback
`node_effect` (the "[cross-half EM groups] equal N" / "beyond-EM band" / "K* =" lines of the production test's step log) is
deleted. Every node written since c7ac3b0bd carries the numbers.

Not done, scoped by reading `run` / `test` / `test_once`: the merge-up's own pair tests at every level still call the
production `test` (`run`, the two `self%test` calls; the `.not. l_equiv` branch is dead since `merge_rule_setup` sets
`l_equiv`). Replacing them by `test_once(cand, size(tst), child, kstar, junk, self%mtests, ugrp=ugrp)` compiles, but `ugrp`
means different things: `test` maps a unit through the board groups and the second `decide`'s groups, so two units share a
group only when `decide` merged them and a unit the EM drained gets 0 (the comment at the `l_equiv` branch relies on that);
`test_once` sets `ugrp(g)` to the mode of the final `child` over the unit's particles, which is almost never 0 and can put
two cannot-linked units in one group through the choose rule. The port is to derive `ugrp` from `decide_once`'s result (EM
group m starts from unit `gcom(m)`; a unit gets the final index of `gof(findloc(plist, m))` when m is in `plist`, else 0,
followed through `choose_members`' renumbering), then delete `test` (with `stash`, `group_of`, `route_members`,
`add_subset_jobs`, `junk_report`), `decide`, `eqn_stats`, `null_upto`, `pair_stat`, `blom_emax`, `rnd_subset`, the
components `nmed` / `nmax` / `nnull`, `shuffle_null%jsc` / `l_df` / `nul`, the constants `MINSEL_JUNK` / `JUNK_CORR`,
and the `.not. l_df` branch of `null_resolve`. The statistic the merge-up tests with then becomes the draw-free one at the
common n, and `mdelta` (from `split_effect` of the same statistic) is matched to it. Gate: OLE's merge-up must keep the small
OLE states and the closed pair apart (memory `equivalence-mergeup-verdict-2026-10-27`). About 1,100 lines leave.

## 12. Pose-change statistics in the node record (E, report only)

From the coordinator's pose-artefact study (`tree_pose_artefact_test_2026_11.md`, `scripts/hetero_tree_pose_study/pose_split.py`).
`state_split` now reads the node's input poses (`../n.simple`, the poses the particles arrived with) beside the stage-B poses it
refines from, and for every particle takes the rotation between the two minimised over the point group (`euler2m`, the
geodesic angle of R_in^T R_B S_k, `sym%apply` on the right). It appends to `split_stats.txt`, gating nothing:

- `# pose node flip30 f nact n nsym N mean_angle a`: the share of the node's particles that turned more than 30 degrees;
- `# pose split pair g|h n na nb flips fa fb dflip d z z null mu sd flag yes|no` for every pair of delivered children (after
  the guard and the floor; pairs with a child under 50 particles are skipped), `# pose states pair ...` over the stage-B
  states where no split is delivered: dflip = |fa - fb|, its z against 200 shuffles of the pair's labels inside the split
  test's own strata (`strata`: 40 view bins x defocus tertiles of the stage-B project), the flag at dflip >= 0.04 and z >= 5.

The port differs from the Python in one convention it states: the strata are the split test's (stage-B view bins, defocus
tertiles per call), not bins of the input poses; the angle, the threshold, the shuffle count and the flag rule are the same.

Replayed on the saved df gate nodes (`lever/splits_de/`, `build_de`) against the study's numbers:

| node | flip30 (study) | pairs: dflip, z (study) | flag |
|---|---|---|---|
| qtoy2 root (5 children) | 0.0009 (0.0009-0.0020 for qtoy2 / r131 / nulls) | 0.0001-0.0014, z -0.9 to 5.0 (<= 0.0011, z <= 4.8 for needed splits) | no |
| qtoy2 r1, r4 | 0.00005, 0.0002 | 0.0001 / -1.1, 0.0005 / 1.2 | no |
| r131 seed 2 root, r1; seed 1 r2 | 0.0017, 0.0011, 0.0019 | 0.0014 / 4.1, 0.0004 / 0.0, 0.0008 / 0.5 | no |
| OLE root (4 children), r1, r2 (c2) | 0.021 (0.021), 0.0011, 0.0035 | root 0.008-0.021, z 0.9-14.7 (<= 0.022, up to 14.5); r1 0.0003 / 0.0; r2 0.0026 / 3.5 | no |
| qtoynull, olenull roots (guard-merged; the stage-B pair) | 0.0020, 0.00006 | 0.0005 / 0.6, 0.0000 / -0.8 (<= 0.0025, z <= 3) | no |
| micnull root (guard-merged; the stage-B pair) | 0.273 (0.27-0.37) | 0.0706 / 20.6 (0.071 / 22) | **yes** |

The port reproduces the study on every node it was read against; on the standard gates the flag fires only on the micelle
null's refused root pair (the pose artefact the study was built on), and the needed splits with truth stay at dflip
<= 0.0014. The user decides whether the flag becomes a gate; nothing here acts on it.
