# hetero_tree: algorithmic audit (2026-11-01)

Scope: the code on branch `hetero-tree` (published 21fac34db) plus the uncommitted distribution work in worktree
`funny-pike-0046ee` (`claude/tree-dist`), which is the version audited here. Defaults are the committed parameter defaults
plus what `exec_hetero_tree` sets (`kfloor 5000, minpop 15000, maxdepth 4, site_disc yes, site_dirs 1,3, node_proj yes`).

## 1. The algorithm as built

One node `u` holds particles `I_u` at its parent's stage-B poses (the root at the input poses). Every step is a queue job.

1. **mean**: `reconstruct3D`, one state.
2. **flex**: `flex_pca` on `I_u`, `K_u = clamp(|I_u| / kfloor, 2, 8)` preimages, lp 15, state maps reconstructed.
3. **cheap children** (below the root only): spherical k-means of the projection directions (~1300 per bin), the latent
   centred per bin, PC1 of the centred latent, 10 quantile bins of PC1 inside every direction bin; `reconstruct3D` of the
   10 bins (halves only). The root uses the flex state maps instead.
4. **node mask** (`falloff_mask`): the maps low-passed to 15 A and least-squares scaled to the weighted mean inside the
   support; `v = sum_j w_j (E_j - E)(O_j - O)` (cross-half, noise cancels in expectation), Gaussian-smoothed 1.5 px;
   core = `v > 0.2 max(v)`, dilated 3 px; weight `exp(-d^2 / 2 6^2)` outside the core; times the molecule sphere.
5. **fflex**: `flex_pca` again with the mask as `pcg_mskfile`. Its states seed the next step.
6. **stA0**: `refine3D prob_neigh`, lp 15, one pass, box 64, nspace 2500, `pcg_focus` on the mask, the fflex states.
7. **recB0**: `reconstruct3D` of the stA0 states (halves).
8. **site discovery** (`discover_sites`, a worker job): the same cross-half covariance over the recB0 state maps; core at
   20 % of the peak; 6-connected components; each component is split into at most two "factors" by a 2-means on the
   per-voxel state-deviation profiles (`|cos| < 0.8`, both shares >= 0.1, 20 iterations); candidates get a dilated falloff
   mask, a reproducibility `r` (masked cross-half correlation of the deviations), a variance share and a molecule fraction;
   non-local components (mask > 30 % of the molecule) are re-cut at a higher threshold; candidates are ranked by share and
   deduplicated by Dice >= 0.9 and by symmetry mates (Dice >= 0.99 after a group rotation).
9. **site ranking** (`site_rank`, one job per site, all sites at once): starts = the 15k-particle lo/hi tails of fflex latent
   directions 1 and 3 (four endpoints). Per endpoint, a two-class EM in the site mask at fixed poses with cross-half
   references (`site_em`), fitted on at most 100k particles and scored on all; the hold-out correlation curve gives `N*`
   and the "equivalent pure members". Then the halves machinery: the same EM is refitted on two halves and four quarters of
   the rows (six more fits per endpoint); the cross-half correlation matrices of the endpoint deviation maps decide which
   endpoints merge into one pole and whether an endpoint is accepted (self-match minus best other partner >= 3 SE, SE from
   the four quarter replicates); the rest (posterior < 0.5 for every pole) is tested as one more endpoint (six fits); the
   **pole rule** refits each pole's seed set on halves and quarters again (six fits) and keeps the pole when its share is
   under 0.5 and both z >= 3. The margin `min(z_row, z_col)` is the site's strength.
10. **seeding**: `strongest_pole=yes` takes the site with the largest kept-pole margin; its kept poles plus the rest are the
    seed states. If no site has a kept pole, the two complete-linkage groups of the fflex states seed. Below
    `kway_maxdepth` the node goes binary (two-pole seeds or the two largest states). `site_states=yes` (route A, off by
    default) adds other sites' minority poles as extra states.
11. **seeds**: `reconstruct3D` of the seed states at box 64.
12. **stage A**: `refine3D`, lp 15, three passes, box 64, nspace 2500, focused on the node mask.
13. **recB**: `reconstruct3D` at box 88.
14. **stage B**: `refine3D`, lp 10, three passes, box 88, nspace 5000, focused.
15. **split test** (`state_split`, see section 4): the stage-B states are grouped into `K*` supported groups and junk
    states; the locality guard then merges a pair whose difference is diffuse.
16. **decide**: groups with >= `minpop` particles and depth < `maxdepth` become children at the stage-B poses; the rest
    become leaves (a multi-state group is one leaf with particle-weighted maps); junk states become junk leaves.

End phase: leaf labels; a complete-linkage tree over the leaf half-maps (`state_tree`, 40-10 A); plain `reconstruct3D` of
every node's leaves; `landscape_sample` for every finished node with two or more accepted poles; the factorial site layer
(`site_layer=yes` by default): every discovered site is registered, climbs the tree level by level (re-ranked at the parent
from its pole labels, accepted while the first-pole Spearman with its home ranking is >= 0.5), copies and duplicates are
folded by Spearman >= 0.5, kept sites get one re-mask round (class maps -> falloff mask -> re-rank -> guard), a sizing run,
pole maps and tables; the merge-up (`state_mergeup` over the leaves along the leaf tree); the final maps (`state_tree` again,
then `tree_wiener` Wiener-weighted maps).

Per node that is 12 queue steps plus one discovery job plus one `site_rank` job per site. A `site_rank` job runs
4 x 7 + 6 + 6 x (kept poles) EM fits, about 45 for a two-pole site. A root with 7 sites runs about 300 EM fits to choose
one seeding site.

## 2. Size

| file | lines | role |
|---|---|---|
| `src/main/strategies/search/simple_hetero_tree.f90` | 7,610 | controller, step machine, masks, discovery, site layer, motion layer |
| `src/main/strategies/search/simple_state_mergeup.f90` | 4,657 | split test, site EM, merge-up, locality guard, junk screen |
| `src/main/strategies/search/simple_landscape_sample.f90` | 1,254 | per-node landscape sampling (end phase) |
| `src/main/commanders/simple/simple_commanders_refine3D.f90` (tree part) | 1,057 | `site_rank` 800, `state_split` 62, `state_tree` 60, `state_mergeup` 43, `landscape_sample` 35, `site_discover` 32, `hetero_tree` 25 |
| `src/main/volume/simple_state_tree.f90` | 709 | complete-linkage leaf tree |
| `src/main/volume/simple_tree_wiener.f90` | 625 | tree-weighted final maps |
| `src/main/ui/simple/simple_ui_heterogeneity.f90` (tree part) | ~250 | key metadata (`hetero_tree` 177, `site_rank` 67) |
| **total in the run path** | **~16,200** | |

Not counted: `simple_landscape_compare.f90` (733, a separate post-hoc program), `flex_pca`, `refine3D`, `reconstruct3D`.

Growth of `simple_hetero_tree.f90`: 1,162 lines on 2026-10-21, 3,637 on 10-23, 4,269 on 10-24, 7,610 on 11-01, over 33
commits. `simple_state_mergeup.f90`: 1,675 on 10-21, 4,657 now. The whole thing is 6.5x its size eleven days ago.

Where the 7,610 lines of the controller go (procedures bucketed by concern, `advance` split at `motion_next`):

| lines | share | concern | on by default |
|---|---|---|---|
| 1,842 | 26 % | motion / gate / floor layer (`motion_stats`, candidates, gates, nulls, inheritance, product cells) | no (`node_motion=no`) |
| 1,511 | 21 % | scheduler, queue, resume, row-local projects, pruning, subtree jobs, key builders | yes |
| 1,191 | 17 % | factorial site layer (`factorial`, `fs_*`, `state_seeds`, pole maps) | yes (`site_layer=yes`) |
| 939 | 13 % | site discovery and the focus masks | yes |
| 821 | 11 % | the node step machine proper (mean .. split, seeding, decide, leaves) | yes |
| 476 | 7 % | pole bookkeeping: parsing `site_rank`'s ten text files back into the controller | yes |
| 217 | 3 % | string and file helpers | |
| 209 | 3 % | end game (plain maps, landscape, merge-up, final maps) | yes |
| 169 | 2 % | header: 70 named constants, 20 step codes, 17 site states | |

The science the user gated (closed corner, nulls, qtoy2, OLE) runs through the 821-line step machine, the 939 lines of
masks and discovery, the 800-line `site_rank` commander, and the split test. Everything else is scheduling, description
output, or switched-off experiments.

Switches still in the code that the repo's own verdicts rejected or never adopted: `node_motion`, `motion_gate`,
`motion_inherit`, `motion_multi`, `motion_floor`, `two_poles`, `pooled_poles`, `node_seed`, `oneflex`, `site_states`,
`subtree_max`, `SIMPLE_TREE_REENTER`, `SIMPLE_TREE_KEEPMAPS`, and in the split test `SIMPLE_SPLIT_RHOMIN`, `_ROUTE`,
`_QSTRATA`, `_ORTHO`, `_KTGAP`, `_HPMASK`, `_EMNULL`, `_EMNULL_SEED`, `_LOCALITY`, `SIMPLE_MERGE_RULE`, `_DISTINCT`,
`_DELTA`, `_TREEDIR`, `SIMPLE_SPLIT_CONTRAST`. Each is a tried variant that stayed as a branch.

## 3. Where the hours go

From the recorded profiles (sacct by WorkDir, memory notes of 10-25 and 10-29):

- 692k RNase P on recipe17d: 5 h 31 min to `STATES_DELIVERED`, 873 core-h, 84 nodes, 67 leaves.
- Root at 692k, per step: flex + cheap 0.23 h, fflex 0.17-0.20 h, stA0 0.10 h, stage B 0.12-0.18 h, split 0.08-0.11 h.
  The seven site ranks: 0.18-0.40 h each, 1.93 h when sequential. Ranking the sites to pick one costs more than every
  refinement of the node together.
- The factorial layer on full692 (V3): 516 core-h, 4.6 h wall, 188 candidates ranked, 88 kept; 145 pooling climbs alone
  were 199 core-h. Its output is a description (labels and posteriors per site); it changes no split.
- Small nodes (10-50k): 0.25-0.5 h for about ten sequential queue steps, each distributed round paying about 24 s of NFS
  latency. That is latency, not compute.
- End phase: merge-up 19-54 min in one job; `final_maps` 18 min serial in the controller.
- The distribution builds (6-11) cut the NOT r2 fixture from 6.5 h to 2.3 h without changing a result. That is
  scheduling: discovery as a job, rank threads sized to the rows, class maps gated on the core budget.

So the wall is set by (a) the number of serial queue steps per node times the depth, (b) the site-rank EM count, and
(c) the factorial climbs. None of these is the refinement.

## 4. The split test (state_mergeup)

`state_split` runs `state_mergeup%test` on a node's stage-B states, then `locality_guard`:

1. The node's particle Fourier planes are loaded once at the crop box (two pixels per 10 A), with a per-particle noise
   sd, one of 40 Fibonacci view bins, the mean defocus and the stage-B `corr`.
2. **Junk screen** of the states: per-particle SSNR beyond the EM band (20 A to Nyquist) relative to the consensus
   below 0.35 sets a state aside (`junk_screen=no` reports only; RNase P runs that way).
3. **Strata** = view bin x defocus tertile.
4. **Shuffle null**: 20 within-strata label permutations per half (50 under the Bonferroni count); for each, `neq`
   particles per group (`neq` = the smallest group in either half), maps, docking onto the other half's consensus, and
   the pair statistic in the coarse band and in the hold-out band beyond the EM band (`kt` = EM limit + 4 shells with a
   focus mask). `thr` = the maximum of the 20 null values; `thr2` = mean + `z_crit` sd with `z_crit` the one-sided
   normal quantile at (1/21) / `ntests`, `ntests = 2 N_root / minpop - 1` (91 at 692k, `z_crit` about 3.7).
5. **Pair statistic** `T(p,q)` = the peak of the 3-shell-smoothed cross-half FSC between the difference maps
   `A_p - A_q` and `B_p - B_q`, each pair at equal N = `min(|p|, |q|)` per half, drawn with a fixed seed.
6. **Board decision**: greedily merge the smallest `T` while `T < thr` in the coarse band; `K*` groups.
7. If `K* >= 2`: a **cross-half EM**, three iterations to 20 A, references = all + FSC-filtered (group - all), docked,
   an E-step with a per-shell scale and noise; the MAP labels replace the stage-B labels; the decision is repeated in the
   hold-out band against `thr2`, extending the null while a statistic sits within 4 SE of the threshold.
8. **Routing** (`choose`): membership is the EM groups unless the stage-B labels win by 3 SE on the per-particle
   hold-out log-likelihood.
9. **Locality guard**: for each live group pair, the smoothed excess residual of one state's map regressed on the
   other's, as the share held by its top 5 % voxels against the consensus's share; a pair with `L / L_ref < 0.5` merges.

The merge-up (`run`) applies the same test along the halving levels of the leaf tree; its equivalence rule re-reads the
nodes' split logs to set the "shown equal" margin and the cannot-link pairs.

Lines by concern in `simple_state_mergeup.f90` (4,657): the statistical test itself about 890; EM and map building about
1,570 (`site_em` alone 403); nulls 177; junk screening about 600; locality guard 268; docking and symmetry 108; pole
merging for `site_rank` 235; I/O, log parsing and helpers about 515; header 131. Of these, about 490 lines are reachable
only through an environment variable or a flag, about 160 compute a quantity that is logged and never used, and about
250 are duplicates (the per-shell scale-and-noise fit is written four times, the softmax three times, `blom_emax` and
`split_zcrit` are the same bisection, three insertion sorts sit beside the imported `hpsort`, `site_em`'s density and
Wiener helpers repeat `band_density` / `band_wiener` / `band_terms`).

What is stacked on a sound statistic (all verified in the source):

- **Maps are built twice in `decide`.** `allm` holds one equal-N map per group per half, docked, and is freed after a
  logged median; `eqn_stats` then rebuilds four maps per pair at `min(|p|, |q|)`: 2K(K-1) reconstructions and dockings
  where 2K already exist.
- **The null and the real statistic are drawn at different N.** The null uses `neq` per group; a real pair uses
  `min(|p|, |q|)`, and after a merge the merged group is re-tested at a larger N against the same `thr`. Every verdict
  also carries the variance of one seeded random subsample.
- **The threshold is the max of 20 x (number of pairs) null draws**, so the per-pair alpha depends on K; the z-mode
  re-encodes "max of 20" as alpha = 1/21 and assumes a max-of-3-shell-means statistic is normal with an sd from 20 draws.
- **The equivalence rule parses printed log text**: `node_effect` reads `T` back from "pair statistics at the start"
  at two decimals and the null mean at four; `shown_equal` then back-estimates the null sd from (max - median) via
  Blom's expected maximum, adds a one-sided 1.645 sd, and converts through `r = 2y / ((1 - y) n)`, a formula written in
  three places. Four approximations on a null whose sd `null_thr` already has.
- **Six z-criteria in one file**: 1.645 (equivalence), 3 (halves), 4 (null extension), 2 (junk gap), 3 (routing),
  `z_crit` (re-test).
- **Routing treats per-particle log-likelihood deltas as independent** (SE = sqrt(N) x sd) although every delta
  shares the same references and docking.
- **The locality guard is asymmetric but mirrored**: `pair_locality(s, t)` regresses `s` on `t`, is computed for
  `s < t` only and copied to `(t, s)`, so the verdict depends on state numbering; its smoothing is in pixels while its
  band is in A, which makes `L / L_ref` box-dependent and is why `LOC_BOX = 88` and the pad-and-rescale in `group_halves`
  exist.
- **Four junk statistics, two applied**: `junk_report` reconstructs `K*` maps per half for a log flag; the amplitude
  mixture in `junk_screen` runs 500 EM iterations with no convergence test and its result is discarded.
- **The hold-out band start is defined three ways**: `self%kt` honours the focus gap; `site_em` hard-codes
  `kb + 2`; `holdout_band` / `half_rel` ignore the gap. So `site_rank`'s halves test and the split test judge different
  bands.
- On the merge-up path `route_members` (two full-band EM passes) and `junk_report` contribute nothing to the result.

## 5. What is necessary

Tier 1, the core that produced the gated result: mean -> flex -> mask -> fflex -> seeds -> stage A -> stage B -> split
test -> decide, plus site discovery and one ranking to choose the seeding site. Keep.

Tier 2, analysis products delivered beside the states. They cost the most and change no label:
- the factorial site layer (climbs, re-mask, sizing, pole maps): 1,191 controller lines plus 516 core-h at 692k;
- `landscape_sample` per node: 1,254 lines and one job per finished node;
- the tree-weighted final maps (`state_tree` + `tree_wiener`, 1,334 lines): the delivered leaf maps are the plain
  maps; the Wiener maps are an extra representation.
These should be opt-in programs run on a finished tree, not steps of it.

Tier 3, switched off and verdict-rejected: the motion/gate/floor layer (1,842 lines, "separates motion types but not
real from spurious splits", floor dropped), route A re-entry, the 14 split-test environment variants, `pooled_poles`,
`two_poles`, `node_seed`, `oneflex`, `subtree_max` (2.1-2.5x the core-h for no wall gain). Delete, or keep on a branch.

Steps that can collapse without a scientific change (each needs one same-seed replay to confirm):
- **stA0 + recB0** exist only to make the state maps that discovery reads. `flex_pca` already reconstructs the fflex
  states; feeding those half maps to discovery removes two queue steps per node. The one-pass stA0 sharpens the states
  by re-assignment, so this is an ablation, not a certainty.
- **cheap children** exist only to make the node mask below the root, because the root uses the flex state maps
  directly (the port spec's "root exception"). One rule for both levels removes a step and 108 lines; the
  view-balancing the cheap children provide can be applied as per-direction-bin weights when the flex states are
  reconstructed.
- **the text-file protocol** between `site_rank` and the controller (ten files per site, re-parsed by `word()` and
  `str2real`, 476 lines, the source of three recorded traps about row order and missing files). `poles.bin` already
  exists as a cache; make it the only interface.

## 6. Where the mathematics can be tightened

The site ranking and the masks are the right place to look, because they are where the knobs stack: `MFRAC 0.2`,
`MDIL 3`, `MFALL 6`, `MSUP 0.1`, `DFACT_COS 0.8`, `DFACT_SHARE 0.1`, `DFACT_ITS 20`, `DMINVOX 200`, `DMINSHARE 0.05`,
`DDICE 0.9`, `DMAXMOL 0.3`, `DMATE 0.99`, `HALVES_KSIG 3`, `FGATE 0.5`, `SITE_NFIT 100000`, `ntail = min(15000, n/6)`,
`site_dirs 1,3`. Each was set by hand on one dataset, and several of them encode the same question answered twice.

### 6.1 The mask core from a null, not from the peak

`v(x) = sum_j w_j (E_j - E)(O_j - O)` is a cross-half covariance: its expectation is the between-state variance, and its
null distribution is known from the same maps. Under no signal `Var v(x) = sum_j w_j^2 sigma_E,j(x)^2 sigma_O,j(x)^2`, and
`sigma^2` per map and voxel is estimated from `(E_j - O_j)^2 / 2` after the same low-pass. A core defined as
`v / sqrt(Var_0 v) > z` (z = 3, say) is node-adaptive and knob-free. The 20 %-of-peak rule depends on the strongest motion
in the node, which is exactly the failure the port spec records (masks covering 72-94 % of the molecule when a global
contrast gradient passed the 20 % core everywhere, and masks at 24-47 % elsewhere). The dilation and the Gaussian falloff
are the user's validated prior and stay, but the falloff length can be set from the correlation length of `v` rather than
6 px.

The same z-map is the locality guard's statistic: a split is "diffuse" when its significant-z voxels are spread over the
molecule, "localized" when they sit in one component. The guard's `L / L_ref < 0.5` and the mask's core are two readings
of one map.

### 6.2 Site discovery as a cross-validated factorization

A component's voxel-by-state deviation matrix `D` (`n_vox x K`, available in both halves as `D_e`, `D_o`) is factorized
today by a 2-means on normalized profiles with three thresholds and 20 iterations, and "one factor or two" is decided by
`|cos| < 0.8` and share >= 0.1. The principled version is the weighted SVD of `D_e`: the number of factors is the number of
left singular vectors whose projection of `D_o` correlates with `D_e`'s above the half-map noise ceiling (the test
`landscape_compare` already implements for pooled components). That gives the factors, their rank, and their
reproducibility in one deterministic pass, with no `DFACT_*` constants, and it generalizes beyond two factors per
component.

### 6.3 Endpoint reproducibility by one subspace test, not thirty refits

`site_rank` refits the two-class EM on halves and quarters to ask whether an endpoint's deviation map is reproducible and
distinct from the others. The answer is already in the first fit: each endpoint's cross-half deviation maps `d_e^{even}`,
`d_e^{odd}` are produced by `site_em`. For the `E` endpoints of a site, the `E` even deviations span a subspace; project the
odd deviations onto it; a direction is reproducible when its odd-half projection exceeds the noise ceiling (from the
half-map FSC inside the mask, with `n_eff` from `mask_neff`, which the code already computes for the hold-out SE). Two
endpoints are the same pole when their deviations are collinear (principal angle below the noise angle). That replaces
`merge_poles_halves`, `halves_cand`, the twin rule, `pooled_rule`, the rest test and the pole rule's refits, cuts the EM
count per site from about 45 to 4, and removes a statistical weakness: the present z is a margin over an SE from four
quarter replicates (three degrees of freedom), the quarters are subsets of the halves so the replicates are not
independent, and the same 3 SE bar is applied five times per endpoint (merge, accept row, accept column, pole row, pole
column) with no multiplicity control.

The hold-out correlation curve (`N*`, the peak, its Fisher SE) is the one principled selection already in `site_rank`.
It should be the acceptance criterion: a pole is real when its hold-out peak exceeds the noise ceiling by its SE, and its
size is `N*`. The halves machinery re-derives a weaker version of this.

### 6.4 "Same variable" in map space, not in posterior space

The factorial layer folds two sites into one variable when the Spearman correlation of their particle posteriors is
>= 0.5 (`FGATE`). A posterior correlation measures how many particles the two classes share, not whether they are the
same structural change. Two coupled motions (the RNase P closed and open arms are conditional on each other) correlate
and get folded; two independent sites with a common quality axis correlate too. The structural question has a
structural test: the angle between the two sites' pole-difference maps on the union of their masks, against the
half-map noise. Both maps exist when the test is made. Mask Dice is computed and written, but never used in a decision.

### 6.5 Climb by scoring, not by refitting

A pooling climb re-runs the whole `site_rank` at the parent (an EM with all the halves refits) to ask whether a child's
site is the same variable one level up. If the site's pole maps are fixed, the question is one E-step: score every
particle of the parent against the two pole maps at its pose (the per-particle log-likelihood ratio that
`landscape_sample` already computes as its coordinate), then test reproducibility with the subspace test of 6.3. That
turns 199 core-h of climbs into a few core-h, and it is the same quantity route A needs for its seeds. Whether climbing
is needed at all is a separate question: the memory record says the pooled sites drove no split.

### 6.6 One site object

Today a "site" is re-derived four ways: `cand_score` (mean of endpoint posteriors), `fisher_score` (a Fisher discriminant
in the flex latent), `landscape_sample`'s `t = s_a - s_b` (a pole log-likelihood ratio), and `fs_sets`' top-EPM labels.
One definition serves the tree, the site layer, the landscape and route A: a mask `M_s`, two pole deviation maps per
half, and a per-particle log-likelihood ratio at fixed pose. With that object in a binary record the 476 lines of pole
bookkeeping and most of the `fs_*` state machine disappear.

### 6.7 The split test: one statistic, one null, one band

The pair statistic (cross-half FSC of the difference maps, equal N, docked) is the right object and it is what passed
the gates; the bloat is around it. One pass does the same job: the cross-half EM once; `allm` built once at `neq` for
every group and half; `T(p, q)` for every pair read off `allm` from one FSC curve that serves both bands; the null built
from the same `allm` construction at the same `neq`; one `z_crit`. That removes `decide`'s second map build, the
null-versus-real N mismatch, the board-then-EM double decision, and the log-parsing merge rule, because the node then
writes `T`, the null mean and sd as numbers in its record and the merge-up reads them.

The junk screen stays (the user's measurement: about 4 s per rank after 9005a8c77, so cost is not the question). What
changes is its use, from a hard gate to two soft ones. First, the zero-contrast posterior `p_empty` from `fit_contrast`
as a per-particle weight `1 - p_empty` in the node's reconstructions and the site EM, instead of removal: the hard
amplitude rule took 9,016 genuine CAHRA particles, and RNase P runs with the screen off for that reason, while a weight
costs nothing when the posterior is wrong by a little. Second, a per-node quality-split flag: the children's mean
contrast and `p_empty` against the within-child spread, written to the node record and never gating. That is the
acquisition flag the motion layer wanted and the cheapest detector of the quality-driven partitions on record (m33's
root site 6). The SSNR unit screen in `test` is unchanged; `junk_report`'s reconstructions are the one part that is
only a log line.

### 6.8 What this does not fix

Two of the standing verdicts are warnings for anyone making the test "more principled": the held-out split test (v1-v5,
19 nodes) was not adoptable because docking, EM refits and best-group noise fits broke it, and true e/o halves with
per-half seeding passed every null but lost the r131 closed corner. A cleaner statistic is only adopted on the same
gates: nulls one leaf, qtoy2 not worse, the closed corner retained, OLE span retained. Every item above is a replay on the
existing fixtures first (maxdepth 0-1, about 40 min), never a new default.

## 7. Suggested order

1. Delete the switched-off layers and env variants (about 2,500 lines, no scientific change, one fast-gate build).
2. Move the factorial layer, landscape and tree-Wiener maps to opt-in post-tree programs (the 692k tree then stops at
   the leaves, labels and merge-up).
3. Replace the halves refits in `site_rank` by the subspace test (6.3); same-seed replay on qtoy2, r131 and the nulls.
4. The single-pass split test (6.7): same replays, plus the merge-up on the 692k leaves.
5. Null-calibrated mask core (6.1) and the SVD discovery (6.2); same replays.
6. Collapse stA0/recB0 and cheap children (section 5); same replays.
7. The site object and binary interface (6.6), which is the refactor that makes the controller small again.

## 8. Step 1 applied (2026-11-01)

Worktree `funny-pike-0046ee` (branch `claude/tree-dist`, HEAD 9005a8c77), uncommitted, not compiled.

Removed:
- the motion / gate / floor layer: the keys `node_motion`, `motion_gate`, `motion_inherit`, `motion_multi`, `motion_floor`
  (params, parser, UI); the steps `S_MOTION`, `S_MBATCH`, `S_GATEW`, `S_OFFP`; `motion_next` .. `gate_rec` inside `advance`;
  `motion_stats`, `motion_factor`, `product_cells`, `cand_*`, `fisher_score`, `lin_solve`, `quint_bins`, `gate_corr`,
  `floor_*`, `anc_*`, `off_*`, `perm_groups`, `groups4`, `post_groups`, `lcg_perm`, `gate_seed`, `write_coords`,
  `cand_mask`, `root_axes`, `pearson`, `quant_of`, `node_motion_summary`, and the orphaned `write_line`, `kept_ends`;
- `SIMPLE_TREE_REENTER` (the re-entry branch of `decide`, `rest_sel`, `reentered_line`);
- in `state_mergeup`: `SIMPLE_SPLIT_QSTRATA` (and the per-particle `corr` column `qv`), `SIMPLE_SPLIT_EMNULL` / `_SEED`
  (`em_null`), `SIMPLE_SPLIT_KTGAP` (`RB_MARGIN_FOCUS` = 4 fixed), `SIMPLE_SPLIT_RHOMIN` (`rel_effect`, the floor),
  `SIMPLE_SPLIT_ORTHO` (`odiff`), `SIMPLE_SPLIT_HPMASK` (the high-pass block of `finish`), `SIMPLE_SPLIT_ROUTE` (membership
  is always the choose rule; the `best` / `band` / `band_direct` / `em` branches are gone), `SIMPLE_MERGE_RULE=off` (the
  equivalence rule is always on), the `SIMPLE_MERGE_DELTA` / `_DISTINCT` overrides, `SIMPLE_SPLIT_LOCALITY` (the guard is
  always on, stage-B directory `../stageB`).

Kept on purpose: `SIMPLE_SPLIT_CONTRAST` (recipe17d exports `maps`; removing it would change the production split test),
`SIMPLE_TREE_KEEPMAPS` (replays), `SIMPLE_MERGE_TREEDIR` (a path override), and the parameter switches `two_poles`,
`pooled_poles`, `node_seed`, `oneflex`, `subtree_max`, `site_states` (not in step 1's scope).

Counts: `simple_hetero_tree.f90` 7,610 -> 5,725; `simple_state_mergeup.f90` 4,657 -> 4,395; 2,260 lines removed and 78
added over five files. Checks: no dangling reference to any removed symbol, block structure intact (a structural checker,
not a compiler), no procedure left without a caller. The build and the fast gate are the user's.

Consequences to know: launchers that export `SIMPLE_SPLIT_RHOMIN`, `SIMPLE_SPLIT_ROUTE=em` or `SIMPLE_MERGE_RULE=off` now
silently run the defaults (floor off, choose, equivalence rule). The step codes kept their numbers (`S_SBATCH` 18,
`S_SDISC` 20), so a node table from the previous build resumes; a tree stopped inside a motion step (15-17, 19) cannot
and has to be restarted at its site step.
