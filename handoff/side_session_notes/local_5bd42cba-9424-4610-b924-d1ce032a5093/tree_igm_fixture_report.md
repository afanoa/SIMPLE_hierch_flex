# hetero_tree IgM-like fixture: prototyping 5 input-free changes (session started 2026-11-05)

Work dir (batch): `/mnt/beegfs/elmlund/afan/igm_cryolo/tree/fixture_igm/` (scripts/, logs/, fix/, null/, trees/).

## Fixture design (scripts/igm_fixture.py, make_fixture.sbatch)
- Core = the REAL IgM consensus (relion/full/Cons/sum_w.mrc, the tb100k pose frame) Fourier-cropped to box 200 / 3.27 A, cut to its envelope inside r 140 A. Keeps the real view distribution of the core.
- 5 synthetic Fab-like arms (pseudo-atoms in two ellipsoid lobes, 50 kDa each on the map's own mass scale: whole IgM = 970 kDa in r < 300 A), pivots at r 145 A in the disc plane (normal = the core's smallest inertia axis), tips ~265 A.
- fix: discrete hinges, independent: a1 out-of-plane -25/+25 (.5/.5), a2 in-plane -30/0/+30 (.45/.35/.20), a3 out -20/+20 (.65/.35), a4 in -20/+20 (.70/.30), a5 out -30/0/+30 (.40/.40/.20) -> 72 combinations.
- null twin: each arm uniform over the same range (2.5 deg steps). Same rows, poses, CTF, quality and noise draws.
- Images: real tb100k poses + CTF (100k), simulate_particles per component (linear), OLE-null quality model, colored noise calibrated per band to the REAL tb100k particles' band correlation with the consensus projection (micelle recipe). Input poses = truth + uniform +-1.5 deg per Euler / +-0.5 px (core-locked frame, no refinement).
- Arm scale: the map's whole-IgM mass estimate was inconsistent (core integral > the r<300 A total), so the core is DEFINED as 1 MDa and each arm = 5 % of its integrated density (the task's 50 kDa on ~1 MDa; real IgM Fab/Fc5J is nearer 10 %, so this fixture is harder than real).

## Progress
- 2026-11-05 13:18 maps-only check (54957625, 3 min): geometry OK, pivots r 143-147 A, tips 258-271 A, 72 combos, counts per arm as designed.
- 13:22 full build 54957741: clean_snr / clean_syn / stk_core (100k, 40 min) done; the arm stacks failed (xargs -L 1 joined lines with trailing blanks). Fixed (blank strip, temp-name writes). The resubmit watcher died with the VPN (ssh 255) -> ~1.5 h lost.
- 15:39 resubmitted as 54960473 (arm + null components, measure, calib, merge, import).
- Parent fact (2026-11-05, CORRECTED): the 3-state de novo on tb100k has NO 6-fold class (s1 43.7k no ring symmetry = side views/junk, s2 29.9k 5-fold with the J-chain Fab rod, s3 26.4k 5-fold); real tb100k = pentamers + side views/junk, not a pentamer/hexamer mix. Do not read a real-data root split as pentamer vs hexamer. The fixture does not depend on it.
- 16:05 real tb100k root, current tree, maps kept: trees/real_root_stock (54961034), so the (a)/(b) prototypes are also read on real IgM.

## Prototype tooling (scripts/, all offline or via existing hooks; no Fortran change)
- proto_regions.py (numpy mirror of half_variance / falloff_mask / discover_sites, checked against the tree's own mask_variance.mrc + mask.mrc on each run):
  - (a) `supp`: the 0.2-of-peak core per connected region of the molecule support (literal); `island`: the stock core, plus every vs peak whose own 0.2-of-peak core touches no kept region (>= DMINVOX voxels, inside the support). No new constants (MFRAC, DMINVOX, MSUP). On diffuse variance no island is added, so it reduces to the stock mask.
  - (b) every region (stock pieces + islands) is a site candidate over the whole support, kept when the discover_sites SVD test passes (even-half factor reproduced by the odd half, z >= 3, share >= 0.2); site masks in the tree's format for a replay.
  - Unit test on synthetic state maps from the fixture volumes: one site per arm (coverage .90-1.00, no cross-arm overlap), 3-angle arms can give 2 sites (one per hinge direction, as discover_sites' factor split).
- fx_replay_sites.sbatch: copies a finished root, puts prototype sites in node_r/sites, sets the root to S_SDISC (20) and resumes (root or full keys) -> the tree's own site ranking, seeds, stage B, split test, site layer + factorial labels on those sites.
- proto_pn.py (e): the effect floor's P / N recomputed (pair_maps mirror) over the mskdiam sphere and inside each site mask.
- fx_subtract.sbatch + subtract.py (c): rest-of-molecule subtraction (simulate_particles projections of consensus x (1 - site mask), per-particle band-limited LSQ scale) -> new project -> focused root (MASK_FILE = site mask; its stA0 re-aligns locally inside the site).
- score_fx.py: NMI per arm / per combination, group purities, null largest-group share, eta^2 of quality draws.

## Findings so far
- TRAP (fixed in the tooling): runs given `infile=` rows are subtree runs (l_subtree), so the ROOT focus mask comes from the 10 CHEAP-children maps, not the 8 flex maps (tb_tree, igm_tree and the 683k run all pass infile). The numpy mirror on the cheap maps matches the tree exactly (corr(vs) 1.0000, mask Dice 1.0000); on the flex maps it does not (0.85 / 0.82).
- Real tb100k root (current tree): the cheap-children cross-half variance is diffuse, not core-peaked: the 0.2-of-peak core is 390,693 voxels (a ~145 A ball, 14x the 26,929-voxel molecule support), so the stock mask = everything; no island is added (one connected variance blob). The core-only mask of the 683k run is a depth-1 (r1) effect; the root problem on real data is a mask that covers everything and a variance led by one diffuse factor (SVD share .87, cross-half r .46).
- 17:03 FIXTURE DONE (54960473, 1 h 24 min, 110 GB peak): fix/src.simple and null/src.simple (100k each, box 200 / 3.27 A). Measured on 2,000 real tb100k particles: qsig .373 (amplitude spread), band noise powers .089/.044/.026/.014. Realism (band corr with the consensus projection, 120-30 / 30-20 / 20-12 / 12-7 A): fix .163/.052/.028/.020, null .156/.049/.026/.019, REAL .144/.050/.028/.021 -> matched within .02 (lowest band) and .002 (others).
- 17:04 current tree (build_igm4 stock, maps kept): fix_root_stock 54962026, fix_full_stock 54962027, null_root_stock 54962028, null_full_stock 54962030. Chains: early region prototypes when each root's mask exists, then a fresh root with the island mask (fix_root_a / null_root_a).
- ROOT MASK, fixture = real (current tree): the 10 cheap children are ordered by a BACKGROUND OFFSET: their solvent means run monotonically +0.23 .. -0.21 of the support's mean density (fixture), -0.20 .. +0.16 (real). The cheap latent's first axis tracks the particles' background level, so the root's cross-half variance is a reproducible offset^2 over the whole 300 A sphere (fixture: 0.2-of-peak core 1.0 M voxels vs 35k molecule; null identical, 1.1 M). Stock root mask = 100 % of the support, core + arms .87-1.00 (null .92-.99); no island is added.
  - With each half map's solvent mean removed first (`--bg yes`, no constant): the variance localizes; the stock rule then gives the CORE-ONLY mask (core .96, arms .00-.06), i.e. the 683k r1 failure. (a) literal per-support-component cores: arms .57/.50/.74/.73/.48 (the arms are separate components of the 0.1 support at this box); islands: arms .00/.51/.74/.43/.49 (5 islands at r 167-251 A, peaks .04-.05 x max).
  - Real tb100k with the offset removed: core 33.6k voxels, mask 38-96 % of the molecule, no island (the smeared real arms are not separate regions).
- 17:31 fix_root_a (54962574): fresh root with MASK_FILE = bg + per-support-component mask from the stock root's own cheap maps. The first, identical-mask fix_root_a (54962467) was cancelled with its 6 workers (by ID) and moved to trees/_aside/.
- Real tb100k root, current tree (real_root_stock, 1 h 58): LEAF. Sites all core (tree's sites.txt: 3 sites, centres 95 / 48 / 82 A). The split test PASSED (K* 2, 52.5k / 47.5k, z 79.7) and the LOCALITY GUARD merged it (L / L_ref 0.479 < 0.50, "diffuse"), so on this real root the guard, not the effect floor, is the stopper. Pose flip at stA0 9.8 % (> 30 deg), not flagged.
  - (e) on that real stage-B pair: P/N 1.40 over the sphere (pair_maps mirror) vs 6.0-11.8 inside the tree's own core sites and the bg-region masks. Needs the null's numbers before it means anything (max(res^2 - nse^2, 0) is biased upward in small masks).
  - recB0 variance (sites input) on real data: core only (regions at r 16 / 96 A), no island with or without the offset removal; recB0 offsets are small (|.02|), the background ordering is a cheap-children effect.
- FIXTURE REPRODUCES "SITES ONLY IN THE CORE" (current tree): fix root 1 site (centre 33 A, 14.7 % of the molecule); null 5 sites (32-68 A). The recB0 states differ reproducibly in the rigid core (SVD factor r .99), a nuisance (quality / alignment), and the 0.2 x max cut + DMINSHARE keep only that.
- (b) prototype on the SAME recB0 maps (proto/<tag>_recB0_bg{no,yes}): stock rule core .79-.93, arms .00; the island rule adds 5-7 islands at r 150-245 A, each with a reproducible SVD factor (z 3.4-5.6, share .33-.52):
  - fix, offsets removed: 6 sites = the core + ONE PER ARM (a1 .63, a2 .50, a3 .59, a4 .59, a5 .50 of each arm envelope; no site touches two arms). Without offset removal a1 is missed (5 sites).
  - null: also the core + one site per arm (its arms move too); discrete vs continuous is for the per-site layer to tell ((d) + its gates).
- 18:20 (c) started on two of those sites (a1 = recB0_bgyes site 06, a2 = site 02): fix_m_a1 / fix_m_a2 (54963618 / 54963619) = focused root with MASK_FILE = the site on the full particles (no subtraction, control); fix_c_a1 / fix_c_a2 = rest-subtracted particles (consensus = the stock root's mean map) + the same focused root.

## Task 2: the current tree on fixture and null (build_igm4 stock, maps kept)
| run | wall | root mask | sites (tree) | split test | guard | result | vs truth |
|---|---|---|---|---|---|---|---|
| fix_root_stock | 1 h 44 | 100 % of support (offset-driven) | 1, core (33 A) | PASS z 56.8, 68,288 / 31,712 | merged, L/Lref 0.275 (diffuse) | 1 leaf | pre-guard groups: NMI 0.000 with every arm, quality eta2 <= .001 |
| null_root_stock | 1 h 38 | 100 % | 5, core (32-68 A) | PASS z 90.2, 54,801 / 45,199 | merged, 0.428 | 1 leaf (null gate PASS) | angle / quality eta2 0.000 |
| fix_full_stock | 1 h 44 | 100 % | 1, core | same as root | merged | 1 leaf | site layer: core site only (fsite_r_s1) |
| real_root_stock (tb100k) | 1 h 58 | 100 % | 3, core (48-95 A) | PASS z 79.7, 52.5k / 47.5k | merged, 0.479 | 1 leaf | - |
- Failure modes reproduced: sites only in the core; no arm split (the pre-guard root split is a pure nuisance partition, correctly merged by the locality guard). Differences from the description: the ROOT mask is not core-only but everything (background-offset ordering of the cheap children; core-only appears once the offset is removed, as at the 683k r1), and the stopper here is the locality guard, not the effect floor (the floor never runs after K* 2 -> 1).
- (e) caution from the null: its stage-B nuisance pair has P/N 1.89 over the sphere but 2.6-3.6 inside each ARM site and 13.8 inside the core site -> a site-mask P/N cannot reuse the sphere's 2.0 floor; it needs its own null calibration.

## (a) focus mask per region (fresh roots with MASK_FILE = offsets removed + per-support-component cores, from each stock root's own cheap maps)
| run | wall | mask | sites (tree) | split test | guard | floor | result | vs truth |
|---|---|---|---|---|---|---|---|---|
| fix_root_a | 1 h 43 | 94 % of support, arms .48-.74 | core only | PASS K* 3 (42.7k/41.3k/16.1k) | kept (L/Lref 1.41-1.51) | P/N 2.20, 1.17, 1.26 -> 1.63 -> merged | 1 leaf | pre-floor groups NMI <= .002 per arm |
| null_root_a | 1 h 38 | 94 % | core only (3) | PASS K* 2 (50.7k/49.3k) | kept (1.471) | P/N 1.67 -> merged | 1 leaf (gate PASS) | angle eta2 <= .001 |
- (a) alone does not reach the arms: discovery still cuts at 0.2 x the node max (core only), so the seeds and stage B stay on core nuisance.
- (a) shifts the null's protection from the locality guard (stock: diffuse, 0.43) to the effect floor (1.67 < 2): a localized focus mask makes any difference look localized. Combined with (e) (a looser in-site floor) this would expose nulls.
- REPLAY NOTE (fix_full_b / null_full_b, 54964439 / 54964440): resuming at S_SDISC without sites/hetero_step.exit re-runs site_discover, which rewrote sites.txt with the tree's own core site(s); the root's seeds therefore came from the core (stock). The site layer had registered all injected sites from the original sites.txt first and ranks the ARM sites from the untouched site_0k.mrc files (md5 checked), so these replays test (d) on the (b) arm sites, with a stock root decision. fx_replay_sites.sbatch now writes the exit file.

## (d) baseline: routes/toy55/site_product_readout.py (gate >= minpop hard, decisiveness >= .5; dependence; product) on the CURRENT full trees
- fix_full_stock: 1 root site (core r_s1), gate MOTION (decisiveness 3.05), product 50,121 / 49,879 -> NMI 0.000 with every arm (a fake factor).
- null_full_stock: 2 core sites, both MOTION (1.91, 2.52), dependent (NMI .56, r .96) -> 1 motion, product 50,607 / 49,393, angle eta2 0.000 -> a product step with this gate SPLITS THE NULL (largest share .51 < .80). The 10-29 decisiveness gate does not reject this fixture's core-nuisance poles.

## (e) P / N inside site masks, every stage-B pair available (fx_pn.sbatch; clipped as pair_maps / unclipped are identical: the nuisance pairs differ everywhere there is density, so positive-part bias is not the issue)
| pair (nature) | sphere | core site | arm sites |
|---|---|---|---|
| fix_root_stock 37.5k/62.5k (nuisance, NMI 0) | 1.57 | 6.15 | 2.32-3.15 |
| null_root_stock (nuisance) | 1.89 | 13.8 | 2.57-3.57 |
| fix_root_a, 3 pairs (nuisance, NMI <= .002) | 1.17-2.18 | 29-68 | 10.1-16.9 |
| null_root_a (nuisance) | 1.67 | 43.7 | 6.8-9.5 |
| real_root_stock (core sites only) | 1.40 | 6.1-11.8 (3 core sites) | - |
- (e) as specified is a NO-GO: the node's half-map noise is spread over the whole sphere while any density-coupled nuisance difference sits on the molecule, so ANY molecule-region mask raises P/N above the 2.0 floor for pure-nuisance pairs (nulls included). What does separate here is RELATIVE: on nuisance pairs the arm sites read 0.2-0.5 x the core site; a site-local floor would need the site's P/N relative to the rest of the molecule (a per-site locality ratio, like the guard's L / L_ref), calibrated on the null. To be checked on real arm pairs (fix_m / fix_c runs).
- (d) gate candidates on the CURRENT trees' fake (core-nuisance) site factors (proto_sitegate.py: the locality guard's L / L_ref and P / N per pole pair, from factorial/maps/halves):
  - fix r_s1: L/Lref 0.419 (diffuse -> out), P/N 0.24. null r_s4: 0.410 (out), P/N 0.26. null r_s2 (4 poles): pairs 0.48-0.55 (borderline: in by max, out by min), P/N 0.24-0.26.
  - every fake pole pair is far under the effect floor (P/N ~0.25 << 2); the arm poles of the replays decide which measure separates.

## Go local: per-arm focused roots (MASK_FILE = one (b) arm site, full particles; current tree otherwise)
| run | wall | sites (tree) | split test | guard | floor (sphere) | vs truth (pre-floor groups) |
|---|---|---|---|---|---|---|
| fix_m_a1 | 1 h 32 | 2, on the arm (r 168 / 184 A) | PASS z 22.7, 52,127 / 47,873 | kept 1.41 | P/N 0.241 -> merged | a1 NMI .128 (each group ~71 % pure in a1), every other arm .000, quality eta2 0 |
| fix_m_a2 | ~1 h 35 | 2, on the arm (165 / 172 A) | PASS z 6.1, 69,869 / 30,131 | kept 1.42 | P/N 0.191 -> merged | a2 NMI .014 (weak; 3-angle in-plane arm) |
- The IgM failure in miniature: a REAL, arm-specific, localized split (a1) is killed by the sphere-wide effect floor.
- (e) on these real arm pairs: P/N 0.40 / 0.36 over the sphere, 0.36 / 0.28 inside the core site, 55.4 / 56.6 INSIDE ITS OWN ARM SITE (x150 the core). Nuisance pairs: arm sites 0.2-0.5 x the core site. -> a RELATIVE site floor (P/N in the site over P/N in the rest of the molecule) separates real arm differences from nuisance by two orders of magnitude, where the absolute in-site floor cannot (core nuisance 6-68).
- Open: the null's continuous arms (null_m_a1/a2, null_c_a1/a2 running): a median cut of a continuous swing is also a real localized difference, so discreteness must come from somewhere else.

## Discrete vs continuous at a site (needed because the null's arms move too)
- mix_test.py / fx_bins.sbatch: a focused root's particles binned in 4 quantiles of their stage-B evidence difference, bin half maps reconstructed, each middle bin fitted inside the site as an affine (optionally + scale + first-order B) combination of the end bins; reproducible residual (cross-half z >= 3) = continuum.
- fix_m_a1 (DISCRETE arm): bins hold a1-state-2 fractions .17/.41/.58/.82 (monotone mixtures; a2 flat). Bin 2 residual z 0.4 (mixture), bin 3 z 4.2 (5.2 with the B terms) -> the naive verdict calls the discrete arm "continuum-like". The arm difference is very reproducible (end-bin cross-half r .986), so any systematic bin effect (quality / pose drift of the focused re-alignment) shows; the null's number decides whether the magnitude separates.

## (b) + (d) replay through the EXISTING site layer (fix_full_b)
- Root: stock decision (core seeds), 1 leaf. Site layer on the 5 injected arm sites: "no accepted pole, not a site" for ALL five. Their halves test: endpoint self z 0.4-1.8 (needs 3), endpoints of 500-40,000 particles: the site ranking's latent (node-level) does not carry the arm states at this SNR.
- Contrast: a focused root on the SAME a1 mask finds the a1 states (NMI .128). -> for (d) the per-site classification must come from a site-focused analysis (its own cheap children / fflex / stA0 / stage B inside the site), not from the existing site ranking. Launched fix_m_a3/a4/a5 (54966307-09) to build the full per-site product.
- (e') RELATIVE site floor = P/N inside the site / P/N in the rest of the molecule support (proto_pn.py):
  - real arm pairs: fix_m_a1 in its a1 site 108, fix_m_a2 in its a2 site 153 (in the core site: 0.03 / 0.04);
  - nuisance pairs (fix/null stock and (a) roots): arm sites 0.09-0.42, core sites 1.20-2.12.
  - ~50x margin, but the existing constant 2.0 sits AT the null's core value (2.12): adopting it needs a new, null-calibrated threshold (design decision). Pending: the null's continuous-arm pairs.

## (c) local re-alignment of the site with the rest subtracted (fix_c_a1: subtract.py + focused root on the subtracted particles)
- Subtraction itself is clean (per-particle LSQ scale positive for 100 %, P5-P95 .00085-.0022, consistent with the fixture's amplitude/B spread).
- The focused root's stA0 re-alignment on the 50 kDa arm alone SCATTERS the poses: flip30 0.754, mean 83 deg (unsubtracted focused root and stock root: 0.057-0.058, 9.8 deg).
- Its split (PASS z 12.8, guard 2.81, floor P/N 0.48 -> merged) has NMI 0.000 with a1 (unsubtracted focused root: .128). -> (c) as specified is a NO-GO at 50 kDa: the arm cannot anchor its own alignment; the core must. The variant worth testing is subtraction at the core-locked poses without re-alignment (fixposes path), not available in build_igm4.
- fix_c_a1 wall 1 h 32 + subtraction 25 min.

## NULL per-arm controls (focused roots on the null twin's a1 / a2 sites, current tree otherwise)
| run | split test | guard | floor | pre-floor groups vs truth |
|---|---|---|---|---|
| null_m_a1 | PASS z 14.3, 65,065 / 34,935 | kept 1.40 | P/N 0.269 -> merged | a1 angle eta2 0.001; view / defocus / quality eta2 0.000 |
| null_m_a2 | PASS z 12.0, 59,051 / 40,949 | kept 1.42 | P/N 0.205 -> merged | a2 angle eta2 0.003; view / defocus / quality 0.000 |
- Every site-level statistic tried reads the null's nuisance split like the fixture's real arm split:
  | statistic | fix_m_a1 (real, a1 NMI .128 / eta2 .172) | null_m_a1 (tracks nothing) |
  |---|---|---|
  | split test z | 22.7 | 14.3 |
  | guard L / Lref | 1.41 | 1.40 |
  | P/N sphere | 0.24 | 0.27 |
  | relative in-site P/N | 108 | 73 |
  | structural share of the state difference / cross-half z | 0.996 / 52.7 | 0.98 / 46.4 |
  (a2 pair: z 6.1 vs 12.0, reversed.) The pose-scattered subtracted pairs (fix_c, NMI 0) read the same too (relative 43-47, structural z 44-46).
- Mechanism: the merged-reference noise-lock of stage B at good poses (memory tree-good-pose-null-fails-halves-and-seeding-2026-10-30): a partition chosen with references built from both halves makes its own group maps reproducible across halves, so post-hoc half-map statistics certify it. In a whole-molecule root the guard (diffuse) or the floor stops it; inside a focused site the guard cannot, and only the sphere-wide floor (which also kills the real arm split) keeps the null.
- Binned-map mixture test on the null a1 site: evidence bins do not sort by angle (bin means -0.25 .. +0.23 deg, sd ~14), end-bin difference cross-half r .982, middle bins = clean mixtures (z 0.8 / 1.2, "discrete-like"); the fixture's real a1 gave z 4.2-5.2. The test cannot separate either: a noise-locked binary partition looks like a two-state mixture.

## SCORECARD (fixture = 5 discrete arms, 72 combos; null = the same arms continuous; 100k, box 200; cost = wall / allocated core-h per run)
| change | run(s) | fixture | null twin | cost | verdict |
|---|---|---|---|---|---|
| current tree | fix/null root + full | 1 leaf; pre-guard root split NMI 0.000 (nuisance); sites core only; site layer: 1 fake core factor (NMI 0) | 1 leaf (gate PASS); site layer: fake core factors (a product would split it 51/49) | root 1 h 44 / 42 core-h; full 1 h 50 / 53 | baseline: IgM failure reproduced |
| (a) mask per region (needs the background-offset removal first) | fix/null_root_a | arms in the mask (.48-.74) but the split is still nuisance (NMI <= .002), floor-merged; 1 leaf | 1 leaf, but now held by the floor (1.67), no longer by the guard | = a root (1 h 43 / 47) | NO-GO alone (no states; weakens the guard's null protection) |
| (b) sites over the whole support (islands) | offline on recB0 + replay | one site per arm (.50-.63, no cross-arm); the existing site layer finds NO pole on any arm site (endpoint z 0.4-1.8) | arm sites too (arms move: correct) | offline minutes; replay 1 h 12 / 31 | GO as the discovery rule, but only with site-focused classification |
| (c) rest subtracted + local re-alignment | fix_c_a1/a2, null_c_a1/a2 | poses scatter (flip30 .75-.77, 83-85 deg), a1 NMI .128 -> .000 | (same scatter) | + 25 min / site subtraction, + a root per site | NO-GO at 50 kDa (fixed-pose subtraction untested: needs a fixposes build) |
| (d) local first: a focused root per site (MASK_FILE = (b) site) | fix_m_a1..a5, null_m_a1/a2 | a1: NMI .128 (groups ~71 % pure), a2: .014, a3-a5: pending; all floor-merged | per-arm splits PASS test (z 14 / 12) + guard, track nothing (eta2 <= .003): accepted -> null FAILS (65 / 59 %) | 1 h 32 / ~40 core-h per site (5 in parallel ~200 core-h) | NOT YET: the merged-reference noise-lock certifies nuisance partitions inside a site |
| (e) effect floor inside the site | offline on all stage-B pairs | real arm pairs: in-site P/N 55-57, relative (site / rest) 108-153 | nuisance pairs in-site 2.3-68 (absolute floor passes nulls); null per-arm pairs relative 73-78 | seconds | absolute: NO-GO; relative: separates whole-root nuisance (<= 2.1), NOT the site-level noise-locks |

## RECOMMENDATION
- Port nothing now (design freeze; no change shows fixture benefit with the null held).
- The blocker for IgM-like states is not the mask or the floor alone: inside a site, stage B's merged-reference classification produces reproducible-looking splits on the null (known open issue, memory tree-good-pose-null-fails-halves-and-seeding). A cross-validated stage B (each half classified against the other half's references, per-half seeds) is the prerequisite; then re-run this battery with (b) + (d) + relative (e) — the expectation to test: null per-arm splits fall to noise, the fixture's a1 survives.
- Separately worth proposing (input-free, no constants), each through the RNase P closed-state and null gates: (b) island discovery (stock on diffuse variance, one site per arm here); and the background-offset removal for the cheap-children mask (on fixture AND real IgM the root mask covers everything because the cheap latent sorts particles by background level).
- Drop (c) at this mass and the absolute (e).
- Real IgM caveat: on real tb100k the recB0 variance gives no arm island (core only, smeared arms; registration still open), so (b) needs consistent arm poses on real data before it can matter.
- null_c_a1 (null, rest subtracted + focused root): split PASS z 18.7, guard kept 2.73, floor P/N 1.03 -> merged, 1 leaf: the subtracted path also certifies a null split up to the floor.

## Cross-validated stage B (user-approved next step, 2026-11-05 ~21:40): Python prototype (scripts/xstageb.py, xsb.sbatch)
- Scoring in Python, reconstructions by SIMPLE reconstruct3D (per-state even / odd maps from the current labels): per particle the residual r_i = x_i - a_i PFULL_i (SIMPLE consensus projection, LSQ scale), scored only on the SITE-MASKED state deviations M (R_k - mean) projected by a Fourier-slice projector with CTF + shift; whitened Gaussian log-likelihood (noise power per shell from the residuals), hard assignment, 12 A disk.
- Projector validated against SIMPLE's own noise-free projections of the consensus: median corr 1.000 on both datasets (SIMPLE: euler2m = Rz(-e3)Ry(-e2)Rz(-e1), fproject samples [h,k,0].M -> Rz(e1)Ry(e2)Rz(e3)(h,k,0); image origin at the box centre; +shift; +CTF). Wrong conventions gave .02-.31.
- Two designs (seeded per half, random labels, 6 iterations):
  - indep = each half classified against ITS OWN half's state maps (gold standard; states matched between halves at the end);
  - swap = each half classified against the OTHER half's maps (as approved). Concern: half A's selected particles inherit half B's map noise, which is also in half B's own test maps -> may still lock; the null decides.
- Runs at the a1 site: fix_a1_indep 54971698, null_a1_indep 54971699, fix_a1_swap 54971700, null_a1_swap 54971701.
- ~22:00 delegated (user asked), notes in the shared folder side-session-notes/local_c1bafd06-1da5-43a7-8173-c2f4899054f5/ (cheap_offset.md, island_gates.md, discrete_vs_continuous.md): cheap-offset (local_aa5835a1: root-mask background-offset mechanism + fix, gate datasets), island-gates (local_8e8f5add: island rule on RNase P closed lineage / nulls / mm3q / continua), discrete-vs-cont (local_4348f89c: oracle study + input-free discreteness statistic at a site). This session keeps the cross-validated stage B, per-half seeds and the scorecard.

### a1 results (2026-11-06; re-scored after a half-alignment fix)
- BUG found and fixed: the halves are seeded separately, so half 1's state 1 can be half 0's state 2. All four runs ended anti-aligned (cross-half deviation corr -0.60 to -0.85), which flipped the sign of the cross-half r, swapped the P/N "signal" and "noise" terms and zeroed the pooled NMI. In swap it also caused a 2-cycle: each iteration relabelled 50-70 % of the particles, because each half takes the other's state numbers. Fix (xstageb.py v2): before every step and in final, half 1's states are renumbered to half 0's (Hungarian on the cross-half correlation of the site-masked state deviations; input-free). The scoring also runs on a 16-process pool (~6 min -> ~1 min per step).
- Aligned results, a1 site, 6 iterations, seed 1 (P/N here = power of the half-averaged group difference over the power of the half difference, ~1 for noise; the tree's own pair_maps P/N from proto_pn.py is pending):

| run | cross-half r / z | P/N site / rest of molecule / relative | per half | groups vs truth |
|---|---|---|---|---|
| fix indep | .766 / 19.2 | 7.33 / 1.66 / 4.41 | a1 NMI .089 / .061 | all a1 NMI .074, groups 66 / 34 % pure; other arms .000 |
| fix swap | .847 / 23.7 | 11.87 / 1.48 / 8.03 | a1 NMI .108 / .079 | all a1 NMI .093, 68 / 32 %; other arms .000 |
| null indep | .595 / 12.5 | 3.84 / 3.64 / 1.06 | a1 angle eta2 .032 / .026 | eta2 .029, group mean angles -2.4 / +2.5 deg |
| null swap | .705 / 16.0 | 5.78 / 3.94 / 1.47 | eta2 .041 / .035 | eta2 .038, -2.8 / +2.9 deg |

- The noise-lock is gone: the null's merged-reference focused root tracked nothing (eta2 .001), its cross-validated split tracks the hinge angle (eta2 .03-.04). The expectation in the RECOMMENDATION ("null per-arm splits fall to noise") is wrong: the null's arm is a REAL continuum, so a cross-validated split of it is a real, reproducible median-like cut (indep, whose cross-half test shares no noise: z 12.5). A reproducibility test cannot hold the null twin; holding it needs a discreteness call.
- The fixture pays some power for cross-validation: a1 NMI .074-.093 (66-68 % pure), vs .128 (71 %) for the merged-reference focused root.
- The relative P/N separates the two at a1 (fixture 4.4-8.0 vs null 1.06-1.47): the null's group difference has as much power in the rest of the molecule as in the site. Checking with the tree's own P/N and a truth-mask breakdown (xsb_pn.sbatch, diag_pn.py) whether that is the arm outside the site (mask coverage) or something else.
- Swap caveat: in swap, half 0's labels come from half 1's maps, so half 0's selected particles carry an imprint of half 1's noise, which half 1's own maps also contain. The final cross-half test can therefore be inflated. indep's test is clean. Control: both designs on a rigid-core patch of the same size (no conformational change there).

### TRAP found on the way: the tree's effect-floor P/N is a lottery past the maps' resolution (pair_maps' FSC match is unbounded)
- `pair_maps` (src_igm `simple_state_mergeup.f90:5208`, mirrored by proto_pn.py) matches group 2 to group 1 per shell with g = 0.5 (c14 sqrt(pw1 pw4) + c23 sqrt(pw2 pw3)) / (fb sqrt(pw3 pw4)), fb = group 2's own half FSC, with no bound. In the 10-40 A band it then forms ex = max(smoothed(res^2) - smoothed(nse^2), 0), p = sum(ex) (the effect floor) and the top-5 % share (the locality guard).
- At shells past the maps' resolution fb is noise (|fb| < .02 at 11-13 A on these 50k-particle group maps), so g is a ratio of two noise numbers: 39.3 / 42.3 / 56.8 / -73.9 at single shells of four a1 pairs (median .95-.98). That shell's group-2 noise enters res and nse amplified ~40-70x.
  - The unclipped sum stays near 0 (zero mean).
  - The clipped sum the tree uses gains a large positive bias that varies from run to run. a1 pairs, clipped vs unclipped P/N:
    - over the sphere: fix indep .36 / .12, fix swap 6.0 / -.08, null indep 8.9 / -.22, null swap 32.2 / .66;
    - in the site: null indep 18.1 / .30, null swap 23.7 / -3.96.
  - Only fix indep (max g 3.1) escaped.
- Consequence for the tree: the effect floor (P/N >= 2) and the locality guard can PASS a split on noise whenever one shell of the 10-40 A band falls past the pair's resolution with |fb| near 0. Small groups (low resolution) are the most exposed. The fixture's floor values above (fix_m_a1 0.24, etc.) happened not to blow up.
- Input-free guard (prototype, proto_pn.py GMODE=sig): match only where fb > 3 / sqrt(N_shell / 2), the 3-sd FSC noise level of the shell, else g = 1 (no matching, so res and nse stay balanced). Not ported; for the tree-distribution owners.
- RE-MEASURED at the gate nodes (2026-11-06, pair-maps-guard session, pair_maps_guard_gates_2026_11_06.md in this folder): the guard is INERT at all 69 gate nodes (231 pairs, bc / df / base / guardfix arms). The tree's stage-B group maps resolve the whole 10-40 A band (fb >= 4.97x the bar, |g| <= 1.15), so no floor or guard verdict changes and the margins are real (r131 corner 3.08 at seed 2, 2.63 at seed 1). The 39-74 g values came from this prototype's box-200 xsb maps, which reach shells the tree's box-164 stage B does not. The trap stays a latent risk for low-resolution pairs, not a current gate problem. The same session notes that on the null roots the floor alone would keep the split (P/N 6-41); only the locality guard merges them.
- With the guard (GMODE=sig), clipped and unclipped agree, and the unclipped "null in-site ~0" above turns out to be the same amplified noise (zero-mean but huge variance), not a real zero. Guarded P/N per a1 pair:

| a1 pair | sphere (clipped / unclipped) | site | rest of the molecule | relative site / rest |
|---|---|---|---|---|
| fix indep | .35 / .12 | 15.4 | .33 | 33.2 |
| fix swap | .37 / .17 | 20.5 | .38 | 41.9 |
| null indep | .32 / .09 | 12.4 | .93 | 12.8 |
| null swap | .36 / .16 | 17.5 | .99 | 17.2 |

- Read: the null's continuum cut is a real in-band difference in the site (12-17, fixture 15-21), so the in-site floor cannot hold the null. Every pair is far below the sphere-wide floor (.3-.4 < 2). The relative floor separates by ~2.5x (null 13-17, fixture 33-42), but both are far above the existing constant 2; using it would need a new, calibrated threshold.
- The 15 A version (xstageb final: fixture 4.4-8.0 vs null 1.06-1.47) separates better because the null's groups also differ at low resolution in the CORE (P/N 4.1-4.5 there, fixture 1.3-1.5; diag_pn.py): part of the null cut follows a global nuisance. The battery's quality eta2 tests that.

### BATTERY (2026-11-06 00:00-00:35): cross-validated stage B at every (b) arm site, fixture and null (xsb_battery_submit.sh, jobs 54979236-49)
Swap design (each half scored against the other half's site-masked state maps), per-half random seeds, half 1 renumbered to half 0 before every step, 10 iterations, K = 2, the 100k particles at their core-locked input poses. Scored by xsb_battery_score.sbatch (54979716). Cost: 20-31 min per site on 16 cores (5.5-8.3 core-h); 14 runs ~90 core-h. Quality eta2 (bfac, snrq) = 0.000 for EVERY split below.
- Columns: r / z = cross-half test of the site difference; "15 A" = xstageb final P/N (site / rest of the molecule / relative); "tree" = pair_maps P/N with the guarded FSC match (site / rest / relative); "core" = diag_pn's P/N in the core outside the site (15 A).

| site | r / z | truth: own arm (per half) | 15 A site / rest / relative | tree site / rest / relative | core |
|---|---|---|---|---|---|
| fix a1 | .964 / 38.1 | NMI .146 (.148 / .144), groups 72 / 28 % pure | 53.6 / 1.77 / 30.3 | 32.4 / .69 / 43 | 1.49 |
| fix a1 seed 2 | .934 / 32.1 | NMI .138 (.146 / .130) | 28.8 / 1.67 / 17.3 | 30.2 / .59 / 45 | 1.35 |
| fix a2 (3 states in plane) | .874 / 24.1 | NMI .047 (.055 / .041); -30/0/+30 fractions .31/.43/.25 vs .59/.27/.14 | 14.6 / 2.61 / 5.61 | 30.0 / 1.22 / 24 | 3.03 |
| fix a3 | .195 / 3.6 | NMI .014 (.021 / .008) | 1.47 / 1.14 / 1.29 | 7.3 / .04 / - | 1.15 |
| fix a4 | .679 / 15.2 | NMI .007 (.013 / .003) | 5.23 / 8.30 / 0.63 | 20.2 / 2.38 / 8.4 | 9.13 |
| fix a5 (3 states out of plane) | .674 / 14.4 | NMI .039 (.057 / .026) | 4.93 / 2.29 / 2.16 | 17.9 / .48 / 32 | 2.58 |
| null a1 | .887 / 25.6 | angle eta2 .069 (.070 / .067), group means -3.7 / +3.9 deg | 16.7 / 6.30 / 2.64 | 30.1 / 1.96 / 15 | 7.44 |
| null a1 seed 2 | .820 / 21.0 | eta2 .061 (.071 / .052) | 10.1 / 4.57 / 2.20 | 19.5 / 1.27 / 15 | 5.41 |
| null a2 | .162 / 3.0 | eta2 .021 (.064 / .001) | 1.38 / 2.54 / 0.54 | .07 / .48 / - | 2.67 |
| null a3 | .497 / 10.1 | eta2 .038 (.059 / .022) | 2.95 / 2.98 / 0.99 | 11.8 / .69 / 17 | 3.28 |
| null a4 | .662 / 14.7 | eta2 .000 | 4.84 / 12.7 / 0.38 | 15.7 / 2.87 / 5.5 | 13.9 |
| null a5 | .815 / 20.2 | eta2 .062 (.060 / .064), -4.3 / +4.4 deg | 9.81 / 3.97 / 2.47 | 24.9 / 1.06 / 22 | 4.74 |
| CONTROL null rigid-core patch, swap | .899 / 26.7 | every arm eta2 .000 | 18.5 / 20.5 / 0.90 | 62.6 / 23.8 / 2.63 | 20.4 |
| CONTROL same patch, indep | .889 / 25.8 | .000 | 16.5 / 19.7 / 0.84 | 59.4 / 23.4 / 2.54 | 19.6 |

- Per-site products (all five sites, 32 groups):
  - fixture: combo NMI .053 (per arm a1 .066, a2 .021, a3 .006, a4 .003, a5 .018). The merged-reference focused roots gave .042 (a1 .128, a2 .014, a3 .004, a4 .043, a5 .014).
  - null: largest group 3.6 % (FAIL; gate >= 80 %).
- Classification: with in-loop half alignment the swap design converges well (fix a1 cross-half deviation corr .16 -> .96 over the iterations). It is the best site classifier tried here (a1 .146, a2 .047, a5 .039). The fixture's a3 and a4 arms (two states at +-20 deg, .65/.35 and .70/.30) are not found (NMI .014 / .007).
- The rigid-core controls settle the swap worry: swap and indep give the same split (z 26.7 vs 25.8), so the swap test is not inflated by shared noise. But a site with NO conformational change still splits at z 26 with in-patch P/N ~60, following no arm and no quality draw. The cross-half test is no gate at all.
- Gates tried on these splits (a pass splits the node):
  - cross-half z >= 3: passes everything (fix a3 3.6, null a2 3.0 marginal).
  - sphere floor (guarded tree P/N >= 2): rejects everything (.27-1.06), the fixture's real a1 / a2 / a5 included.
  - tree relative >= 2: passes everything, the rigid-core controls included (2.5-2.6).
  - 15 A relative >= 2:
    - fixture: a1 30 / 17, a2 5.6 and a5 2.16 pass; a3 1.29 and a4 .63 are missed;
    - null: a1 2.64 / 2.20 and a5 2.47 PASS, a2 .54, a3 .99 and a4 .38 held;
    - controls: rejected (.84-.90).
    The null's two most strongly moving arms (+-4 deg group means) pass, so the null fails: 4 groups, largest ~25 %.
- Bottom line of the battery: cross-validation removes the noise-lock, and the 15 A relative floor removes global splits (rigid-core controls, both a4 sites). What remains is the null's REAL continuum cuts at a1 / a5, which pass every map-level gate the fixture's real arms pass. Only a discreteness call can separate them. The sibling session's GAPw is being applied to these runs.
- What the rigid-core split follows (pose_split.py, truth vs input poses; input error median 1.43 deg / 4.0 A): the particles' SHIFT ERROR. R2 of the label on the shift error .32 / .33 (swap / indep); rotation error .01; view and defocus .000. On the arm sites R2 is <= .035 for every run; it is highest for the two a4 splits (.023 / .035), which track no arm. Mechanism: with nothing moving inside a site, a site-focused classifier learns two copies of the site displaced by +-(sub-pixel to 1 px). Both halves share the shift-error distribution, so the split reproduces (z 26) and is global (15 A relative .84-.90). Real data have the same sub-pixel shift errors.
- GAPw on the battery (discrete-vs-cont session, discrete_vs_continuous.md section 6; its constants unchanged):
  - fixture a1 DISCRETE in both seeds and both readouts (six-bin .09-.10; direct .08-.10);
  - every null run and both rigid-core controls: no false DISCRETE. The controls read clearly not discrete (1.45-1.54); the null arms are "no call": fewer than 10 voxels between the lobes (a1, a4, a5) or the factor gates (a2, a3);
  - fixture a2 / a4 / a5: no call (B); a3: no call (G1).
  - The .28 threshold decides nothing here (calls at .08-.10 vs 1.45-1.54). The outcome rests on the geometric gates (B >= 10 voxels between lobes, L >= 10 A, lobe axes within 30 deg, r1 >= .5).

## SCORECARD UPDATE (2026-11-06; supersedes the (d) and (e) rows above)
| change | fixture | null twin | cost | verdict |
|---|---|---|---|---|
| (d) + cross-validated stage B (swap + per-step half alignment) | best site classifier: a1 NMI .146 (seed 2 .138), a2 .047, a5 .039; a3 / a4 not found (.014 / .007); five-site product combo .053 (merged-reference .042) | REAL cuts of the continuous arms (a1 / a5 eta2 .06-.07, a3 .04); product 3.6 % largest group (FAIL without a further gate) | ~25 min, 6-8 core-h per site (16 cores; Python scoring) | GO as the per-site classifier (the noise-lock is gone); NOT sufficient alone |
| (e) site floors | 15 A relative >= 2 keeps a1 / a2 / a5, misses a3 / a4; the tree-band (10-40 A) relative keeps everything | 15 A relative passes a1 (2.2-2.6) and a5 (2.47), rejects the rigid-core shift-error splits (.84-.90); tree-band relative passes all, controls included (2.5-2.6); sphere floor rejects everything, the fixture too | seconds | 15 A relative: a useful global-split filter (shift error, nuisance), not a discreteness gate; tree-band versions: NO-GO. pair_maps' FSC match needs a guard (TRAP section) |
| discreteness: GAPw (sibling) on the cross-validated splits | a1 DISCRETE (.08-.10); a2 / a5 (3 states) and a3 / a4 no call | no false DISCRETE (9 runs); controls not discrete | ~8 min per site (8 cores) | GO as the split condition, with its limits: misses 3-in-a-row states and lobes closer than ~10 voxels; needs >= ~50k particles at a 50 kDa arm; its constants come from this fixture |
| COMBINED: (b) islands + cross-validated stage B + split only when GAPw says DISCRETE | a1 split (72 / 28 % pure), the other four arms unsplit | 1 leaf (all arm sites and both controls held) | ~35 min, ~10 core-h per site | the FIRST configuration in this study with a fixture gain AND the null held. The gain is one arm of five, at modest purity. |

## RECOMMENDATION (2026-11-06; supersedes the one above)
- Still port nothing. The combined path is the first with a fixture gain and the null held, but it has not seen the RNase P closed-state gate or the null fixtures (qtoynull / micnull / olenull). Its constants are this fixture's (GAPw), and it recovers one of five arms.
- Next step, if approved: a gate plan for the combined per-site path as an ADDITION to the current tree (it can only add splits, and only on DISCRETE). Replay it on saved nodes (Python, no Fortran):
  - the r131 closed-corner nodes, both seeds, and the 692k closed-lineage nodes;
  - the three null fixtures;
  - qtoy2 / mm3q.
  Island-gates already shows the island rule adds sites only at r131 seed-1 r1 and qtoy2 r4 (692k reruns pending).
- Real IgM is not ready for any of this: real tb100k has no arm islands (smeared arms; registration still open), so (b) has nothing to work on there yet.
- Separate, tree-wide: pair_maps' FSC match (TRAP section) is unbounded; with a 3-sd guard on fb, re-measure the effect floor and guard at the gate nodes (r131 corner 2.6-3.2) before trusting their margins.
- Dropped: (a) mask per region, the cheap-children offset removal (cheap-offset session: NO-GO, the mask is a symptom of a root latent with no arm signal), (c) subtraction plus local re-alignment, and the absolute and tree-band (e) floors.

## GATES: the combined per-site path as an add-only layer on saved nodes (2026-11-06, user-requested; Python, no Fortran)
What was replayed, at every site of each saved node (plus its islands):
- the cross-validated stage B (scripts/xsbg.py: swap design, per-half seeds, half alignment every step, 10 iterations, K = 2), then GAPw (sibling, discrete_vs_continuous.md section 7, made dataset-generic: maps Fourier-cropped to ~5.1 A/px, constants in A);
- the split rule: split only on DISCRETE.
Saved nodes:
- node projects reference the full 932k stacks, so each node is first exported to its own stack;
- poses, CTF and the even/odd halves are kept on re-import;
- the projector conventions are refit per node.
Positive control through the same path: the fixture's a1 site is DISCRETE (GAPw .09), as in the battery.

Wave 1: 56 runs, jobs 54980902-05 and 54980988-1039 (mean 8 min, max 27 min on 16 cores, ~130 core-h with the preps). Tables:
- gates_x/runs/*/ (readout: scripts/gate_table.sh);
- GAPw per run: discrete_vs_continuous.md section 7.

| gate | nodes (particles, sites) | runs | DISCRETE | what the per-site splits do | verdict |
|---|---|---|---|---|---|
| null fixtures | qtoynull root (60k, 20), olenull root (81k, 4, c2); micnull root has no site and no island, so there is nothing to add | 24 | 0 | every run fails GAPw's first gate (z1 <= 3.4, r1 <= .57); the one six-bin r1 >= .5 (qtoynull s07) fails G2 / G3 / B (lobes 4 A apart) | PASS: both nulls stay one leaf |
| RNase P closed corner | bc_r131 seed 2: r (44k, 4), r1 (25k, 5), r2 (19k, 2); seed 1: r (44k, 2), r1 (23k, 2 + 1 island), r2 (21k, 4) | 20 | 0 | closed reference rows never separated: the most enriched group holds 1.39x the node's share (r131s1 r2 site 4), usually .9-1.2x. The closed-arm occupancy of a run's two group maps differs by at most .19 (.81 vs 1.00), at r131 root site 2 in both seeds and r1 site 1 | PASS: nothing is added; the tree's own corner split is untouched |
| 692k closed leaves | r1312 (12k, 3), r1321 (12k, 4), r2112 (13k, 5) | 12 | 0 | both groups keep the leaf's closed arm (r1312 1.07 / 1.07, r1321 .92 / .93, r2112 .91 / .92); closed rows split evenly | PASS: no closed-substate split |

(The leaf runs are tagged _iNN by a submit-pattern slip; they are the leaves' own tree sites, not islands.)

Without GAPw, the other gates would have split nulls and closed leaves:
- every run has tree-band P/N >= 5, so any effect floor passes all 56;
- 15 A relative >= 2 passes 7 of 20 qtoynull sites (up to 3.97), olenull s04 (2.01), r131s1 r2 s04 (3.24), r131s2 r s04 (2.61), and lf_r2112 s02 / s04 (5.20 / 2.12; a closed-substate split);
- cross-half z >= 3 passes r131s1 r s01 / s02, r131s1 r2 s04, lf_r2112 s02 and qtoynull s07.

The strongest real-data factors are reproducible but carry no two-state motion: r131s1 r s01 (r1 .61-.64), r131s1 r2 s04 (.65-.71) and lf_r2112 s02 (.73-.80). Their lobes are 7-8 A apart with no density between them, like the rigid-core shift-error control, and G2 / B stop them.

What the gate does not show:
- Sensitivity. The RNase P gate nodes hold 12-44k particles. GAPw needed >= ~50k to call the fixture's own a1 arm (at 25k it read .32, a miss). So these nodes show the layer is safe, not that it would find a real discrete RNase P motion.
- The closed state is not a per-site motion here. Per-site stage B never isolates it, so at r131 it is the tree's own stage-B corner, which spans sites.
- Untested: the 692k root (4 sites; a 692k Python replay).
- The other ancestors of the closed lineage are wave 2 below.

Wave 2 covers the rest of the 692k closed lineage: 33 runs, jobs 54982741-57 and 54982828-43 (mean 9 min, ~87 core-h with the preps).
- Nodes:
  - the closed leaf r1311;
  - the ancestors r13, r131, r132, r21 and r211.
- Site masks come from the island-gates reruns, each node rerun as a root.
- The r13, r131 and r21 reruns stopped before stage B. For those three, the replay takes the production node's stage-B project together with the rerun's sites.
- GAPw per run: discrete_vs_continuous.md section 7b. 7b covers all 33 runs; lf_r1311_s01, the last to finish, is no call (fails G1-G3).
- "Agreement" below compares a site's two groups with the tree's own child split of that node: production rows.txt of the children against the run's labels_final.txt.

| node | particles, closed share | sites | DISCRETE | what the per-site splits do |
|---|---|---|---|---|
| r1311 (closed leaf) | 12.5k, .060 | 6 | 0 | All fail G1 (z1 <= 2.8). Closed rows 44-56 % per group. Closed arm .76-.82 in every group (difference <= .06). No closed-substate split. |
| r131 | 25k, .068 | 3 | 0 | All fail G1 (r1 <= .42). Closed rows split evenly. |
| r132 | 22k, .057 | 7 | 0 | All fail G1 (r1 <= .44). Nearest call: s01 (GAPw .30-.32), which also fails G2. |
| r13 | 47k, .063 | 5 | 0 | s02 is reproducible (z1 6.9-8.1, r1 .58-.65), but its lobes are 6-7 A apart (G2). No site follows the tree's r131 / r132 split (agreement <= .61, ARI <= .04). |
| r21 | 53k, .029 | 8 | 0 | **s01-03** are strong (z1 9-17, r1 .82-.93) and are the tree's own r211 / r212 split: agreement 82-86 %, ARI .41-.51. The r211 side holds ~90 % of the closed rows (1.7-1.8x) at closed arm .71-.73; the other side reads open (.28). GAPw stops them at G2 / G3: lobes 6-8 A apart, half axes 22-74 deg. **s06 / s08** are reproducible but fail G2 and do not follow r211 / r212 (ARI .02). |
| r211 | 27k, .052 | 3 + 1 island | 0 | **s01** is the strongest factor of either wave (z1 23-25, r1 .97-.98, lobes 14 A, half axes 2-6 deg). It passes every gate, then reads not discrete (GAPw .45). It is 83 % the tree's own r2111 / r2112 split (ARI .45), the split that makes the closed leaf r2112: closed arm .88 vs .56, 90 % of the closed rows in one group. The island fails G1. |

Both waves together: 89 runs, 0 DISCRETE.

| gate | runs | DISCRETE | verdict |
|---|---|---|---|
| null fixtures (qtoynull, olenull; micnull has no site) | 24 | 0 | PASS: each null stays one leaf (100 %) |
| RNase P closed corner (bc_r131, both seeds) | 20 | 0 | PASS: the tree's corner split is untouched |
| 692k closed lineage: leaves r1311, r1312, r1321, r2112; ancestors r13, r131, r132, r21, r211 | 45 | 0 | PASS: no split added anywhere; no closed leaf split into substates |

What this means for a port:
- Both gates pass. The layer adds nothing on the nulls or on the closed lineage. It is also inert on RNase P: its only gain so far is the fixture's discrete a1 arm, and only at >= ~50k particles.
- It must stay add-only. On the closed lineage, the strongest per-site factors (r21 sites 1-3, r211 site 1) are the tree's own splits seen through one site, and GAPw calls none of them DISCRETE. With GAPw as the only split rule, the tree would stop at r21 and never reach the closed leaf r2112. Why they read this way:
  - at r21 the lobes are close and their axes unstable;
  - at r211 the density between the lobes is filled.
  The open-to-closed change is a broad range here, with closed at one end, not a two-state switch at one site.
- The leaf results rest on low power. At 8-13k particles GAPw cannot call even the fixture's arm, so the closed leaves pass by lack of power as much as by design.
- The untested risk is where GAPw has power on RNase P:
  - the 692k root: 4 sites; needs a 692k Python replay;
  - r5 (120k).
- Not covered:
  - the closed leaves r1322, r5111 and r5121, and the r5 lineage: production prunes site masks and no island-gates rerun covers them, so each needs a site-discovery rerun first.
- Real IgM: registration is closed (user, 2026-11-06; this supersedes "registration still open" above). The tree input is the 683k consensus as is:
  - ~1/3 of particles sit in the wrong 60-degree register;
  - the J slot is ~60 % filled.
  That error is itself a two-state contrast at the J slot and the arm sites. So on real IgM, a DISCRETE call there must be checked against register before it counts.
- Recommended next, before any port proposal: the generalization check on the saved 60k fixture roots. Cost: 3 preps + 7 site runs, ~20 core-h, under 1 h wall. Expected readings:
  - ctun (uniform continuum): must read no DISCRETE;
  - ctpk (dwell continuum): may read DISCRETE;
  - mm3q (three independent two-state domains): should read DISCRETE at its motion sites.
  - ctpk and mm3q stopped before stage B, so their preps take the root's input project.
- If mm3q's motion sites do not read DISCRETE, the layer's value is one fixture arm and porting it is not worth it. If they do, replay the 692k root next, before a port proposal.

## PORT PROPOSAL (2026-11-06, user-requested; the generalization check skipped at the user's call)
port_proposal_site_split.md in this folder: the layer acts only at leaves the tree's own tests stop above the depth cap (12 of the 66 f692 leaves; no closed leaf); islands are for the layer only; the direct GAPw readout; off by default (site_split=no). About 1,300 lines on fork/hetero-tree 72e5aaba0. Gates 0-1 come after approval; ctun / ctpk / mm3q run for the first time at gate 5.
