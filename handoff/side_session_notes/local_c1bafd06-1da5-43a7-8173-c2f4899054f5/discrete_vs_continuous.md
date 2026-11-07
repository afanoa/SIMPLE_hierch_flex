# Discrete vs continuous at a site: IgM-like fixture (2026-11-05)

Work dir (batch): `/mnt/beegfs/elmlund/afan/igm_cryolo/tree/fixture_igm/discrete/` (scripts/, labels/, rec/, crops/, out/, logs/, xsplit/). Python only, no Fortran, nothing committed.

## Bottom line

1. **Bin-map rank / shape tests cannot certify discreteness at this SNR.** The affine residual and the number of reproducible bin-map factors separate a two-state mixture from the continuous hinge only with near-perfect orderings. Per-bin peak positions do differ at d' 1.1 (the fixture's jump between two fixed positions, the null's stay central), but they are path-based and inherit point 5. On the null, the second factor's cross-half z is 13-23 with exact-angle bins, 8.5 at d' = 4, 1.9 at d' = 2 and 0.3 at d' = 1.1, and d' ~ 1.0-1.1 is what any classifier reaches here (focused root NMI .128 = 71 % purity; the cross-validated site score below: corr .42-.46 with the a1 state). These tests can reject "two states" with a strong ordering, never confirm them. On a noise-locked ordering they also fail the other way: the fixture's merged-reference evidence bins read z2 4.9 / affine z 5.0, i.e. "continuum".
2. **What separates: GAPw, the consensus density between the two lobes of the site's first factor.** It is path-free and measured on the node consensus, so the ordering only has to place the lobes, and the lobes must come from a signal-dominated factor (cross-half r1 >= .5).
   - Fixture a1: .02-.09 (oracle bins), .02 (noise-locked evidence bins), .20 / .23 / .32 (input-free at 100k / 50k / 25k: needs >= ~50k).
   - Fixture a3 / a4 (two states +-20 deg, populations .65/.35 and .70/.30; oracle only): .13-.15.
   - Null: >= .37 on every gated bin set (26 sets: four arms, exact-angle to d' = 1.1 orderings, 5.6-20k per bin).
   - The threshold is .28, taken from the null.
   - Without the r1 gate a weak null factor (25k, z1 3.5) reads .21, a false call.
3. **It misses a2.** The fixture's a2 has three states in a row (-30/0/+30); the middle state fills the gap between the outer lobes (oracle GAPw .36-.62), so a2 reads "not discrete". Input-free it is no call (GAPw 1.34 on a 5-voxel gap region). That is a miss, not a false positive. Three-state vs continuum needs the third factor (z3), which shows only with exact orderings.
4. **The ordering must be cross-validated.** Each half is ordered by the other half's site direction. A merged-reference (stage-B) ordering of the null gives a reproducible fake first factor (z1 38.9). It was stopped here only because the locked pattern has no dipole (lobe half-separation 2 A); it is not a guarantee.
5. **1-D profile dips along the factor path are not robust.** On the null, a 10 deg path tilt gives a false dip in up to 4 of 12 random directions. Input-free paths differ by 12-18 deg between halves, and the input-free fixture's dip vanished (.015, z 0.5). Do not use them.

## Setup

- Sites (box 200 masks, resampled to 128): fixture a1 = `proto/fix_root_stock_recB0_bgyes/sites/site_06.mrc`, a2 = `site_02`; null a1 = `proto/null_root_stock_recB0_bgyes/sites/site_03.mrc`, a2 = `site_04`. The site covers .57-.62 of the arm's fan mass.
- Bin maps: SIMPLE `selection` on a copy of `<ds>/src.simple` (core-locked input poses) + `reconstruct3D nstates=K mskdiam=600 pgrp=c1 box_crop=128` (5.11 A/px; ~1 min per 100k), `_unfil` half maps, 15 A Gaussian low-pass (as mix_test.py), cut to a box around each arm (`crops/`).
- Oracle bin sets (`dc_setup.py`; truth builds bins only):
  - truth atoms (fixture: arm state x 12 replicates; null: 5-deg angle groups x 3), combined into graded mixtures / angle bins;
  - K = 6 quantile bins of a noisy truth score s = d' g(angle) + N(0,1), where g is the arm's 3D discriminability profile (from the noise-free component maps), so d' = the state separation in noise sd; d' = 4, 2, 1.1 and exact angle;
  - quality controls: bins by particle B, amplitude (`quality.txt`) and input-pose angular error (truth vs perturbed).
- Each bin set is also scored at the arm it was not made for, as a no-signal control.

### Statistics (all cross-half; every nuisance regressor is built per half from that half's own maps)

Shared regressors carry both halves' noise and drive cross-half z negative: the first nuisance version gave affine z down to -17.

- **FAC**: factors of the centred bin maps in the site (constant removed; `nuis`: also each half's consensus V, its Laplacian (B blur) and its rotational blur sum_a L_a^2 V). Even-half SVD loadings are applied to the odd half and back; z per factor = Fisher z of the cross-half map correlation over n_eff (mask voxels x (4/3) pi (5.11/15)^3 / 2, ~370 at a1).
- **AFF**: middle bins fitted from the end bins (plain: affine, as mix_test.py), max cross-half residual z.
- **TUBE**: cross-validated 1-D dip of the consensus along the segment joining the first-factor lobes (25 A tube). Rejected, see below.
- **GAPw** (proposed): see section 3.

## 1. Oracle table

a1 site (fixture: 2 states +-25 deg out of plane, .5/.5; null: uniform over the same range). n/bin = particles per bin.

| bins (truth) | n/bin | FAC z1 | FAC z2 | AFF z | TUBE dip (z) | GAPw |
|---|---|---|---|---|---|---|
| fix, p(state 2) 0/.5/1 | 8.3k | 25.7 | 0.7 | 0.8 | .14 (2.0) | .05 |
| fix, 0/.5/1 | 16.7k | 31.3 | 0.3 | 0.3 | .16 (4.3) | .05 |
| fix, 0/.5/1 | 33.3k | 38.1 | 0.3 | 0.3 | .15 (4.7) | .04 |
| fix, 0/.25/.5/.75/1 | 16.7k | 33.6 | -0.2 | 0.0 | .18 (5.3) | .05 |
| fix, noisy score d' 4 / 2 / 1.1 (K 6) | 16.7k | 41.4 / 36.6 / 28.2 | 0.0 / -0.7 / 0.2 | 0.6 / 0.5 / 1.3 | .12 (3.9) / .14 (4.9) / .17 (5.4) | .05 / .04 / .06 |
| null, exact angle K 3 | 6.7 / 13.4 / 20k | 20.7 / 26.3 / 29.9 | 12.3 / 17.6 / 20.4 | 12.9 / 18.1 / 21.2 | .03 / .01 / -.04 | .48 / .45 / .42 |
| null, exact angle K 5 | 6.7 / 13.3 / 20k | 23.2 / 29.3 / 33.5 | 12.5 / 17.7 / 20.5 (z3 6-11) | 12.9 / 18.1 / 21.2 | .04 / -.01 / .00 | .51 / .45 / .44 |
| null, exact angle K 6 | 16.7k | 33.8 | 20.9 (z3 12.8) | 20.2 | -.03 | .44 |
| null, noisy score d' 4 / 2 / 1.1 | 16.7k | 26.9 / 18.8 / 10.2 | **7.1 / 1.7 / 0.3** | **6.2 / 1.5 / 1.0** | .03 / -.06 / -.03 | .47 / .52 / .55 |
| random bins (fix, null) | 16.7-33k | <= 0.2 | - | - | gated | gated |

a2 site (fixture: 3 states -30/0/+30 in plane, .45/.35/.20; null: uniform):

| bins (truth) | n/bin | FAC z1 | FAC z2 | FAC z3 | AFF z | GAPw |
|---|---|---|---|---|---|---|
| fix, 3 pure states | 8.6 / 17.1k | 23.9 / 30.8 | 17.2 / 22.3 | 0.0 / 0.0 | 17.8 / 23.4 | .36 / .41 |
| fix, graded 3-state K 5 | 16.4k | 31.1 | 20.9 | 0.3 | 22.7 | .50 |
| fix, noisy score d' 4 / 2 / 1.1 | 16.7k | 32.2 / 26.5 / 19.2 | 19.0 / 4.7 / -0.2 | -0.2 / -0.2 / 0.0 | 17.6 / 4.1 / 0.2 | .47 / .62 / .51 |
| null, exact angle K 3 | 5.6 / 11.1 / 16.7k | 19.9 / 25.4 / 29.0 | 12.7 / 18.2 / 21.7 | - | 13.4 / 19.1 / 22.7 | .42 / .40 / .39 |
| null, exact angle K 6 | 5.6 / 11.1 / 16.7k | 24.3 / 30.6 / 33.9 | 13.9 / 19.0 / 22.8 | 7.8 / 11.1 / 14.1 | 13.3 / 18.3 / 22.1 | .46 / .46 / .47 |
| null, exact K 6; d' 4 / 2 / 1.1 | 16.7k | 34.2; 26.0 / 18.8 / 11.3 | 23.0; 6.9 / 1.1 / -0.4 | 14.5; -0.3 / 0.4 / -0.8 | 22.4; 6.6 / 1.4 / 0.4 | .47; .50 / .54 / .45 |

Other arms (oracle, K 6, 16.7k per bin; sites identified by fan-envelope overlap: a3 = fixture site_04 / null site_06, a4 = fixture site_03 / null site_05). These test unequal populations and a smaller separation (40 deg):

| arm, bins | FAC z1 | FAC z2 | AFF z | GAPw tau .4 / .5 / .6 |
|---|---|---|---|---|
| a3 (out of plane +-20, .65/.35): fix d' 4 / 1.1 | 34.1 / 23.3 | 0.0 / -0.2 | 1.2 / -0.1 | .09/**.13**/.19, .10/**.15**/.22 |
| a3 null: exact / d' 4 / d' 1.1 | 33.0 / 26.8 / 11.3 | 20.0 / 7.4 / 0.1 | 18.8 / 6.7 / 0.8 | .26/**.38**/.65, .32/**.37**/.69, .12/**.42**/.61 |
| a4 (in plane +-20, .70/.30): fix d' 4 / 1.1 | 32.7 / 21.1 | 0.1 / -0.7 | 0.9 / 0.5 | .09/**.15**/.23, .09/**.14**/.24 |
| a4 null: exact / d' 4 / d' 1.1 | 32.5 / 26.6 / 11.5 | 18.9 / 7.2 / 0.4 | 17.5 / 6.8 / 0.7 | .35/**.41**/.73, .37/**.46**/.84, .18/**.83**/1.06 |

- At tau .5 the 2-state arms read .13-.15 against .37-.83 for their continuous twins.
- tau .4 gives false gaps on every weak-ordering null (d' 1.1: a1 .29, a2 .30, a3 .12, a4 .18), so use tau .5.

Second-factor share of the cross-half covariance (a1):

- null: .18-.26 with exact orderings, .08 / .03 / .015 at d' 4 / 2 / 1.1;
- fixture: <= .006 everywhere (noise-free expectation 0.000).

Per-bin density peak along the path (oracle path, 25 A tube):

| bins | peak t per bin (1 -> 6) |
|---|---|
| null, exact angle | -.8 -.6 -.2 +.2 +.6 +.8 (drifts) |
| null, d' 4 | -.8 -.2 -.2 +.2 +.2 +.6 |
| null, d' 1.1 | -.2 -.2 -.2 .0 .0 +.2 (barely moves) |
| fix, d' 1.1 | +.8 +.8 +.8 -.8 -.8 -.8 (jumps between two fixed positions) |

A peak that jumps between two fixed positions is a visible discreteness signature even at d' 1.1, but it is path-based and inherits the TUBE fragility (section 2).

Quality / pose controls (bins by particle B, amplitude, input-pose error; fixture and null; both sites):

- Amplitude bins give a strong raw first factor (z1 11.6-14.7) that is pure scale. The per-half nuisance projection removes it (FAC nuis z1 -1.1 to -0.3).
- B bins (raw z1 <= 0.4) and pose-error bins (raw <= 2.3, nuis <= 0.8) show none.
- All are gated out (no call) by FAC nuis z1 >= 3.

Particles per bin:

- Two states vs continuum by rank / AFF: the null's z2 is already 12-13 at 6.7k per bin with exact orderings. The limit is the ordering, not N.
- GAPw: the bins only place the lobes (z1 >= 20 at 6.7-8.3k per bin). The gap is read on the consensus.
  - Oracle: .05 with 25k particles in the whole set; .09 at 25k with a d' 1.1 ordering.
  - Input-free (cross-validated ordering): .20 / .23 / .32 at 100k / 50k / 25k. It needs >= ~50k particles at the node (>= 8k per bin, K 6) for this 50 kDa arm.
  - The null side needs the G1 gate: at 25k with a d' 1.1 ordering the null factor is weak (r1 .18) and reads .21.

## 2. Rejected candidates (why)

- **Affine residual, factor count, peak positions**: no power at realistic d' (table, rows d' 2 / 1.1). With d' >= 4, rank 1 (fixture) vs rank >= 2 (null) is clean, and rank 2 (3 states, z3 ~0) vs >= 3 (null, z3 8-14) also separates. Realistic orderings never get there.
- **TUBE dip** (1-D consensus profile along the first-factor path, cross-validated dip between the lobes):
  - oracle paths: fixture .11-.18 (z 2.0-5.4), null <= .04 (z <= 1.6);
  - stress test (`out/tilt.txt`): tilting the path 10 deg about 12 random axes gives a false null dip (> .07 and z >= 3) in 1-4 of 12 directions, at 20 deg up to 5 of 12. The fixture's dip survives 10-20 deg and is gone at 30 deg;
  - searching the path in the site dilated by 20 A gave false null a2 dips of .30-.37 (z 7-16), present even in the noise-free arm-only maps: a tilted path crosses the arm's own inner / outer-lobe waist;
  - input-free paths differ by 12-18 deg between halves; the input-free fixture's dip vanished (.015);
  - on noise-defined paths (quality orderings, off-site controls) the tube "dips" .3-.9 with z 9-76, so it is meaningless without the factor gates.
- **Slabs over the site's hull**: robust to tilt but integrate the overlapping inner arm parts, so the fixture's dip dilutes to .00-.09.
- **GAPx** (lobes from one half, density read on the other): at weak orderings the single-half lobes are noise-shaped. The between-region then reaches into solvent and the null d' 1.1 drops to .22 (a1) / .27 (a2). Use the cross-half lobes (GAPw).

## 3. Proposed statistic: GAPw (input-free)

Inputs: the node's particles at the site, the site mask, half-map reconstructions.

1. **Ordering (cross-validated, one pass)** (`dc_xorder.py`):
   - per half, a site direction d_h = M (R1 - R2) from that half's own two-group split at the site (xstageb.py `indep`: per-half seeds, own-half references);
   - each half's particles are scored against the other half's direction: s_i = Re sum_f conj(r_i) p_i / n(f) / sqrt(sum_f |p_i|^2 / n(f)), where r_i = the residual to the consensus projection and p_i = CTF x projection of d;
   - bins: K = 6 quantiles of s within each of S = 3 strata of the particle's information along d (a_i^2 sum |p_i|^2 / n(f)), pooled, per half. Every bin then has the same quality / defocus / view mix: B and amplitude eta2 0.000 on both datasets.
2. **Gate G1 (signal-dominated first factor)**: the first factor map's cross-half correlation r1 >= 0.5 (signal power >= noise power in each half's factor map; z1 >= ~10 at n_eff ~370), with FAC nuis z1 >= 3 as the heterogeneity test.
   - Why not z1 >= 3 alone: the oracle null at 25k particles with a d' 1.1 ordering has z1 3.5 (r .18). Its noise-shaped lobes passed G2 / G3 (L 13 A, 26 deg) and gave GAPw .21, a false DISCRETE.
   - r1 >= .5 is the half-map "signal = noise" point; it was not tuned on the truth.
3. **Lobes**: w = sign(F1e + F1o) sqrt(max(F1e F1o, 0)), where F1h is the first cross-half factor applied to half h's bins. Selection imprints ("Einstein from noise" along each half's scoring direction) differ between halves and cancel in w. M+ = w > 0.5 max w, M- = w < 0.5 min w, in the site.
4. **Gate G2 (a dipole: density moves)**: half-separation of the lobe centroids L >= 10 A.
5. **Gate G3 (lobes reproducible)**: the even-half and odd-half lobe axes agree within 30 deg. Randomly oriented axes average 57 deg; real factors here read <= 23 deg, and the one marginal factor (z1 3.1) read 51 deg.
6. **GAPw**: H = convex hull of M+ and M-; B = H minus the lobes dilated by 2 voxels (10 A); require |B| >= 10 voxels, else no call. GAPw = mean consensus density in B / the lower of the two lobes' mean densities, per half, averaged.
7. **Call DISCRETE (two states)** if G1-G3 pass and GAPw < 0.28.

Threshold from the null, not the truth:

- The 26 null bin sets passing G1-G3 (a1 10, a2 10, a3 3, a4 3; exact to d' 1.1 orderings) read .37-.83 at tau .5: median .455, MAD-based sd .059.
- .28 = median - 3 robust sd (= .75 x the null minimum .37). With a1 / a2 alone it was .30.
- The input-free null reads 1.14-1.26 or has B < 10 (no call). The fixture values were not used.
- r1 = tanh(z1 / sqrt(n_eff - 3)). The weakest passing null sets are the d' 1.1 ones (z1 10.2-11.5, r1 .51-.55, GAPw .42-.83).
- Constants: tau = .5, 10 A dilation, K = 6, S = 3, r1 >= .5, z1 >= 3 (as discover_sites), L >= 10 A, 30 deg, |B| >= 10, threshold .28.

## 4. Input-free values

| site / ordering | z1 | L (A) | lobe axis even/odd | GAPw | call |
|---|---|---|---|---|---|
| fix a1, cross-validated (xsb indep directions), 100k | 26.9 | 16 | 12 deg | .20 | DISCRETE |
| fix a1, same, 50k (8.3k per bin) | 26.3 | 16 | 10 deg | .23 | DISCRETE |
| fix a1, same, 25k (4.2k per bin) | 24.6 | 16 | 12 deg | .32 | not discrete (miss) |
| null a1, cross-validated, 100k / 50k / 25k | 17.1 / 16.9 / 16.2 | 14-15 | 14-18 deg | 1.22 (B 11) / 1.14 (B 3) / 1.26 (B 9) | not discrete / no call / no call |
| oracle fix d' 1.1, 25k (4.2k per bin) | 15.5 | 19 | 10 deg | .09 | DISCRETE |
| oracle null d' 1.1, 25k | 3.5 (r1 .18) | 13 | 26 deg | .21 | gated by G1 (r1); a false DISCRETE under z1 >= 3 alone |
| fix a1, merged-reference evidence bins (fix_m_a1, K 4) | 49.4 | 16 | 3 deg | .02 | DISCRETE (but z2 4.9 / AFF 5.0) |
| null a1, merged-reference evidence bins (null_m_a1) | 38.9 (noise-locked) | 2 | 23 deg | .41 | gated (G2) |
| null, a1 ordering read at the a2 site | 3.1 | 10 | 51 deg | .20 | gated (G1, G3) |
| fix / null a1, data-derived quality orderings: bins by information along d | raw 15.0 / 16.6 (scale), nuis 1.6 / 1.2 | - | - | - | gated (G1) |
| fix / null a1, bins by abs(score) (confidence) | raw 1.6 / 1.3, nuis 1.0 / 1.2 | - | - | - | gated (G1) |
| fix a2, cross-validated (own per-half splits at site_02) | 13.6 (r1 .64) | 16 | 8 deg | 1.34 (B 5) | no call (B < 10): the three-state miss |
| null a2, cross-validated (own per-half splits at site_04) | -1.1 | 7 | 144 deg | - | gated (G1): no reproducible site factor |

- What the input-free orderings track (truth only scores):
  - fixture: bins' state-2 fraction .82 -> .18, per-half eta2 .17 / .20, score-state corr .42-.46;
  - null: bin mean angle +5.0 -> -5.6 deg, eta2 .05-.06;
  - the per-half xsb indep splits themselves: fixture NMI .06-.09 per half; null angle eta2 .026-.032 per half (real, weak); the halves' labels are anti-aligned (cross-half r -.56 / -.39, flipped);
  - a2 per-half splits (`xsplit/`: xstageb.py indep, 6 iterations; the a1 prep reused, since the consensus and residuals are site-independent):
    - fixture: NMI .026 / .009 per half, cross-half r -.47; bins' 0-deg-state fraction .23 -> .46, eta2 .03-.07;
    - null: cross-half r -.05, per-half directions r .01; bin mean angle -2.4 -> +2.4 deg, eta2 .00 / .04;
  - quality eta2 of the bins 0.000 throughout.
- The input-free fixture's GAPw (.20) is worse than the oracle's (.05): residual shared imprint along the masked direction pulls the lobes inward (L 16 vs 19 A, B 24 vs ~150 voxels).

## 5. Known failure modes

- **Three or more states in a row** (a2): the middle state fills the hull, so it reads not discrete. Rank 2 vs >= 3 (z3) needs exact orderings.
- **Weak factors and small nodes**:
  - r1 < .5 means no call; noise-shaped lobes otherwise give false gaps (null at 25k: .21);
  - the input-free fixture needs >= ~50k particles at the node: at 25k it reads .32, a miss.
- **Close states**: lobes nearer than ~1 resolution element plus the lobe size leave B empty (no call). The fixture's a1 outer lobes are ~78 A apart; real IgM Fabs may be closer or more flexible.
- **Noise-locked orderings** (merged-reference stage B) give reproducible fake factors. G2 / G3 caught this case only because the locked pattern had no dipole. Use cross-validated orderings only.
- **tau sensitivity**:
  - at tau .4 every weak-ordering null drops (d' 1.1: a1 .29, a2 .30, a3 .12, a4 .18);
  - at tau .6 the input-free fixture rises to .38-.40 (null >= .36), so the call flips;
  - tau = .5 is the only setting with margin both ways (input-free .20-.23, oracle 2-state arms <= .15, vs null >= .37).
- **Calibration is this null's**: GAPw of a continuum (.4-.55 here) is a property of the fan's geometry (the noise-free arm-only maps give the same values). A real node needs its own null level; that is not designed yet.
- **Bimodal continua** (a hinge mostly at its two ends) read discrete, which is arguably right.
- **Population balance and separation**: .65/.35 and .70/.30 two-state arms at 40 deg (a3, a4) are detected at the oracle level (.13-.15). Not tested: minor states below ~30 %, separations below 40 deg, and a3 / a4 with input-free orderings.
- **Blur**: quality and pose blur fill gaps, making the call conservative (fewer discrete calls), never false.

## 6. GAPw on the cross-validated stage-B battery (xstageb v2, 2026-11-06)

Inputs: the coordinating session's runs under `xsb/`:

- swap design, half 1 renumbered to half 0 every step, 10 iterations;
- fixture and null at a1-a5, plus a1 seed 2;
- the rigid-core patch controls `null_x_core` (swap) and `null_i_core` (indep).

Constants unchanged: tau .5, 10 A lobe dilation, G1 z1 >= 3 and r1 >= .5, G2 L >= 10 A, G3 axis <= 30 deg, |B| >= 10, DISCRETE if GAPw < .28.

The site is taken from each run's SITEMASK. The crop is the mask dilated by 30 A (`scripts/dc_battery.py`), so a5 and the core patch need no truth crops. Two readouts per run:

- **six-bin**: the validated call (section 3): a one-pass cross-validated ordering scored against the run's own per-half maps (`dc_xorder2.py`, K 6, S 3), reconstructed and run through GAPw.
- **direct**: the run's final two groups (`labels_final.txt`, `rec/`) used as K = 2 bins.

The code reproduces the earlier input-free a1 values (.20 / 1.22).

| run | what the split tracks (truth, from the run's log) | six-bin: z1 / r1, L, axis, B, GAPw -> call | direct: GAPw -> call |
|---|---|---|---|
| fix_x_a1 | a1 NMI .146 | 41.4 / .97, 16 A, 2 deg, 30, .10 -> **DISCRETE** | .10 -> **DISCRETE** |
| fix_x_a1_s2 | a1 NMI .138 | 37.5 / .96, 16 A, 5 deg, 49, .09 -> **DISCRETE** | .08 -> **DISCRETE** |
| fix_x_a2 (3 states) | a2 NMI .047 | 27.1 / .91, 16 A, 8 deg, B 3 -> no call (B) | B 1 -> no call (B) |
| fix_x_a3 | a3 NMI .014 | 4.6 / .24, axis 72 deg -> no call (G1, G3) | no call (G1, G3) |
| fix_x_a4 | a4 NMI .007 (tracks nothing of a4) | 15.7 / .69, 13 A, 18 deg, B 5 -> no call (B) | B 7 -> no call (B) |
| fix_x_a5 (3 states) | a5 NMI .039 | 17.7 / .76, 14 A, 9 deg, B 0 -> no call (B) | B 1 -> no call (B) |
| null_x_a1 | a1 angle eta2 .069 | 28.2 / .91, 14 A, 12 deg, B 1 -> no call (B) | B 1 -> no call (B) |
| null_x_a1_s2 | a1 eta2 .061 | 23.5 / .86, 13 A, 9 deg, B 1 -> no call (B) | B 1 -> no call (B) |
| null_x_a2 | a2 eta2 .021 | 3.2 / .17 -> no call (G1, G3, B) | no call (G1, G3) |
| null_x_a3 | a3 eta2 .038 | 11.5 / .55, axis 50 deg, B 3 -> no call (G3, B) | no call (G1, G3) |
| null_x_a4 | a4 eta2 .000 (tracks nothing) | 15.6 / .69, 11 A, 22 deg, B 4 -> no call (B) | B 4 -> no call (B) |
| null_x_a5 | a5 eta2 .062 | 23.1 / .86, 13 A, 7 deg, B 0 -> no call (B) | B 0 -> no call (B) |
| null_x_core (swap) | no arm (eta2 0); split cross-half r .899 | 26.5 / .90 (z2 6.1), 15 A, 10 deg, B 111, 1.45 -> **not discrete** | 1.44 -> **not discrete** |
| null_i_core (indep) | no arm; r .889 | 26.6 / .90 (z2 5.9), 15 A, 10 deg, B 98, 1.54 -> **not discrete** | 1.54 -> **not discrete** |

- **Null side: no false DISCRETE** on any of the 9 null / control runs in either readout. The fixture a1 split is DISCRETE in both seeds and both readouts. It is cleaner than the earlier indep-based input-free value (.20): the swap split is stronger (NMI .146 vs .06-.09), so the bins' a1 eta2 is .25 (vs .17-.20).
- **Core controls: G1 does not gate them.** The swap and indep stage B find a strongly reproducible split in the rigid patch (cross-half r .89-.90, FAC nuis z1 ~26, r1 .90; P/N site / rest .84-.90, so not site-specific). It is not scale or blur: it survives the per-half nuisance projection (z2 ~6). GAPw stops it: the patch has density between the lobes (1.45-1.54, B ~100), so not discrete.
- **"No call (B)" dominates.** On most arm sites the two lobes leave fewer than 10 voxels between them, for the null's continuous splits (a1, a5) and the fixture's 3-state arms (a2, a5) alike. That is the safe outcome (no split), but GAPw then cannot tell them apart.
- **Fixture a3 / a4 are not split by the battery:**
  - a3 NMI .014, G1 fails;
  - a4 NMI .007, and the split tracks something outside the site (P/N relative .63);
  - both were DISCRETE at the oracle level (.13-.15), so the limit there is the stage B, not GAPw.
- **Caveat (swap design):** the per-half maps used as six-bin directions were themselves labelled against the other half's maps, so the two directions are not independent. Any lock of the swap loop would be inherited. The core controls show that such reproducible splits exist (r .9), and GAPw, not G1, is what stops them here.

Outputs: `discrete/out/battery_<run>.txt` (truth diagnostics of the six-bin ordering included); bins in `rec/<ds>_b_<run>/`; jobs 54979475-54979709 (8 cores, 16 GB, ~8 min each).

## 7. GAPw on the gate-replay runs (RNase P, null fixtures; xsbg.py, 2026-11-06)

The call is now dataset-generic (`scripts/dc_gates.py`, `dc_gorder.py`, `dc_gates.sbatch`). Conversions, all physical constants unchanged:

- **Analysis grid:** every map and mask is Fourier-cropped to the even box nearest to 5.11 A/px, the IgM fixture's validated grid. Fixture 200 x 3.27 A -> 128 x 5.11 A; RNase P and qtoynull 160 x 1.464 A -> 46 x 5.09 A; olenull 180 x 1.908 A -> 68 x 5.05 A.
- **Re-derived on each grid:**
  - low-pass 15 A (Gaussian, half-max at 15 A);
  - crop = the site mask dilated by 30 A (6 voxels);
  - lobe dilation 10 A = round(10 / smpd) = 2 voxels;
  - |B| >= 1,330 A^3 = round(1330 / smpd^3) = 10 voxels;
  - n_eff = mask voxels x (4/3) pi (smpd / 15)^3 / 2;
  - L >= 10 A measured in A.
- **Dimensionless, unchanged:** tau .5, axes <= 30 deg, r1 >= .5, z1 >= 3, threshold .28.
- **Ordering:** in xsbg.py's own prep conventions:
  - the projector Proj(v, vsmpd) imported from xsbg.py;
  - |CTF| for phase-flipped stacks (sg 0);
  - a_i as xsbg fitted it, so no projector scale c.
- **Six-bin reconstruction:** on a copy of the node's stageB/b.simple with the run's mskdiam / pgrp (c2 for olenull), at the native box.
- **Checks:**
  - the battery runs fix_x_a1 / null_x_core reproduce exactly (.10 / .10 / 1.44);
  - the gate-path validation v_fix_a1 is DISCRETE: direct .09 (z1 38.5, r1 .96, L 16 A, axis 3 deg, B 40); six-bin .09 (z1 41.4, bins' a1 state-2 fraction .86 -> .13, other arms flat). This reproduces fix_x_a1 (.08-.10).

RNase P runs also carry, from the run's own log (logs/igmfx_gxsb.<job>.out): SEPARATES(ref), the closed reference rows' recall / purity per group, and the arm readout of its two group maps.

Totals over the 56 gate runs (the coordinating session's xsbg.py swap stage B, 10 iterations, at every site of each saved node; island masks as i<NN>):

- **Null fixtures** (qtoynull_r 20 sites, olenull_r 4): **0 DISCRETE**, so no false positive. Every run fails G1 (z1 <= 3.4, r1 <= .57); the one six-bin run with r1 >= .5 (qtoynull_r_s07, r1 .52) fails G2 / G3 / B (lobes 4 A apart).
- **RNase P** (r131 s2 / s1 nodes r, r1, r2; the 692k closed leaves lf_r1312, lf_r1321, lf_r2112; 32 runs): **0 DISCRETE**. No call either way, because no split carries a closed state:
  - SEPARATES(ref) is no on every run, and the closed reference rows split evenly (recall .32-.68 per group, purity at the node share);
  - the arm readout calls both group maps intermediate on every run.
  - So the per-site stage B produced no closed-state split for GAPw to accept or lose.
- **Factors that pass G1 but are no two-lobe motion:** r131s1_r_s01 (r1 .61 / .64), r131s1_r2_s04 (.65 / .71), lf_r2112_i02 (.73 / .80). G2 / B stop them: lobes 7-8 A apart, no between-region.
- **GAPw "n/a"** marks noise lobes lying on solvent (lobe density <= 0 makes the ratio meaningless); G1 had already failed in every such run.
- **Combined path:** with these runs, the add-only per-site layer (islands + cross-validated stage B + split only on DISCRETE) splits nothing: it holds both nulls and keeps every RNase P node as it was.

| run | node, site | direct: z1 r1, L, axis, B, GAPw -> call | six-bin: same -> call | RNase P (run log): closed reference and arm readout |
|---|---|---|---|---|
| v_fix_a1 | igmfix_root site_06 | z1 38.5 r1 .96, L 16 A, axis 3, B 40, GAPw .09 -> **DISCRETE** | z1 41.4 r1 .97, L 16 A, axis 2, B 40, GAPw .09 -> **DISCRETE** | - |
| lf_r1312_i01 | lf_r1312 site_01 | z1 1.6 r1 0.19, L 3 A, axis 108, B 4, GAPw 1.24 -> no call (G1+G2+G3+B) | z1 1.8 r1 0.20, L 4 A, axis 99, B 3, GAPw 1.01 -> no call (G1+G2+G3+B) | SEPARATES no; recall/purity g1 0.54/0.0781 g2 0.46/0.0753 ; arm g1 intermediate g2 intermediate |
| lf_r1312_i02 | lf_r1312 site_02 | z1 0.8 r1 0.09, L 5 A, axis 65, B 13, GAPw -0.45 -> no call (G1+G2+G3) | z1 0.8 r1 0.10, L 2 A, axis 150, B 15, GAPw -0.26 -> no call (G1+G2+G3) | SEPARATES no; recall/purity g1 0.49/0.0762 g2 0.51/0.0773 ; arm g1 intermediate g2 intermediate |
| lf_r1312_i03 | lf_r1312 site_03 | z1 2.0 r1 0.26, L 2 A, axis 88, B 12, GAPw n/a -> no call (G1+G2+G3) | z1 2.4 r1 0.32, L 2 A, axis 67, B 1, GAPw nan -> no call (G1+G2+G3+B) | SEPARATES no; recall/purity g1 0.47/0.0716 g2 0.53/0.0820 ; arm g1 intermediate g2 intermediate |
| lf_r1321_i01 | lf_r1321 site_01 | z1 0.8 r1 0.09, L 9 A, axis 42, B 51, GAPw 1.31 -> no call (G1+G2+G3) | z1 0.8 r1 0.09, L 8 A, axis 102, B 36, GAPw 1.42 -> no call (G1+G2+G3) | SEPARATES no; recall/purity g1 0.49/0.0616 g2 0.51/0.0611 ; arm g1 intermediate g2 intermediate |
| lf_r1321_i02 | lf_r1321 site_02 | z1 0.1 r1 0.02, L 7 A, axis 59, B 0, GAPw nan -> no call (G1+G2+G3+B) | z1 0.2 r1 0.03, L 7 A, axis 40, B 0, GAPw nan -> no call (G1+G2+G3+B) | SEPARATES no; recall/purity g1 0.45/0.0541 g2 0.55/0.0689 ; arm g1 intermediate g2 intermediate |
| lf_r1321_i03 | lf_r1321 site_03 | z1 0.0 r1 0.01, L 3 A, axis 115, B 0, GAPw nan -> no call (G1+G2+G3+B) | z1 0.1 r1 0.02, L 4 A, axis 98, B 0, GAPw nan -> no call (G1+G2+G3+B) | SEPARATES no; recall/purity g1 0.51/0.0625 g2 0.49/0.0602 ; arm g1 intermediate g2 intermediate |
| lf_r1321_i04 | lf_r1321 site_04 | z1 0.3 r1 0.05, L 7 A, axis 43, B 0, GAPw nan -> no call (G1+G2+G3+B) | z1 -0.1 r1 -0.02, lobes < 4 vox -> no call (G1) | SEPARATES no; recall/purity g1 0.49/0.0604 g2 0.51/0.0623 ; arm g1 intermediate g2 intermediate |
| lf_r2112_i01 | lf_r2112 site_01 | z1 1.1 r1 0.12, L 13 A, axis 7, B 14, GAPw 1.81 -> no call (G1) | z1 1.2 r1 0.13, L 14 A, axis 9, B 11, GAPw 2.00 -> no call (G1) | SEPARATES no; recall/purity g1 0.43/0.0846 g2 0.57/0.1056 ; arm g1 intermediate g2 intermediate |
| lf_r2112_i02 | lf_r2112 site_02 | z1 5.3 r1 0.73, L 7 A, axis 27, B 0, GAPw nan -> no call (G2+B) | z1 6.3 r1 0.80, L 7 A, axis 24, B 0, GAPw nan -> no call (G2+B) | SEPARATES no; recall/purity g1 0.44/0.0859 g2 0.56/0.1043 ; arm g1 intermediate g2 intermediate |
| lf_r2112_i03 | lf_r2112 site_03 | z1 0.6 r1 0.10, L 6 A, axis 58, B 1, GAPw nan -> no call (G1+G2+G3+B) | z1 0.8 r1 0.13, L 7 A, axis 59, B 0, GAPw nan -> no call (G1+G2+G3+B) | SEPARATES no; recall/purity g1 0.48/0.0926 g2 0.52/0.0981 ; arm g1 intermediate g2 intermediate |
| lf_r2112_i04 | lf_r2112 site_04 | z1 2.8 r1 0.39, L 12 A, axis 26, B 3, GAPw 1.91 -> no call (G1+B) | z1 3.0 r1 0.41, L 12 A, axis 28, B 2, GAPw nan -> no call (G1+B) | SEPARATES no; recall/purity g1 0.53/0.1018 g2 0.47/0.0890 ; arm g1 intermediate g2 intermediate |
| lf_r2112_i05 | lf_r2112 site_05 | z1 0.3 r1 0.05, L 10 A, axis 20, B 1, GAPw nan -> no call (G1+B) | z1 0.4 r1 0.07, L 8 A, axis 34, B 0, GAPw nan -> no call (G1+G2+G3+B) | SEPARATES no; recall/purity g1 0.46/0.0896 g2 0.54/0.1009 ; arm g1 intermediate g2 intermediate |
| olenull_r_s01 | olenull_r site_01 | z1 1.4 r1 0.12, L 1 A, axis 180, B 10, GAPw 1.18 -> no call (G1+G2+G3) | z1 1.5 r1 0.12, L 1 A, axis 180, B 10, GAPw 1.31 -> no call (G1+G2+G3) | - |
| olenull_r_s02 | olenull_r site_02 | z1 2.4 r1 0.20, L 8 A, axis 41, B 69, GAPw 6.02 -> no call (G1+G2+G3) | z1 2.6 r1 0.21, L 9 A, axis 17, B 62, GAPw 7.22 -> no call (G1+G2) | - |
| olenull_r_s03 | olenull_r site_03 | z1 0.7 r1 0.09, L 9 A, axis 62, B 3, GAPw n/a -> no call (G1+G2+G3+B) | z1 0.9 r1 0.11, L 9 A, axis 57, B 5, GAPw n/a -> no call (G1+G2+G3+B) | - |
| olenull_r_s04 | olenull_r site_04 | z1 1.3 r1 0.20, L 2 A, axis 47, B 0, GAPw nan -> no call (G1+G2+G3+B) | z1 1.3 r1 0.20, L 3 A, axis 45, B 1, GAPw nan -> no call (G1+G2+G3+B) | - |
| qtoynull_r_s01 | qtoynull_r site_01 | z1 1.9 r1 0.14, L 8 A, axis 80, B 62, GAPw n/a -> no call (G1+G2+G3) | z1 2.6 r1 0.19, L 11 A, axis 108, B 53, GAPw n/a -> no call (G1+G3) | - |
| qtoynull_r_s02 | qtoynull_r site_02 | z1 1.3 r1 0.08, L 5 A, axis 111, B 97, GAPw n/a -> no call (G1+G2+G3) | z1 -1.4 r1 -0.09, L 3 A, axis 119, B 625, GAPw n/a -> no call (G1+G2+G3) | - |
| qtoynull_r_s03 | qtoynull_r site_03 | z1 1.5 r1 0.12, L 12 A, axis 72, B 39, GAPw n/a -> no call (G1+G3) | z1 2.0 r1 0.16, L 12 A, axis 67, B 69, GAPw n/a -> no call (G1+G3) | - |
| qtoynull_r_s04 | qtoynull_r site_04 | z1 -0.1 r1 -0.02, L 11 A, axis 34, B 14, GAPw n/a -> no call (G1+G3) | z1 0.1 r1 0.01, L 10 A, axis 107, B 56, GAPw n/a -> no call (G1+G2+G3) | - |
| qtoynull_r_s05 | qtoynull_r site_05 | z1 1.4 r1 0.18, L 9 A, axis 25, B 2, GAPw nan -> no call (G1+G2+B) | z1 1.6 r1 0.21, L 8 A, axis 23, B 2, GAPw nan -> no call (G1+G2+B) | - |
| qtoynull_r_s06 | qtoynull_r site_06 | z1 1.2 r1 0.22, L 4 A, axis 150, B 0, GAPw nan -> no call (G1+G2+G3+B) | z1 1.2 r1 0.21, L 4 A, axis 140, B 0, GAPw nan -> no call (G1+G2+G3+B) | - |
| qtoynull_r_s07 | qtoynull_r site_07 | z1 3.1 r1 0.48, L 4 A, axis 77, B 0, GAPw nan -> no call (G1+G2+G3+B) | z1 3.4 r1 0.52, L 4 A, axis 72, B 0, GAPw nan -> no call (G2+G3+B) | - |
| qtoynull_r_s08 | qtoynull_r site_08 | z1 1.7 r1 0.32, lobes < 4 vox -> no call (G1) | z1 1.8 r1 0.35, lobes < 4 vox -> no call (G1) | - |
| qtoynull_r_s09 | qtoynull_r site_09 | z1 2.2 r1 0.40, L 6 A, axis 70, B 0, GAPw nan -> no call (G1+G2+G3+B) | z1 2.3 r1 0.40, L 7 A, axis 56, B 0, GAPw nan -> no call (G1+G2+G3+B) | - |
| qtoynull_r_s10 | qtoynull_r site_10 | z1 0.2 r1 0.04, L 5 A, axis 78, B 0, GAPw nan -> no call (G1+G2+G3+B) | z1 0.1 r1 0.01, L 5 A, axis 74, B 0, GAPw nan -> no call (G1+G2+G3+B) | - |
| qtoynull_r_s11 | qtoynull_r site_11 | z1 0.4 r1 0.08, L 5 A, axis 81, B 0, GAPw nan -> no call (G1+G2+G3+B) | z1 0.5 r1 0.11, L 6 A, axis 80, B 0, GAPw nan -> no call (G1+G2+G3+B) | - |
| qtoynull_r_s12 | qtoynull_r site_12 | z1 1.4 r1 0.30, L 7 A, axis 40, B 1, GAPw nan -> no call (G1+G2+G3+B) | z1 1.5 r1 0.32, lobes < 4 vox -> no call (G1) | - |
| qtoynull_r_s13 | qtoynull_r site_13 | z1 0.5 r1 0.11, L 6 A, axis 152, B 0, GAPw nan -> no call (G1+G2+G3+B) | z1 0.5 r1 0.10, lobes < 4 vox -> no call (G1) | - |
| qtoynull_r_s14 | qtoynull_r site_14 | z1 1.0 r1 0.22, L 6 A, axis 43, B 0, GAPw nan -> no call (G1+G2+G3+B) | z1 1.0 r1 0.22, L 6 A, axis 48, B 0, GAPw nan -> no call (G1+G2+G3+B) | - |
| qtoynull_r_s15 | qtoynull_r site_15 | z1 1.7 r1 0.31, L 6 A, axis 79, B 0, GAPw nan -> no call (G1+G2+G3+B) | z1 1.8 r1 0.33, L 6 A, axis 78, B 0, GAPw nan -> no call (G1+G2+G3+B) | - |
| qtoynull_r_s16 | qtoynull_r site_16 | z1 1.8 r1 0.40, L 7 A, axis 56, B 0, GAPw nan -> no call (G1+G2+G3+B) | z1 1.9 r1 0.41, L 7 A, axis 60, B 0, GAPw nan -> no call (G1+G2+G3+B) | - |
| qtoynull_r_s17 | qtoynull_r site_17 | z1 1.9 r1 0.40, lobes < 4 vox -> no call (G1) | z1 2.0 r1 0.41, L 9 A, axis 42, B 0, GAPw nan -> no call (G1+G2+G3+B) | - |
| qtoynull_r_s18 | qtoynull_r site_18 | z1 0.8 r1 0.19, L 7 A, axis 78, B 0, GAPw nan -> no call (G1+G2+G3+B) | z1 0.9 r1 0.20, L 7 A, axis 77, B 0, GAPw nan -> no call (G1+G2+G3+B) | - |
| qtoynull_r_s19 | qtoynull_r site_19 | z1 2.5 r1 0.56, lobes < 4 vox -> no call (G1) | z1 2.6 r1 0.57, lobes < 4 vox -> no call (G1) | - |
| qtoynull_r_s20 | qtoynull_r site_20 | z1 0.5 r1 0.13, lobes < 4 vox -> no call (G1) | z1 0.7 r1 0.19, lobes < 4 vox -> no call (G1) | - |
| r131s1_r1_i01 | r131s1_r1 site_01 | z1 2.3 r1 0.28, L 9 A, axis 39, B 0, GAPw nan -> no call (G1+G2+G3+B) | z1 2.8 r1 0.33, L 9 A, axis 22, B 0, GAPw nan -> no call (G1+G2+B) | SEPARATES no; recall/purity g1 0.60/0.0819 g2 0.40/0.0533 ; arm g1 intermediate g2 intermediate |
| r131s1_r1_s01 | r131s1_r1 site_01 | z1 0.6 r1 0.09, L 7 A, axis 65, B 1, GAPw nan -> no call (G1+G2+G3+B) | z1 0.6 r1 0.09, L 8 A, axis 66, B 2, GAPw nan -> no call (G1+G2+G3+B) | SEPARATES no; recall/purity g1 0.57/0.0743 g2 0.43/0.0601 ; arm g1 intermediate g2 intermediate |
| r131s1_r1_s02 | r131s1_r1 site_02 | z1 0.2 r1 0.03, L 8 A, axis 63, B 1, GAPw nan -> no call (G1+G2+G3+B) | z1 -0.3 r1 -0.06, L 5 A, axis 127, B 0, GAPw nan -> no call (G1+G2+G3+B) | SEPARATES no; recall/purity g1 0.54/0.0725 g2 0.46/0.0624 ; arm g1 intermediate g2 intermediate |
| r131s1_r2_s01 | r131s1_r2 site_01 | z1 1.7 r1 0.21, L 6 A, axis 82, B 8, GAPw 2.60 -> no call (G1+G2+G3+B) | z1 1.3 r1 0.17, L 6 A, axis 95, B 10, GAPw 2.42 -> no call (G1+G2+G3) | SEPARATES no; recall/purity g1 0.54/0.0734 g2 0.46/0.0694 ; arm g1 intermediate g2 intermediate |
| r131s1_r2_s02 | r131s1_r2 site_02 | z1 0.2 r1 0.02, L 8 A, axis 129, B 0, GAPw nan -> no call (G1+G2+G3+B) | z1 0.1 r1 0.02, L 8 A, axis 104, B 0, GAPw nan -> no call (G1+G2+G3+B) | SEPARATES no; recall/purity g1 0.56/0.0797 g2 0.44/0.0630 ; arm g1 intermediate g2 intermediate |
| r131s1_r2_s03 | r131s1_r2 site_03 | z1 0.4 r1 0.08, L 5 A, axis 58, B 0, GAPw nan -> no call (G1+G2+G3+B) | z1 0.6 r1 0.11, L 6 A, axis 44, B 0, GAPw nan -> no call (G1+G2+G3+B) | SEPARATES no; recall/purity g1 0.48/0.0693 g2 0.52/0.0736 ; arm g1 intermediate g2 intermediate |
| r131s1_r2_s04 | r131s1_r2 site_04 | z1 4.5 r1 0.65, L 7 A, axis 34, B 0, GAPw nan -> no call (G2+G3+B) | z1 5.0 r1 0.71, L 7 A, axis 24, B 0, GAPw nan -> no call (G2+B) | SEPARATES no; recall/purity g1 0.32/0.0444 g2 0.68/0.0995 ; arm g1 intermediate g2 intermediate |
| r131s1_r_s01 | r131s1_r site_01 | z1 5.8 r1 0.61, L 8 A, axis 29, B 0, GAPw nan -> no call (G2+B) | z1 6.1 r1 0.64, L 8 A, axis 26, B 0, GAPw nan -> no call (G2+B) | SEPARATES no; recall/purity g1 0.49/0.0695 g2 0.51/0.0693 ; arm g1 intermediate g2 intermediate |
| r131s1_r_s02 | r131s1_r site_02 | z1 2.7 r1 0.40, L 12 A, axis 28, B 0, GAPw nan -> no call (G1+B) | z1 3.0 r1 0.45, L 12 A, axis 25, B 0, GAPw nan -> no call (G1+B) | SEPARATES no; recall/purity g1 0.65/0.0834 g2 0.35/0.0527 ; arm g1 intermediate g2 intermediate |
| r131s2_r1_s01 | r131s2_r1 site_01 | z1 0.8 r1 0.09, L 10 A, axis 2, B 1, GAPw nan -> no call (G1+G2+B) | z1 1.2 r1 0.14, L 13 A, axis 5, B 1, GAPw nan -> no call (G1+B) | SEPARATES no; recall/purity g1 0.38/0.0570 g2 0.62/0.0845 ; arm g1 intermediate g2 intermediate |
| r131s2_r1_s02 | r131s2_r1 site_02 | z1 0.9 r1 0.10, L 8 A, axis 9, B 5, GAPw 0.70 -> no call (G1+G2+B) | z1 1.1 r1 0.12, L 9 A, axis 10, B 6, GAPw 1.11 -> no call (G1+G2+B) | SEPARATES no; recall/purity g1 0.50/0.0715 g2 0.50/0.0711 ; arm g1 intermediate g2 intermediate |
| r131s2_r1_s03 | r131s2_r1 site_03 | z1 1.0 r1 0.16, L 7 A, axis 62, B 0, GAPw nan -> no call (G1+G2+G3+B) | z1 0.8 r1 0.12, L 8 A, axis 63, B 0, GAPw nan -> no call (G1+G2+G3+B) | SEPARATES no; recall/purity g1 0.57/0.0862 g2 0.43/0.0578 ; arm g1 intermediate g2 intermediate |
| r131s2_r1_s04 | r131s2_r1 site_04 | z1 0.4 r1 0.07, L 7 A, axis 13, B 0, GAPw nan -> no call (G1+G2+B) | z1 0.5 r1 0.09, L 7 A, axis 17, B 0, GAPw nan -> no call (G1+G2+B) | SEPARATES no; recall/purity g1 0.55/0.0759 g2 0.45/0.0665 ; arm g1 intermediate g2 intermediate |
| r131s2_r1_s05 | r131s2_r1 site_05 | z1 0.3 r1 0.05, L 6 A, axis 71, B 0, GAPw nan -> no call (G1+G2+G3+B) | z1 0.4 r1 0.07, L 6 A, axis 72, B 0, GAPw nan -> no call (G1+G2+G3+B) | SEPARATES no; recall/purity g1 0.49/0.0677 g2 0.51/0.0751 ; arm g1 intermediate g2 intermediate |
| r131s2_r2_s01 | r131s2_r2 site_01 | z1 -0.8 r1 -0.10, L 12 A, axis 62, B 1, GAPw nan -> no call (G1+G3+B) | z1 0.4 r1 0.06, L 6 A, axis 68, B 14, GAPw 0.35 -> no call (G1+G2+G3) | SEPARATES no; recall/purity g1 0.44/0.0594 g2 0.56/0.0740 ; arm g1 intermediate g2 intermediate |
| r131s2_r2_s02 | r131s2_r2 site_02 | z1 1.3 r1 0.22, L 6 A, axis 71, B 0, GAPw nan -> no call (G1+G2+G3+B) | z1 1.7 r1 0.29, L 6 A, axis 67, B 0, GAPw nan -> no call (G1+G2+G3+B) | SEPARATES no; recall/purity g1 0.43/0.0568 g2 0.57/0.0770 ; arm g1 intermediate g2 intermediate |
| r131s2_r_s01 | r131s2_r site_01 | z1 1.8 r1 0.19, L 7 A, axis 65, B 0, GAPw nan -> no call (G1+G2+G3+B) | z1 2.0 r1 0.22, L 7 A, axis 64, B 0, GAPw nan -> no call (G1+G2+G3+B) | SEPARATES no; recall/purity g1 0.49/0.0680 g2 0.51/0.0707 ; arm g1 intermediate g2 intermediate |
| r131s2_r_s02 | r131s2_r site_02 | z1 1.8 r1 0.25, L 13 A, axis 29, B 4, GAPw 0.93 -> no call (G1+B) | z1 2.2 r1 0.30, L 13 A, axis 24, B 4, GAPw 0.83 -> no call (G1+B) | SEPARATES no; recall/purity g1 0.34/0.0518 g2 0.66/0.0844 ; arm g1 intermediate g2 intermediate |
| r131s2_r_s03 | r131s2_r site_03 | z1 2.3 r1 0.29, L 6 A, axis 72, B 0, GAPw nan -> no call (G1+G2+G3+B) | z1 2.5 r1 0.31, L 6 A, axis 72, B 0, GAPw nan -> no call (G1+G2+G3+B) | SEPARATES no; recall/purity g1 0.48/0.0672 g2 0.52/0.0715 ; arm g1 intermediate g2 intermediate |
| r131s2_r_s04 | r131s2_r site_04 | z1 2.6 r1 0.50, L 5 A, axis 66, B 0, GAPw nan -> no call (G1+G2+G3+B) | z1 2.8 r1 0.52, L 6 A, axis 60, B 0, GAPw nan -> no call (G1+G2+G3+B) | SEPARATES no; recall/purity g1 0.54/0.0737 g2 0.46/0.0649 ; arm g1 intermediate g2 intermediate |

## 7b. GAPw on gate-replay wave 2: the 692k closed lineage's other nodes (2026-11-06)

Same call, constants and conversions as section 7: RNase P, box 160 x 1.464 A -> analysis grid 46 x 5.09 A. Six-bin bins are reconstructed on each prep's node.simple. For r13 / r131 / r21 the NODE is the production node (integration/f692_lf_s1) and the SITEMASK the island_gates rerun's site mask, as each manifest states.

Closed reference share per node (run logs / dc_gorder.py):

| node | particles | closed share |
|---|---|---|
| r13 | 47k | .063 |
| r131 | 25k | .068 |
| r1311 | 12.5k | .060 |
| r132 | 22k | .057 |
| r21 | 53k | .029 |
| r211 | 27k | .052 |

The 692k root rate is .0126.

- **0 DISCRETE in 33 runs.**
- **SEPARATES(ref) is no on every run.**
- **Arm readouts:** both groups "intermediate" on all but three runs; lf_r21_s01-03 have one group "open".
- **Real splits that the gates or GAPw stop:**
  - **lf_r211_s01:**
    - the strongest factor of the wave: z1 23-25, r1 .97-.98, lobes 14 A apart, half axes within 2-6 deg, B 10;
    - it passes every gate and reads **not discrete** (GAPw .45 in both readouts);
    - the groups are ~50/50, and 90 % of the node's closed rows fall in g1 (purity .090 vs .011, the node share .052);
    - so the split carries the closed rows as a minority (9 %) of a broad group, the density between its lobes is filled, and both group maps read intermediate.
  - **lf_r21_s01-03:**
    - strong factors (z1 9-17, r1 .82-.93); each splits ~90 % of the node's closed rows into one group (purity .050 vs .006, node share .029);
    - the arm readout calls the other group "open";
    - stopped as **no call** at G2 / G3: lobes 6-8 A apart, even / odd axes 22-74 deg.
  - **lf_r13_s02** (r1 .58-.65), **lf_r21_s06** (.60-.62) and **lf_r21_s08** (.72-.73) are reproducible but fail G2 (lobes 6-7 A apart).
  - Every other run fails G1.
- **Reading for the gate:** no closed state is gained or lost by the add-only layer here.
  - No per-site split isolates the closed rows (SEPARATES no throughout).
  - The splits that do enrich them (lf_r211_s01, lf_r21_s01-03) are an open / not-open contrast with closed as a minority of one side.
  - A DISCRETE call on them would have produced an "intermediate" group, not a closed one.

| run | node, site | direct: z1 r1, L, axis, B, GAPw -> call | six-bin: same -> call | RNase P (run log): closed reference and arm readout |
|---|---|---|---|---|
| lf_r1311_s02 | lf_r1311 site_02 | z1 2.4 r1 0.29, L 6 A, axis 31, B 2, GAPw nan -> no call (G1+G2+G3+B) | z1 2.8 r1 0.33, L 6 A, axis 27, B 4, GAPw 1.60 -> no call (G1+G2+B) | SEPARATES no; recall/purity g1 0.53/0.0628 g2 0.47/0.0571; arm g1 intermediate g2 intermediate |
| lf_r1311_s03 | lf_r1311 site_03 | z1 1.2 r1 0.15, L 3 A, axis 135, B 1, GAPw nan -> no call (G1+G2+G3+B) | z1 1.4 r1 0.18, L 3 A, axis 109, B 1, GAPw nan -> no call (G1+G2+G3+B) | SEPARATES no; recall/purity g1 0.50/0.0593 g2 0.50/0.0607; arm g1 intermediate g2 intermediate |
| lf_r1311_s04 | lf_r1311 site_04 | z1 1.6 r1 0.24, L 4 A, axis 54, B 0, GAPw nan -> no call (G1+G2+G3+B) | z1 1.7 r1 0.26, L 6 A, axis 39, B 0, GAPw nan -> no call (G1+G2+G3+B) | SEPARATES no; recall/purity g1 0.53/0.0631 g2 0.47/0.0568; arm g1 intermediate g2 intermediate |
| lf_r1311_s05 | lf_r1311 site_05 | z1 1.0 r1 0.17, L 3 A, axis 53, B 0, GAPw nan -> no call (G1+G2+G3+B) | z1 1.0 r1 0.17, L 3 A, axis 67, B 0, GAPw nan -> no call (G1+G2+G3+B) | SEPARATES no; recall/purity g1 0.51/0.0595 g2 0.49/0.0605; arm g1 intermediate g2 intermediate |
| lf_r1311_s06 | lf_r1311 site_06 | z1 1.1 r1 0.18, L 5 A, axis 105, B 0, GAPw nan -> no call (G1+G2+G3+B) | z1 1.3 r1 0.21, L 6 A, axis 98, B 1, GAPw nan -> no call (G1+G2+G3+B) | SEPARATES no; recall/purity g1 0.44/0.0528 g2 0.56/0.0671; arm g1 intermediate g2 intermediate |
| lf_r131_s01 | lf_r131 site_01 | z1 3.1 r1 0.35, L 13 A, axis 27, B 5, GAPw 1.24 -> no call (G1+B) | z1 3.8 r1 0.42, L 13 A, axis 28, B 5, GAPw 1.33 -> no call (G1+B) | SEPARATES no; recall/purity g1 0.49/0.0682 g2 0.51/0.0685; arm g1 intermediate g2 intermediate |
| lf_r131_s02 | lf_r131 site_02 | z1 1.4 r1 0.19, L 6 A, axis 55, B 1, GAPw nan -> no call (G1+G2+G3+B) | z1 1.5 r1 0.21, L 6 A, axis 67, B 0, GAPw nan -> no call (G1+G2+G3+B) | SEPARATES no; recall/purity g1 0.53/0.0709 g2 0.47/0.0657; arm g1 intermediate g2 intermediate |
| lf_r131_s03 | lf_r131 site_03 | z1 0.7 r1 0.11, L 5 A, axis 65, B 0, GAPw nan -> no call (G1+G2+G3+B) | z1 0.7 r1 0.11, L 6 A, axis 71, B 0, GAPw nan -> no call (G1+G2+G3+B) | SEPARATES no; recall/purity g1 0.53/0.0714 g2 0.47/0.0653; arm g1 intermediate g2 intermediate |
| lf_r132_s01 | lf_r132 site_01 | z1 3.2 r1 0.34, L 9 A, axis 23, B 15, GAPw 0.32 -> no call (G1+G2) | z1 3.4 r1 0.36, L 8 A, axis 21, B 13, GAPw 0.30 -> no call (G1+G2) | SEPARATES no; recall/purity g1 0.48/0.0564 g2 0.52/0.0575; arm g1 intermediate g2 intermediate |
| lf_r132_s02 | lf_r132 site_02 | z1 3.3 r1 0.37, L 9 A, axis 71, B 3, GAPw 1.99 -> no call (G1+G2+G3+B) | z1 3.6 r1 0.40, L 9 A, axis 63, B 1, GAPw nan -> no call (G1+G2+G3+B) | SEPARATES no; recall/purity g1 0.53/0.0589 g2 0.47/0.0550; arm g1 intermediate g2 intermediate |
| lf_r132_s03 | lf_r132 site_03 | z1 1.2 r1 0.17, L 5 A, axis 93, B 2, GAPw nan -> no call (G1+G2+G3+B) | z1 1.2 r1 0.17, L 5 A, axis 77, B 0, GAPw nan -> no call (G1+G2+G3+B) | SEPARATES no; recall/purity g1 0.48/0.0545 g2 0.52/0.0595; arm g1 intermediate g2 intermediate |
| lf_r132_s04 | lf_r132 site_04 | z1 0.7 r1 0.10, L 3 A, axis 110, B 0, GAPw nan -> no call (G1+G2+G3+B) | z1 0.6 r1 0.08, L 2 A, axis 104, B 1, GAPw nan -> no call (G1+G2+G3+B) | SEPARATES no; recall/purity g1 0.50/0.0567 g2 0.50/0.0572; arm g1 intermediate g2 intermediate |
| lf_r132_s05 | lf_r132 site_05 | z1 3.2 r1 0.42, L 7 A, axis 76, B 2, GAPw nan -> no call (G1+G2+G3+B) | z1 3.4 r1 0.44, L 6 A, axis 83, B 3, GAPw 0.90 -> no call (G1+G2+G3+B) | SEPARATES no; recall/purity g1 0.51/0.0573 g2 0.49/0.0566; arm g1 intermediate g2 intermediate |
| lf_r132_s06 | lf_r132 site_06 | z1 0.3 r1 0.05, L 8 A, axis 12, B 0, GAPw nan -> no call (G1+G2+B) | z1 0.3 r1 0.06, L 8 A, axis 14, B 0, GAPw nan -> no call (G1+G2+B) | SEPARATES no; recall/purity g1 0.48/0.0546 g2 0.52/0.0594; arm g1 intermediate g2 intermediate |
| lf_r132_s07 | lf_r132 site_07 | z1 0.6 r1 0.10, lobes < 4 vox -> no call (G1) | z1 -0.6 r1 -0.11, lobes < 4 vox -> no call (G1) | SEPARATES no; recall/purity g1 0.49/0.0555 g2 0.51/0.0584; arm g1 intermediate g2 intermediate |
| lf_r13_s01 | lf_r13 site_01 | z1 2.9 r1 0.28, L 12 A, axis 46, B 50, GAPw 1.00 -> no call (G1+G3) | z1 3.6 r1 0.35, L 14 A, axis 44, B 44, GAPw 1.04 -> no call (G1+G3) | SEPARATES no; recall/purity g1 0.54/0.0683 g2 0.46/0.0577; arm g1 intermediate g2 intermediate |
| lf_r13_s02 | lf_r13 site_02 | z1 6.9 r1 0.58, L 6 A, axis 55, B 2, GAPw nan -> no call (G2+G3+B) | z1 8.1 r1 0.65, L 7 A, axis 30, B 4, GAPw -0.10 -> no call (G2+G3+B) | SEPARATES no; recall/purity g1 0.54/0.0674 g2 0.46/0.0585; arm g1 intermediate g2 intermediate |
| lf_r13_s03 | lf_r13 site_03 | z1 0.5 r1 0.06, L 2 A, axis 54, B 1, GAPw nan -> no call (G1+G2+G3+B) | z1 0.5 r1 0.06, L 6 A, axis 40, B 1, GAPw nan -> no call (G1+G2+G3+B) | SEPARATES no; recall/purity g1 0.49/0.0626 g2 0.51/0.0634; arm g1 intermediate g2 intermediate |
| lf_r13_s04 | lf_r13 site_04 | z1 1.0 r1 0.11, L 5 A, axis 102, B 10, GAPw 1.69 -> no call (G1+G2+G3) | z1 1.1 r1 0.13, L 5 A, axis 114, B 10, GAPw 3.02 -> no call (G1+G2+G3) | SEPARATES no; recall/purity g1 0.51/0.0643 g2 0.49/0.0617; arm g1 intermediate g2 intermediate |
| lf_r13_s05 | lf_r13 site_05 | z1 0.8 r1 0.15, L 6 A, axis 30, B 0, GAPw nan -> no call (G1+G2+B) | z1 1.1 r1 0.21, L 6 A, axis 29, B 0, GAPw nan -> no call (G1+G2+B) | SEPARATES no; recall/purity g1 0.42/0.0531 g2 0.58/0.0728; arm g1 intermediate g2 intermediate |
| lf_r211_i02 | lf_r211 site_02 | z1 4.3 r1 0.45, L 4 A, axis 85, B 2, GAPw nan -> no call (G1+G2+G3+B) | z1 4.9 r1 0.49, L 4 A, axis 69, B 2, GAPw nan -> no call (G1+G2+G3+B) | SEPARATES no; recall/purity g1 0.50/0.0508 g2 0.50/0.0530; arm g1 intermediate g2 intermediate |
| lf_r211_s01 | lf_r211 site_01 | z1 23.0 r1 0.97, L 14 A, axis 6, B 10, GAPw 0.45 -> not discrete | z1 25.4 r1 0.98, L 14 A, axis 2, B 10, GAPw 0.45 -> not discrete | SEPARATES no; recall/purity g1 0.90/0.0899 g2 0.10/0.0106; arm g1 intermediate g2 intermediate |
| lf_r211_s02 | lf_r211 site_02 | z1 2.0 r1 0.28, L 5 A, axis 93, B 2, GAPw nan -> no call (G1+G2+G3+B) | z1 2.1 r1 0.29, L 3 A, axis 76, B 5, GAPw 0.95 -> no call (G1+G2+G3+B) | SEPARATES no; recall/purity g1 0.49/0.0537 g2 0.51/0.0503; arm g1 intermediate g2 intermediate |
| lf_r211_s03 | lf_r211 site_03 | z1 1.7 r1 0.28, L 5 A, axis 74, B 0, GAPw nan -> no call (G1+G2+G3+B) | z1 2.0 r1 0.33, L 5 A, axis 54, B 0, GAPw nan -> no call (G1+G2+G3+B) | SEPARATES no; recall/purity g1 0.46/0.0485 g2 0.54/0.0552; arm g1 intermediate g2 intermediate |
| lf_r21_s01 | lf_r21 site_01 | z1 15.0 r1 0.82, L 7 A, axis 74, B 70, GAPw 1.74 -> no call (G2+G3) | z1 17.4 r1 0.87, L 6 A, axis 67, B 68, GAPw 1.69 -> no call (G2+G3) | SEPARATES no; recall/purity g1 0.10/0.0056 g2 0.90/0.0514; arm g1 open g2 intermediate |
| lf_r21_s02 | lf_r21 site_02 | z1 10.9 r1 0.88, L 8 A, axis 38, B 1, GAPw nan -> no call (G2+G3+B) | z1 12.8 r1 0.93, L 7 A, axis 22, B 2, GAPw nan -> no call (G2+B) | SEPARATES no; recall/purity g1 0.89/0.0500 g2 0.11/0.0066; arm g1 intermediate g2 open |
| lf_r21_s03 | lf_r21 site_03 | z1 9.3 r1 0.83, L 8 A, axis 47, B 5, GAPw 1.20 -> no call (G2+G3+B) | z1 10.6 r1 0.87, L 7 A, axis 43, B 7, GAPw 0.83 -> no call (G2+G3+B) | SEPARATES no; recall/purity g1 0.10/0.0062 g2 0.90/0.0497; arm g1 open g2 intermediate |
| lf_r21_s04 | lf_r21 site_04 | z1 3.1 r1 0.43, L 6 A, axis 60, B 1, GAPw nan -> no call (G1+G2+G3+B) | z1 3.3 r1 0.45, L 6 A, axis 49, B 0, GAPw nan -> no call (G1+G2+G3+B) | SEPARATES no; recall/purity g1 0.51/0.0292 g2 0.49/0.0279; arm g1 intermediate g2 intermediate |
| lf_r21_s05 | lf_r21 site_05 | z1 1.1 r1 0.18, L 9 A, axis 37, B 0, GAPw nan -> no call (G1+G2+G3+B) | z1 1.2 r1 0.19, L 9 A, axis 47, B 0, GAPw nan -> no call (G1+G2+G3+B) | SEPARATES no; recall/purity g1 0.51/0.0293 g2 0.49/0.0278; arm g1 intermediate g2 intermediate |
| lf_r21_s06 | lf_r21 site_06 | z1 4.7 r1 0.60, L 7 A, axis 64, B 0, GAPw nan -> no call (G2+G3+B) | z1 5.0 r1 0.62, L 7 A, axis 58, B 0, GAPw nan -> no call (G2+G3+B) | SEPARATES no; recall/purity g1 0.43/0.0246 g2 0.57/0.0324; arm g1 intermediate g2 intermediate |
| lf_r21_s07 | lf_r21 site_07 | z1 1.4 r1 0.24, L 10 A, axis 38, B 1, GAPw nan -> no call (G1+G3+B) | z1 1.5 r1 0.25, L 10 A, axis 48, B 1, GAPw nan -> no call (G1+G3+B) | SEPARATES no; recall/purity g1 0.55/0.0308 g2 0.45/0.0262; arm g1 intermediate g2 intermediate |
| lf_r21_s08 | lf_r21 site_08 | z1 5.6 r1 0.72, L 7 A, axis 26, B 0, GAPw nan -> no call (G2+B) | z1 5.6 r1 0.73, L 7 A, axis 30, B 0, GAPw nan -> no call (G2+B) | SEPARATES no; recall/purity g1 0.43/0.0242 g2 0.57/0.0331; arm g1 intermediate g2 intermediate |
| lf_r1311_s01 | lf_r1311 site_01 | z1 1.0 r1 0.11, L 4 A, axis 81, B 8, GAPw 1.22 -> no call (G1+G2+G3+B) | z1 1.5 r1 0.16, L 3 A, axis 92, B 22, GAPw 0.81 -> no call (G1+G2+G3) | SEPARATES no; recall/purity g1 0.51/0.0592 g2 0.49/0.0608; arm g1 intermediate g2 intermediate |

## Files

All under `discrete/`.

- `scripts/`:
  - `dc_setup.py`: labels and truth composition;
  - `dc_rec.sbatch`: reconstruct per partition;
  - `dc_extract*.py`: crops;
  - `dc_stats2-5.py`: FAC / AFF / TUBE;
  - `dc_tilt.py`: path stress test;
  - `dc_gap.py`, `dc_gap2.py`: GAP / GAPw;
  - `dc_xorder.py` and `dc_xrun.sbatch`: cross-validated ordering;
  - `dc_xsplit.sbatch`: per-half split, then the ordering;
  - `dc_sub.sbatch`: particle-count subsamples;
  - `dc_qlab.py` and `dc_qrun.sbatch`: data-derived quality orderings;
  - `dc_arms.py` and `dc_arms.sbatch`: a3 / a4 site identification (`crops/arm_sites.json`), orderings, crops, statistics;
  - `dc_battery.py`, `dc_battery.sbatch` and `dc_xorder2.py`: the GAPw call on any mask-defined site / xsb run (section 6);
  - `dc_gates.py`, `dc_gorder.py`, `dc_gates.sbatch` and `gates_table.sh`: the dataset-generic call on xsbg.py gate runs (section 7).
    Outputs are in `out/gates/<run>.txt`, bins in `grec/<run>/`. The 13 pre-fix six-bin attempts are in `grec/_aside/` and `out/gates/_aside/`: they failed at selection because the node's b.simple refers to the full project's stacks, and they were redone on the prep's node.simple.
- `out/`:
  - `v3_*`, `v4_*`, `v5_base.txt`: summaries;
  - `gap2*.txt`: GAPw at tau .5 / .4 / .6, tabulated with `scripts/gap2_table.awk`;
  - `tilt.txt`.
- `rec/<ds>_<stem>/`: bin half maps; `labels/`: per-particle labels and scores.
- Cost and disk: ~50 core-h in 16-core norm jobs (the 6-iteration per-half split is ~45 min per site, a K = 6 bin reconstruction ~1 min). 20 GB on beegfs; nothing deleted; one partial script moved to `scripts/_aside/`.
