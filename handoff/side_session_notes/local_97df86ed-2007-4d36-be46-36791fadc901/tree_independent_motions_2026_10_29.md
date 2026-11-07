# tree-independent-motions (2026-10-29): several independent motions vs the hetero_tree

## Fixture mm3 (batch2 routes/toy55)
- `make_mm_maps.py`: base = RNase P conf01 (maps160/docked_all), plus 3 rigid textured domains (18x13x12 A) on the rim, 120 deg apart in the
  plane of the long axes. Each domain has 2 hinge states (+-25 deg in-plane), independent populations P(state 2) = 0.50 / 0.65 / 0.80, so
  8 combinations at product weights (1,890 to 14,040 of 54k genuine).
  - Per-domain state difference / map power: lp10 0.12-0.15, 40-15 A 0.22-0.29, 15-10 A 0.08-0.13.
  - For scale, qtoy2 conformer pairs: c1-c22 0.128 (F1 0.78 in r17b), c1-c28 0.216 (F1 0.94).
  - Optional `--cont DELTA:NBIN` adds a continuous domain.
- `run_mm_sim.sbatch` (QT=mm3): the qtoy simulator with generated maps, consensus = `maps160/mm3/consensus_mm.mrc`.
  - Settings: N 60k, JF 0.05 (two junk classes of 3k each), SNR0 0.0012, refine3D from 10-deg-perturbed poses.
  - Outputs: `simmm3/src.simple`, truth `simmm3/truth_qtoy.txt` (combo = 3-digit id, e.g. 122 = f1 state 1, f2 state 2, f3 state 2), and
    `maps160/mm3/factors.txt` plus `domain{1,2,3}_mask.mrc` (masks for scoring only).
  - Trap: Python's csv writer adds `\r` at line ends, and the `while read` loop then fails on every map path. Fixed with lineterminator.
- `mm_score.py TREE_DIR` reports:
  - per-leaf factor purity;
  - best-leaf F1 per combination and its fragmentation (leaves holding >= 10 %);
  - which factor each split follows (Cramer V >= 0.3), and repeats (deepest nodes splitting the same factor = disjoint branches);
  - per factor state: best leaf / node / unit F1 and leaf-majority accuracy;
  - merge-up dedup vs mixed fusions.
- `mm_leaf_factors.py TREE_DIR`: post-hoc factorization of the leaf even/odd maps.
  - Method: cross-half rank, spatial ICA, leaf codes, truth overlap with the domain masks.
  - Checked on a synthetic 4-leaf tree: 3 components, each 0.89-0.93 inside one domain mask.

## Python toy (mm_toy.py; scratchpad, copy to follow)
Setup: 2D maps, 60k particles, the tree emulated with a K-way root to depth 2, binary below, held-out z > 3, minpop 900, maxdepth 4.
Four representations are compared:
- TREE: the emulated tree.
- FLAT: one flat K-way partition.
- TPF: tree, then factorization of the unit maps, held-out pruning, and a factorial E-step with a free joint prior.
- DF: the same factorial model started from the particles' PCA.

Three independent binary factors (populations 0.5 / 0.65 / 0.8), grid 1:
| d' | tree leaves | per-factor acc, tree | per-factor acc, TPF | map contrast, tree | map contrast, TPF |
|---|---|---|---|---|---|
| 1.0 | 3 (f3 lost) | .66/.66/.80 | .68/.67, f3 missed | 39/34/3 | 44/41/0 |
| 1.5 | 5 | .74/.77/.83 | .77/.77/.75 | 48/41/36 | 65/60/43 |
| 2.5 | 7 | .88/.89/.89 | .89/.90/.92 | 73/61/40 | 95/90/73 |
| 4.0 | 9 | .86/.96/.97 | .98/.98/.98 | 76/83/51 | 116/111/92 |

- Map contrast = the factor's purity difference between two delivered maps x sqrt of the harmonic N.
- At d' 1.0, DF finds f3 (contrast 27); TPF misses it because no unit varies in f3.
- At d' 4.0 the tree splits f2 in 4 disjoint branches and loses f1 in some of them.

Controls:
- One motion, and three domains moving together: one factor, identical to the tree.
- 4-of-8 combinations (coupled populations): the free joint prior concentrates on the 4.
- One long curved continuous motion splits into 4-5 components (|r| up to 0.94). Grouping dependent components (label NMI > 0.05 or
  |r| > 0.2) gives 1 motion at the tree's eta^2 / tail F1 (0.824 vs 0.825; 0.59/0.63 vs 0.62/0.63).

Grid 2: the final version, with held-out pruning (drop a factor unless it raises the held-out likelihood, z >= 3) and
dependence grouping.
- 3 independent binary motions, d' >= 1.5:
  - TPF finds exactly 3 factors in 9/9 runs (3 seeds x d' 1.5 / 2.5 / 4).
  - Total correlation ~0.0002 nats.
  - Combination NMI, TPF vs tree: 0.198-0.199 vs 0.190-0.194; 0.487-0.490 vs 0.456-0.463; 0.826-0.833 vs 0.706-0.801.
  - Per-motion map contrast, TPF vs tree: 65/60/43 vs 47-52/41-50/36-38; 95/90/73 vs 72-73/53-61/40; 116/111/92 vs 76-86/71-83/50-51.
- mm2c (2 binary + 1 continuous):
  - Contrast 95/86 vs 59/55 at d' 2.5, and 116/106 vs 62/55 at d' 4.
  - The continuous motion's 10 % tails F1: .29/.30 vs .23/.27 (d' 2.5); .40/.42 vs .28/.28 (d' 4).
  - The continuous motion takes 2 components, grouped back into one motion.
- partial (4 of 8 combinations): grouped into 1 motion at d' >= 2.5 (total correlation .53-.67), identical to the tree. At d' 1.5 it
  loses (combination NMI .243 vs .285).
- cont1 (one curved continuum): 1 motion at the tree's eta^2; tails slightly lower (.55/.58 vs .63/.64 at d' 4). At d' 1.5 the grouping
  misses (2 motions) and eta^2 drops .209 -> .127.
- one, coupled: 1 factor, identical to the tree.
- d' 1.0: the tree never splits f3, so TPF misses it too. DF (the particles' own PCA) finds it (contrast 27).

## The production site layer already does part of this (site_layer=yes, hetero_tree factorial())
- A site found at a deep node climbs to the parent level when its pole ranking reproduces there (Spearman >= FGATE 0.5).
- Sites found in different branches are marked copies (F_COPY) or duplicates (F_DUP): DAG sharing, for sites.
- Output: factorial/table.txt, per kept site the pole labels and posteriors over its level (often the root), plus pole maps.
- The difference from TPF: one-sided poles (pole vs rest), no combination maps, no re-assignment against the factor model.
- mm_score.py section 6 scores the poles against the truth factor states (checked on qtoy2: r3_s1 = closed, F1 .85).

## Real tree on mm3q
- Sims: mm3q (SNR0 0.0048 = qtoy2's noise; poses 2.60 deg median, p90 4.71) and mm3 (SNR0 0.0012 = the null's noise; 3.59 / 8.28).
- Tree: recipe17d floor-free via sst_gate4, qtoy2 keys + site_layer=yes, maxdepth 4, nodepar 4, seed 1.
  - Job 54667841, run dir hetero_tree_20261028/multimotion/runs/mm3q_r17d_s1, root z_crit 3.38 (132 tests).
- Helpers: mm_node_score.py (per-node split vs the factors while the tree runs).
- mm_fassign.sbatch (scratchpad, not yet on batch2): one stage-B-keyed refine3D pass against the 2^r combination maps from
  mm_leaf_factors.py synth=yes.

## mm3q root (2026-10-29 ~18:40)
- Site discovery: 5 root sites. Site 1 = d2, sites 2 and 3 = d1, sites 4 and 5 = d3. The node focus mask covers all three domains.
- Strongest pole: site 5 (11,990 = 0.20).
- Route A K* = 3. Cramer V: f2 .89, f3 .60, f1 .26 (the 50/50 motion is not split at the root).
  - r1 (27,088): 122 + 222.
  - r2 (19,736): 112 + 212.
  - r3 (13,176): f3=1 .65 -> 221 / 222 / 121.
- So f1 must be rediscovered in r1, r2 and r3.

## IgM armsigm1, offline, on the recipe17b igm_s1 tree (83 leaves; no new compute)
- mm_leaf_factors.py rank=16: one component per arm angle for all 10 angles. |r| of the leaf code with the leaf-mean angle .73-.98;
  leaf-mean R^2 .83-.98.
- Trap: the cross-half rank rule gives 82/82 on real trees. Selection imprint keeps every even x odd eigenvalue positive, so the factor
  count has to come from a particle-level held-out test.
- Per-particle ordering of each arm angle (igm_compare.py, genuine rows):
  - leaves (leaf-mean predictor) .54-.74;
  - one leaf-map factor code (leaf-constant) .42-.73, all 10 angles;
  - best root-level site-layer pole posterior .48-.76 on 8 angles, missing a3_th (.13) and a4_ph (.14).
- Reading: the leaves carry every arm. Offline, the factor layer adds completeness and per-motion naming; accuracy needs the
  re-assignment.

## mm3q depth 1 (~19:00)
- f1 is rediscovered in all 3 branches. Cramer V: r1 .96, r2 .96, r3 .90.
- r3 (K-way, 4 groups) also splits f3 (V .68) and isolates junk: r34 = 917 rows, 100 % junk.
- Groups now: r11 122 82 %, r12 222 89 %, r21 212 70 %, r22 112 75 %, r31 221 51 %, r32 121 74 %, r33 222 87 %.
- 222 leads two groups in different branches (r12, r33).
- 111 and 211 (1,890 each) are leaves nowhere; they sit at 10-14 % inside r21, r22, r31 and r32. In r2, f3 is not split (V .04).

LATTICE EVIDENCE, from maps only: stage-B child difference maps (lp 12, envelope), cross-half cosines (even of one x odd of the other).
- f1's three splits: r1:2-1 vs r2:2-1 -0.98, vs r3:2-1 -0.96; r2 vs r3 +0.95. Each has 78-81 % of its energy inside domain 1.
- Splits that follow other motions: |cos| .02-.40 against these.
- r3's f3 split (3-2, d3 .61) vs the root's f3-heavy contrast (1-2, d3 .62): -0.74.
- So equivalent sub-splits in different branches can be recognized without truth (|cos| >= .95), and could be tied: the DAG/lattice.

## mm3q FULL TREE (job 54667841: 1 h 31 min, 69.2 core-h; run dir multimotion/runs/mm3q_r17d_s1; score in scratchpad mm3q/score_tree.txt)
- 16 leaves: 14 genuine (median 4,242, smallest 1,061) and 2 junk. The merge-up fused nothing (16 units).
- Combinations: 1 of 8 is a clean leaf (112, F1 .89); the others reach .51-.65. 222 is spread over 4 leaves, 111 over 4, 122 over 3,
  121 over 3.
- Repeated splits: f1 in 3 disjoint branches (r1, r2, r3); f3 in 3 (r3, r21, r22); f2 only at the root.
- 5 accepted splits follow no motion (r11, r12, r32, r111, r2111).
  - Each cuts a single combination in two: 122 into r1111 / r1112 at ~88 % each; 222 into r121 / r122; 212 into r2111_s1 / s2.
  - Not quality: snrq / bfac AUC .47-.51.
  - The merge-up kept all of them.
- Per-motion state: best leaf F1 .28-.54; best node .50-.88; leaf-majority accuracy .94-.98.
- Site layer, 7 kept sites: 3 root-level sites, one per motion. r_s1 = f2 state 1 (F1 .92), r_s3 = f1 state 1 (.85),
  r_s4 = f3 state 1 (.82). Their poles are mutually independent (NMI .001, |r| <= .05). The 4th root-level site, r33_s1, is a weak
  view of f2 (r -.20 with r_s1).

## Post-hoc representations on mm3q (genuine rows; truth only scores)
Every "pass" below is one stage-B-keyed refine3D pass (mm_fassign.sbatch, ~1.5-2.5 core-h, 4-7 min).
Columns: combination NMI | per-combination best F1 (mean / smallest / combinations at F1 >= .7) | per-factor accuracy / contrast.
- tree leaves (15 states): .647 | .65 / .51 / 1 | .975/54 .948/58 .938/39
- control, one pass vs the 16 leaf maps: .673 | .72 / .57 / 4 | .978/56 .962/60 .954/42
- site layer, 3 root poles combined, NO pass: .732 | .81 / .73 / 8 | .943/103 .972/104 .959/81
- site product + one pass (mm_site_combos.py): .780 | .82 / .67 / 7 | .976/111 .965/102 .957/81
- site product with the redundant 4th site + one pass: .700 | .75 / .53 / 5 (a redundant site hurts: the duplicate filter matters)
- leaf-map ICA rank 3, pooled, no pass: .748 | .75 / .59 / 4 | .975/110 .948/98 .938/75
- leaf-map ICA rank 3 + one pass: .801 | .82 / .64 / 7 | .978/111 .971/104 .962/82

Cost so far: sims 39.2 + tree 69.2 + passes 7.4 = 115.8 core-h.

## RECOMMENDATION
- A partition tree alone is the wrong representation for independent motions: fragmentation, repeated sub-splits, and no leaf for a
  motion's state.
- But hetero_tree already holds a second, factorial structure: the site layer finds each motion once, at the root.
- So keep the tree, and add an opt-in post-hoc product step over the site layer (default off, so RNase P runs stay byte-identical; not
  proposed or run on RNase P):
  1. Take the kept sites pooled at the root.
  2. Drop or group dependent sites. Pole NMI / posterior correlation, or a cross-half cosine of their difference maps (the lattice test:
     same motion >= .95, others <= .40).
  3. Product labels: per particle, its pole bits give the combination.
  4. Optionally one stage-B pass against the 2^K combination maps synthesized from the leaf maps (+.05 NMI, ~2 core-h at 60k).
  5. Deliver: per-motion state maps (the site layer's pole maps already exist) and combination maps for combinations >= minpop.
- Leaf-map spatial ICA gives the same 3 factors, and all 10 IgM arm angles where the site layer has 8. Keep it as a diagnostic only:
  region-basis-like.

## CONTROLS: site-product readout on every non-RNase-P tree with a site layer (2026-10-29 evening, 2.8 core-h new)
Tool: routes/toy55/site_product_readout.py RUN tag=.. [truth= cols=] [fields=] [labels_out=]. Logic:
1. Kept root-level sites.
2. GATE, from the site layer's own factorial/sites.txt: >= minpop hard members (posterior >= .5) and decisiveness = hard / soft pole
   size >= .5.
3. Dependence: pole NMI > .05 or |posterior r| > .2 -> one motion.
4. Product of the motions' pole bits.

Gate decisiveness: real and continuum sites .81-3.29; fakes <= .40 (qtoy2 quality pole .05, broken-junk pole .22, mm3q duplicate .12,
IgM junk-leaning .40, nulls 0). The .5 threshold was set AFTER seeing these controls; the gap is wide, but it needs a held-out check.

| tree | root sites | motions after gate | expected | notes |
|---|---|---|---|---|
| mm3q | 4 | 3 | 3 | duplicate rejected; product NMI .732 vs leaves .647 |
| qtoy2 r17b / md2 / inherit | 6 | 1 | <= 1 | arm swing (4 sites grouped, closed_arm rho -.50 to -.61); quality + junk poles rejected; product NMI with conformer .241 vs leaves .706 (closed/holo are leaf-only states: the product is no substitute for leaves) |
| ctpk, ctun (base, inherit s1/s2: 6 runs) | 1 | 1 | 1 | rho with t -.41 / +.57; no view / defocus association |
| nulls with a kept site: qtoynull 17b/inherit, olenull r17d/17b/hier, micnull fp_order/gate_order | 1 | 0 | 0 | 0 hard members: WITHOUT the gate, 7 null runs would get one fake factor |
| nulls without root sites: qtoynull r17d/hier, N000, N00_12k | 0 | 0 | 0 | |
| micelle mic3 fp nomask / order / gate, micnull fp nomask | 0 | 0 | 0 | the protein states the tree misses stay missed |
| rare state L20_6k, L10_12k (md1 and hier md2) | 0 | 0 | - | no fake factor, no help: rare states are leaves or nothing |
| IgM armsigm1 | 11 | 8 | up to 10 arm angles | every motion is a real arm angle (|rho| .48-.76); a3_th and a4_ph missed; 0 fake |
| OLE r14 and r17b, s1 and s2 | 3 | 2 per seed | ? | both reproduce across seeds (A |r| .89, B 1.00/.78); not views/defocus (|rho| <= .07, view-grid NMI <= .004) |

OLE cross-seed reproducibility (80,599 rows):
- product (4 combinations) .536; root children (5 vs 5) .523; leaves (21-26) .324-.328.
- One stage-B pass against the site-product maps (2 x ~1.4 core-h): .462. Motion A .402 -> .516; motion B .743 -> .405 (29 % relabelled,
  minority pole 18 % -> 43 %).
- So: no re-assignment pass on real data; the site layer's own poles are the reproducible labels.

## HELD-OUT CHECK of the gate (2026-10-29 late; offline, 0 core-h)
- Trees not used to pick .5: qtoy2 (r14, r16, ancestral, dsnr root on/off, the qbp_r1211 fixture), 14 more nulls (r14/r16/r17e,
  route9 s2-4, dsnr on/off), IgM (r17f, r14, the m69M r111/r12 fixtures, dsnr igm_r121 on/off), micelle non-fp (mic3/micnull
  nomask/order/umask), 16 more detlim fixtures, the OLE r32 node fixtures.
- 67 root sites in all: 24 repeat seen sites (same numbers, all called right); 43 are new.
- Truth class fixed in advance (gate_heldout.py):
  - qtoy2 conformer NMI >= .05 or |rho arm| >= .3, with pole junk < .5;
  - IgM max |rho angle| >= .3, junk < .5;
  - detlim rare F1 >= .3;
  - mic3 state NMI >= .05;
  - nulls fake.
- New sites: REAL and PASS 25; FAKE and rejected 14; FAKE but PASS 3; REAL but rejected 1.
  - FALSE PASS igm_r17f r_s4 (dec .91) and igm_r14 r_s4 (.76): weak a3_ph sites (|rho| .23). Both are DEPENDENT on the strong arm-3 site
    and merge into that motion, so the product gets no fake motion from them.
  - FALSE PASS igm_r14 r4_s3 (dec .71, 23,800 hard members, pole 94 % junk): a decisive JUNK pole. It is a 7th "motion" in that
    product, and the one product-level false factor in the held-out set.
    - The tree's own junk screen found the same thing: the root log reports an amplitude-rule junk component of 24,716 particles
      ("not applied").
    - So an in-tree step can drop a site whose pole is mostly amplitude-rule junk (node site_junk.txt).
    - It can't be verified offline: finished trees prune site_junk.txt.
  - MISS igm_r17f r1_s3 (dec .47, a3_th |rho| .32): a weaker copy of a3_th, which r2_s7 still carries, so the product loses nothing.
- Product level, held-out:
  - IgM r17f: 6 motions, all real.
  - IgM r14: 6 real + 1 junk.
  - qtoy2 trees: 1 motion (the arm), with the quality and junk poles rejected.
  - Nulls, micelle, detlim: 0 motions. Every null site has 0 hard members; mic3 / L025 / P_L10 sites have 0 hard members too.
- Also fixed: site_product_readout.py now parses an overflowed halves z ("*******") instead of silently rejecting the site.
- recipe18c wide-test trees: still running at write time. Readout command:
  site_product_readout.py RUN tag=.. truth=.. cols=.. (or gate_heldout.py with a runs file).

## FINAL REVISED STEP (opt-in key, default off; RNase P runs byte-identical; never proposed or run on RNase P)
Runs after the site layer (factorial()) in the tree controller, from files the tree already writes:
1. Candidates: the kept sites pooled at the root (factorial/table.txt level r).
2. Gate:
   - >= minpop hard members (posterior >= .5);
   - decisiveness = hard members / soft pole size >= .5;
   - pole NOT mostly amplitude-rule junk (the root junk screen's site_junk.txt, read before the end-of-run prune).
3. Dependence: pole NMI > .05 or |posterior r| > .2 = one motion; the strongest site (halves z) represents it.
4. Product labels: per particle, the motions' pole bits give the combination id. Write hetero_tree_product.txt and a motion table.
5. Maps: the site layer's per-motion pole maps (already written); combination maps only for combinations >= minpop (optional).
6. NO re-assignment pass. It hurt OLE reproducibility (.536 -> .462).
- Leaves stay the primary deliverable: states that are not a product of motions (closed, holo, junk, composition) live only in leaves.
- Evidence:
  - mm3q combination F1 .81 vs leaves .65.
  - OLE: 2 motions reproduced across seeds; cross-seed NMI .536 vs root children .523 vs leaves .33.
  - IgM: 6-8 real arm motions per tree.
  - Controls (nulls, continua, micelle, rare state): correct counts.
- Open:
  - a held-out confirmation of the junk rule;
  - recipe18c readouts;
  - the factor count is only as complete as the site layer (IgM misses 2-4 of 10 angles).
