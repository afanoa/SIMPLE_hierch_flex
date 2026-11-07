# Island site-discovery rule (change (b)) on the hetero_tree gate and control datasets — INTERIM (reruns still running)

Session island-gates, 2026-11-05. Offline, Python only, no Fortran, nothing committed. Cluster work dir:
`/mnt/beegfs/elmlund/afan/igm_cryolo/tree/fixture_igm/island_gates/` (scripts/, logs/, out/<node>_bg{no,yes}/{log.txt,summary.json,vs.mrc,sites/}, reruns/, fixtures/, truth/).

## Method
- `scripts/island_regions.py` = a copy of `fixture_igm/scripts/proto_regions.py` (region finding and factor test unchanged; original not edited) plus:
  the dataset's `mskdiam` (proto hard-codes 600), the tree's own candidate mask (add(): 6-connected dilation by MDIL, falloff MFALL),
  counts from the node's `sites/disc_labels.txt` (the labels the tree used), a mirror check against the node's own SITE_DISCOVER log,
  truth placement (RNase P closed / open arm masks by arm_readout.py's rule; mm3q domain masks; continuum swing region; mic_het TM-rock / Fv regions).
- Mirror check: stock core components match the tree's logged components to <= 1 voxel on every node (OLE c2: <= 28 voxels of 17k-122k, float precision at the cut);
  on the IgM fixture the copy reproduces the fixture session's proto run exactly (fix: 6 / 7 islands, 4 / 5 island sites without / with offsets removed).
- Two readings of "every region is a site candidate, kept when the SVD factor test passes":
  - ADDITIVE = the tree's own discovery unchanged + every island whose factor test keeps >= 1 factor;
  - LITERAL = every region, stock pieces included, gated by the factor test.

## LITERAL reading: NO-GO (hard gate)
At RNase P's 1.46 A/px most stock regions hold n_eff 1-10 (15 A elements), so the factor test keeps no factor even at cross-half r .5-.97.
From the production 692k tree's own discovery logs (f692_lf_s1, build_lf; maps pruned but logs kept):

| node | rows | tree sites | sites whose component keeps 0 factors (literal drops) |
|---|---|---|---|
| r | 692,392 | 4 | 0 |
| r1 | 162,926 | 4 | 0 |
| r13 | 46,991 | 5 | 3 |
| r131 | 24,897 | 1 | 0 |
| r1311 | 12,470 | 6 | 6 (all) |
| r1312 (closed leaf r1312_g1_2) | 12,427 | 3 | 3 (all) |
| r132 | 22,094 | 4 | 3 (incl. SITE 1) |
| r1321 (closed leaf r1321_g1_2) | 12,044 | 5 | 4 |
| r2 | 145,632 | 4 | 0 |
| r21 | 53,412 | 2 | 0 |
| r211 | 27,417 | 4 | 2 |
| r2112 (closed leaf r2112_g1_2) | 13,202 | 3 | 3 (all) |

The literal form would leave two of the three closed leaves' nodes with no site at all. Only the ADDITIVE form is a candidate; everything below is that form.

## ADDITIVE reading on kept-map nodes (gate replays, build_bc / build_head; bg no / bg yes = offsets kept / removed)
- Nulls (md0 roots): qtoynull (60k, 20 tree sites), olenull (80.6k), micnull (56.8k): 0 islands, with and without offset removal -> no change.
- r131 fixture (the r17d production tree's closed node, 44,334 rows), seed 2: roots, r1 (corner, 1.00/0.49) and r2: 0 islands -> no change.
- r131 seed 1: root 0 islands; r2 (corner, 0.97/0.53): 1 island, factor test fails (0 with offsets removed) -> no change;
  r1 (23,449 rows, 1,582 of the 3,075 closed reference rows, split currently not separating 0.78/0.59 vs 0.92/0.62): 1 island site
  (3,615 voxels, peak 0.18 x max, 16.5 % of the node variance, centre 60 A, molecule .106, kept by factor 2: share .42 / r .93 / z 3.4; on neither arm mask);
  with offsets removed a second (centre 30 A, molecule .276, z 3.3). -> CHANGE at a non-corner node of the closed lineage.
- qtoy2 (60k, md1): root, r1: 0 islands; r2 / r3 / r5: 2-3 islands each, none passes; r4 (9,461): 1 island site (8,058 voxels, 29 % of the node variance,
  centre 36 A, molecule .302, z 4.5 / 4.9, off the closed arm) — over DMAXMOL 0.3, so the tree's own non-locality filter would drop it if a port keeps that filter.
- mic_het (2 seeds x root + 2 children): 0 islands everywhere; the stock core already covers both moving parts (root: TM rock .94, Fv 1.00 in one component).
- OLE (real, 80.6k, md1): root, r1, r4: 0 islands; r2 / r3: 2-3 islands, none passes -> no change.
- IgM fixture (validation): fix root 4 island sites (5 with offsets removed) on 4-5 arms; null root 5 island sites (its arms move too).

## Pending (root-only reruns, maxdepth 0, maps kept, build_lf, site_layer=no, seed 1)
mm3q, ctun, ctpk roots; production 692k nodes r13, r21, r131, r132, r211, r1311, r1312, r1321, r2112 (+ r1, r2); the 692k root (approved).
