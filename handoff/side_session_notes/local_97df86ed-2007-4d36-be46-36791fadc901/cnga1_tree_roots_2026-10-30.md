# CNGA1 root-only hetero_tree (session tree-micelle-synth, 2026-10-30; W = hetero_tree_20261028/cnga1_tree)

## Build: cnga1_tree/build (norm SLURM build, batch2 login down; fast gate 15/15)
- Source: fork hetero-tree 21fac34db + halves_halfseed_r17d_v2.diff (md5 4f7724ed; = route9/halves/v2/src, hetero_tree md5 f2c35513).
- On top of it, cnga1_tree_on_halves_v2.diff (md5 5ced43f7, 2 comment lines), every part opt-in:
  - tree_mskfile (user support mask: node focus / site masks and the root flex stay inside it).
  - pcg_bgmskfile + user_mask_keys: unmasked reconstructions solve inside the support; focused ones solve their consensus background inside it.
  - SIMPLE_TREE_FLEX_C1=yes (test switch): flex_keys drops pgrp.
- With no new key or env set, the build is exactly v2.

## Why SIMPLE_TREE_FLEX_C1 exists
- hetero_tree's common_keys sets pgrp on every step, and flex_keys calls it, so all three flex_pca calls (root fflex, the per-half flex, the gate refits) get the tree's pgrp.
- flex_pca with pgrp=cN silently sets nsym=N and symmetrises every accumulation (memory cnga1-is-a-continuum-not-discrete-states).
- So every OLE tree (pgrp=c2) has run flex_pca with nsym=2. That is unchanged here (the switch defaults off) and it is what the gates measure.
- For CNGA1 the C4 arms set the switch.

## Inputs (cnga1_prep, job 54728642; masks/)
- cnga1_test.simple: 245,424 active rows, eo 122,712 / 122,712, e1 in [0,90) 99.97 %, defocus 0.42-3.43 um.
- minpop 3,681 (1.5 %).
- Masks (dilated 4 A, edge 6 A), as automask kept by cylindrical shell from the C4 axis; columns are 0-10 / 10-20 / 20-30 / 30-40 / 40-50 / 50-60 / 60-70 / 70-80 A:

| Mask | Kept by shell | Consensus power kept | Dice with locorder |
|---|---|---|---|
| support_c2 (band 6.8-4.5 A) | 1.00 / 1.00 / 1.00 / 1.00 / 0.98 / 0.74 / 0.18 / 0 | 0.923 | 0.960 |
| support_c3 | 0.99 / 0.93 / 0.94 / 0.91 / 0.66 / 0.17 / 0 / 0 | 0.792 | 0.764 |
| locorder (Sept. lcc 6-4.5 A, Otsu 0.511) | 1.00 / 1.00 / 1.00 / 1.00 / 0.99 / 0.74 / 0.14 / 0 | 0.920 | 1.000 |

- Choice for the support arms: support_c3.
  - c2 and the local-order cross-check keep 74 % of the 50-60 A shell, where the halo variance is 3.5x enriched.
  - c3 keeps 17 % there and none beyond. Its cost: 34 % of the 40-50 A periphery and ~4 % of TM consensus power.
  - MASK=masks/support_c2.mrc switches the support arms to c2.

## Launch (one command; waits for the coordinator's go)
- `sbatch $W/scripts/cnga1_launch.sbatch` submits all four roots: c1_nomask, c1_support, c4_nomask, c4_support.
  - Keys: nparts=10 mskdiam=180 site_layer=yes site_states=yes kway_maxdepth=2 minpop=3681, maxdepth 0, nodepar 4, seed 1, norm.
  - Env: CONTRAST=maps, RHOMIN unset, HALVES=yes, HALFSEED=yes, KEEPMAPS=yes.
  - ARMS='c1_nomask c1_support' submits a subset.
- Readout after: `sbatch $W/scripts/cnga1_check.sbatch`. It prints the split lines, where the reproducible stage-B group difference sits by shell (pore / TM / periphery / belt / outside), belt / TM density per group, and AUCs for defocus, astigmatism and input corr.
- Price, scaled from the 692k RNase P root (134 core-h) with the measured halves + halfseed overhead (+10-50 %): ~45-75 core-h per root, ~250 for all four.

## OPEN ISSUE (whole tree, not only CNGA1): the tree hands its pgrp to flex_pca
- Every OLE tree and OLE gate (pgrp=c2) has run flex_pca with nsym=2, a symmetrised basis that nobody chose.
- Default left as it is (coordinator 2026-10-30). The coordinator raises it with the user.
- A default change needs the OLE gates rerun. Memory: hetero-tree-passes-pgrp-to-flex-pca.
- Status: all four roots approved at ~250 core-h with support_c3; launch waits for the coordinator's "go" after the split/merge gates.

## Launch history (2026-10-30)
- GO from the coordinator with SIMPLE_TREE_HALVES / HALFSEED UNSET: they failed the r131 closed-corner hard gate.
- Launch 1 failed in 11 s: "particle 8 has different halves in ptcl2D and ptcl3D".
  - The Sept project's ptcl2D holds the 2D halves; they agree with refine3D_auto's ptcl3D on 49.4 % of rows. No CNGA1 project has a congruent ptcl2D.
  - Fix: test switch SIMPLE_EO_ADOPT3D=yes (eo_adopt3d.diff, adopt_legacy).
- Launch 2 failed in seconds: "nptcls_stk not present". The 2022 stack has 274,333 images for 245,424 rows, and fix_projfile distrusts legacy indstk.
  - Proof that the stored indstk is the physical image:
    - each row's CTF equals import row indstk(i) exactly (100 % vs 0.02 % for row i);
    - Thon rings correlate with the CTF at 0.60 / 0.34 when read by indstk, and -0.13 / -0.34 when read by row.
  - Fix: test switch SIMPLE_FIX_TRUST_INDSTK=yes (fix_trust_indstk.diff) -> input/cnga1_ready.simple. 2 repairs, 0 errors, all fields identical.
- Launch 3 (job 54733188): running on cnga1_ready.simple with env FLEX_C1=yes, EO_ADOPT3D=yes, CONTRAST=maps, KEEPMAPS=yes.
- Build: incremental rebuilds, fast gate 15/15 each. Combined diff on v2: cnga1_tree_on_halves_v2.diff (md5 d52b2153; 4 opt-in pieces, one comment line per switch).

## Results (launch 3; readouts cnga1_check / cnga1_bmatch; regions on stage-B maps at box 88 / 2.76 A)
- Region definitions:
  - pore/S4-S6: r<22 A in the slab z -40..20;
  - VSD: 22-45;
  - belt: 45-80;
  - CNBD side: z>20, the cytoplasm by protein mass;
  - extracellular side.
- Enrichment = share of the reproducible cross-half between-state variance / share of the consensus power.

| root | K* (sizes) | z / z_crit | guard L/L_ref | where | nuisance | verdict |
|---|---|---|---|---|---|---|
| c1_nomask | 4 (80k/64k/53k/49k) | 13.10 / 3.38 | 0.93-1.25 | belt 3.33x at 40-15 A; 15-8 A at noise | corr AUC 0.35/0.59 | micelle halo + quality |
| c1_support | 4 (67k/66k/61k/50k) | 4.97 | 1.24-1.52 | spread like the protein (VSD 1.31x, CNBD 1.04x, pore 0.83x) | tilt AUC 0.255 / 0.63, view NMI 2x | view imbalance |
| c4_nomask | 2 (185k/60k) | 24.81 | 1.233 | 40-15 A pore 2.2x, survives B-match; r<20 share 0.27 (gating 0.71) | 60k: corr 0.24, pose change 0.71, dB +160 | worse particles |
| c4_support | 2 (159k/87k) | 13.85 | 1.213 | 40-15 A pore 2.1x, survives B-match; r<20 0.24; 15-8 A removed 2/3 by B-match | 87k: corr 0.315, pose 0.664, dB +193 | worse particles |

- Defocus, astigmatism, CTF res and contrast (ice proxy) are at chance everywhere. No maps delivered.
- Cost: roots 408.6 core-h, OVER the approved 250: not tracked in flight. The C4 post-split steps were cancelled by ID once noticed. Builds, prep and checks ~6.
