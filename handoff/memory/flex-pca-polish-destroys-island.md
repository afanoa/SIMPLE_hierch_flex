---
name: flex-pca-polish-destroys-island
description: "10028 joint-fit (polish) forensics, 2026-09-07: ONLINE windows, the plain-prior warm-up, and polish depth >1 all erode the compositional island; recipe that wins on the 20k harness is no windows + mixture prior from it 1 + rank cut + ONE full-N joint iteration (d' 2.48, kNN purity 0.50)"
metadata:
  node_type: memory
  type: project
  originSessionId: 7e382bce-970c-4cdb-bc32-4a1a1c1801a9
  modified: 2026-09-07T07:01:10.731Z
---

Harness: 10028 cryoDRGN poses, first 20 000 particles (1 552 of auto4's 7 881-particle
compositional island), paired recipe neigs=20 MERGE_PRUNE=50, 6 half-fit iterations,
nthr=8 per arm, ~17-25 min each; `~/run_*arms20k.sh` on the box, scored by
`~/island_score.py` (d' = Mahalanobis separation of the island in the delivered latent,
affine-invariant; kNN purity k=15 in the whitened latent; "detached" = fraction of
island particles with >=80% island neighbours). Single-fit arms need
`~/embed_basis20k.sh` because PROBE_ONLY skips the embedding on the single path.

| arm (20k) | d' | purity | detached |
|---|---:|---:|---:|
| single 20-comp, old warm-up | 2.89 | 0.39 | 0.15 |
| single 7-comp, new defaults | 2.76 | 0.58 | 0.39 |
| windowed halves (7600/10k) + windowed polish 2 (the pre-fix recipe) | 0.75 | 0.12 | 0.00 |
| windowed halves, no polish | 1.07 | 0.19 | 0.02 |
| full-pass halves, no polish, old warm-up | 1.46 | 0.30 | 0.13 |
| full-pass halves, polish 2 or 6, old warm-up | 1.25 / 1.23 | 0.21 / 0.19 | 0.03 / 0.01 |
| + mixture prior from it 1, no polish | 1.97 | 0.39 | 0.17 |
| + mixture from it 1, polish 1 (SHIPPED as v3 recipe) | 2.48 | 0.50 | 0.29 |
| + mixture from it 1, polish 6 | 1.96 | 0.38 | 0.16 |
| + mixture from it 1, polish 2 WITHOUT rank cut (rank 20) | 1.43 | 0.16 | 0.00 |

Full N (105 247, nparts=10 nthr=8, `~/run_paired_full_v3.sh`, outputs
`$R28/paired_rankcut_br360_v3`): rank cut 20 -> 17, delivered d' 3.09 / purity 0.54 /
detached 0.34 (auto4 single-fit reference 3.43 on its own labels; pre-fix paired
recipe 1.47). No-polish full-N arm (`paired_nopolish_v3`): 2.98 / 0.53 / 0.35, so at
full N the joint iteration is neutral on the gates and only concentrates the island
onto one axis (best-axis 2.11 vs 1.54). GMM AUTO: 3 macro-clusters; state 6 = 6 582
particles, 86% auto4-island members (72% of the island), 8.0 A; a 25k isolated state
at 4.6 A; a 4-state continuum at 4.0-4.2 A. Wall: paired 39 min + states 14 min.

Spliceosome 10180 (46 783 particles, `$S18/9_flex_pca_paired_br320_v3`, mskdiam 400
lp 20 box_rec 320, sigma2 from the abinitio3D star): the CV cut keeps only 2 of 20
axes with both the v2 (windowed) and v3 recipes; no island, a continuum tiled into
6 GMM states of 2.3k-11.2k particles at 8.4-9.4 A. The 23k-particle halves at 6
iterations may simply be short of budget for more reproducible axes; untested.
Downloads for both datasets: `~/Downloads/flex_paired_20260906/{ribo10028_v3,splice10180_v3}`.

Excluded as causes of the polish loss (each within 0.02 of baseline): merge weights
2c/(1+c) (MERGE_WEIGHT=0), the cross-FSC ridge (XFSC_REG=0), polish depth (6 = 2),
and the ONLINE running average (ONLINE=1 is windows-only; blend is level 2).

Commits on box branch `auto4-rhofix` (~/SIMPLE-h0901-int2): 4fd44e5c2 joint fit is a
full pass (`full_pass` dummy on probe_subspace_iteration; SIMPLE_COV_POLISH_WINDOWS=1
restores for A/B), 93117b344 `n_mix_warm` default 3 -> 1 (SIMPLE_COV_EM_MIX_WARM
overrides). User decision 2026-09-07: ONLINE windows are gone for good
(SIMPLE_COV_ONLINE=0 in every recipe); the half-set -> joint structure stays.

**How to apply:** score any change to the paired chain with island d' AND kNN purity
on this harness before a full run; UMAP alone misleads (the 2.89 single fit shows no
detached cluster, the 1.46 merged basis does). Polish 1 is a transient optimum, not a
converged one; deeper joint fits regress to the merged level. Never drop the rank cut
before the joint fit. Related: [[flex-pca-split-half-checkerboard-bias]],
[[flex-pca-transport-rank-is-not-the-barrier]] (variance ordering, consistent with
rank 7 beating rank 20 on detachment).
