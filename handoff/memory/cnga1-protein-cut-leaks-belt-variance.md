---
name: cnga1-protein-cut-leaks-belt-variance
description: "2026-09-24 MEASURED: a hard protein-only support (cyl45_n10, r<45 A) pushes the belt's variance onto the TM surface -- 35-45 A share 18.5% -> 41.8%, absolute ~1.6-1.9x, edge/inner 1.43 at belt height; also order is only weakly a variance measure (voxelwise rho 0.24-0.34 in protein, belt variance follows density +0.58 not disorder -0.36)"
metadata:
  node_type: memory
  type: project
  originSessionId: 2590d10f-7aff-4aff-b2cc-d9d85f962ad3
  modified: 2026-09-24T18:38:38.744Z
---

Script `varmap.py` in `/mnt/beegfs/elmlund/afan/cnga1/solvent_w/` (job 53338085): V(r) = sum lambda_k u_k^2 of the merged
basis (unit-norm, orthonormal) for `flex_pcg_env/sweep_n10` (automask) vs `flex_pcg_env/cyl45_n10`
(`mask_cyl/protein_cyl45_state01.mrc`), identical recipe otherwise (pcg maxits 8, neigs 10, box_crop 128, lp 8,
DEFLATE_BG + EM_DEFLATE=4). Maps in `solvent_w/varmap/`. Sum of eigenvalues 1842 (automask) vs 1416 (cyl45).

**SUPERSEDED 2026-09-24 by [[cnga1-micelle-leak-is-negligible]]: the edge excess below was the cylinder mask / older
binary, not micelle leak (latent regression: ~1 %).** Leak (as first read): % of variance at 35-40 / 40-45 A: automask 10.2 / 8.3 (enrichment over consensus power 1.08x / 1.06x),
cyl45 20.9 / 20.9 (2.22x / 2.66x). Relative variance V/rho^2 at 40-45 A 0.017 -> 0.032, i.e. the TM edge now looks like
the belt (0.031 at 50-55 A). Edge(r40-46)/inner(r30-38) per-voxel variance 0.82 -> 1.13 overall, 1.10 -> **1.43 at belt
height (z -40..20)** -- the pile-up sits where the micelle overlaps in projection, which is what a leak predicts and a
real protein motion would not. So a hard cut is worse than useless on detergent data; the belt needs its own block.

**Order vs variance:** shell-wise, relative variance rises 0.007 -> 0.059 as local 6-4.5 A e/o corr falls 0.93 -> 0,
but voxelwise Spearman(V, 1-lcc) is +0.24 in the protein and **-0.36 in the belt** (there V tracks |rho| at +0.58);
beyond 60 A lcc keeps falling while V/rho^2 plateaus -> weak contrast, not variance. Order is a label, not a variance
prior. Related: [[cnga1-local-order-separates-micelle]], [[cnga1-solvent-prior-weight-cannot-separate-micelle]],
[[cnga1-continuum-is-the-micelle-halo]].
