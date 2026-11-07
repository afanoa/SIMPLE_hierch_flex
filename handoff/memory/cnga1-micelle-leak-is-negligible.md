---
name: cnga1-micelle-leak-is-negligible
description: "2026-09-24 DECISIVE: on real cnga1 a protein-only support (order-derived mask) leaks only ~1 % micelle variance into its components (latent regression on a jointly fitted belt block); the belt block is unnecessary, the cyl45 'leak' was its mask, and the linear-LS 20 % prediction does not survive the real M-step"
metadata:
  node_type: memory
  type: project
  originSessionId: 2590d10f-7aff-4aff-b2cc-d9d85f962ad3
  modified: 2026-09-24T18:38:24.334Z
---

Runs `/mnt/beegfs/elmlund/afan/cnga1/flex_region/v2_R{1_prot,2_block10,3_block20}` (region binary ~/SIMPLE_region2
90c90a30b, identical recipe, masks from [[cnga1-local-order-separates-micelle]]); `regress.py` / `regress_report.txt`.

**Regression of R1 (protein-only, 10 cols) per-particle latents, cross-validated over eo halves:** on R2's belt block
R^2 0.003-0.022 (variance-weighted **0.010**), on R3's 20-col belt 0.013; on R2's protein block 0.40-0.90 (weighted 0.669,
the noise-limited ceiling); null 0.000; within-R2 protein~belt coupling 0.004. So ~1 % (~1.5 % of the reproducible
content) of the protein-only components is micelle. R1 vs R2 protein canonical correlations 0.97 ... 0.59.

**Edge scores:** R1 edge/inner at belt height 1.07 (TM-edge share 0.31 for 0.30 of voxels) -- no pile-up; R2/R3 0.73/0.74
but that drop is the region branch's 20 A window taper (protein window mean 0.61 inside; relvar 45-50 A 0.030 -> 0.014),
not removed leak. Belt block columns carry lambda 250-360 (2x the protein's 100-170), 1 % energy on protein. R3's extra
belt rank cost protein reproducibility (8/10 -> 5/10 at min|cos| 0.7).

**Why it matters:** the original complaint (micelle-dominated PCs) was the automask support enclosing the belt; the fix is
a protein-only support from local order, no solvent/NU/symmetry prior and no belt block. The cyl45 1.43 edge excess was its
cylinder mask. The linear toy's 20 % leak ([[micelle-linear-leak-toy]]) is an ideal-LS upper bound the regularised
4-iteration M-step does not realise. Supersedes the leak claim in [[cnga1-protein-cut-leaks-belt-variance]].
