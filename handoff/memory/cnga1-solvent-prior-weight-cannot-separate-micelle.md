---
name: cnga1-solvent-prior-weight-cannot-separate-micelle
description: "2026-09-24 MEASURED: master's pcg_solvent Otsu weight w on cnga1 calls 77% of the envelope solvent incl. half the protein; penalty (1-w) is only ~1.4-1.9x larger on the micelle belt than on the protein core, so the solvent prior cannot define/suppress the micelle"
metadata:
  node_type: memory
  type: project
  originSessionId: 2590d10f-7aff-4aff-b2cc-d9d85f962ad3
  modified: 2026-09-24T14:22:36.046Z
---

Replica of `build_solvent_prior_weight` (recipe unchanged 34e69f0c4 -> origin/master bd77d4fec: |rho| of the prior-free
half, `bp(0,lp)` cosine edge 10 idx, lp=max(8, 2*FSC0.143)=8 A, 256-bin Otsu inside support>=0.5, logistic width = sd of
solvent class) on `10_refine3D_auto/recvol_state01_{even,odd}_unfil.mrc` with `automask3D_state01.mrc` as support
(FSC0.143 3.30 A). Script + outputs: `/mnt/beegfs/elmlund/afan/cnga1/solvent_w/` (w.py, report.txt, w_even/odd.mrc),
job 53336790. No batch2 build has bd77d4fec (~/SIMPLE and ~/SIMPLE_masterchk stop at 34e69f0c4) -- none of the cnga1
runs ever used pcg_solvent.

Result: threshold 0.92, 76.9 % of the support w<0.5, even/odd weight corr 0.999. By cylindrical shell (frac w<0.5 / mean w):
0-10 A 0.63/0.38, 10-20 0.56/0.46, 20-30 0.51/0.50, 30-40 0.67/0.36, 40-50 0.81/0.30, 50-60 0.95/0.26, 60-70 1.00/0.15,
70-80 1.00/0.07. So the belt IS in the solvent class, but so is half the protein; ridge (1-w) is 0.50-0.64 on protein vs
0.74-0.93 on belt. **Why:** inside a density envelope there is almost no true solvent, so Otsu splits by density LEVEL
(Wang statistic is for bulk solvent), not protein-vs-detergent.

**How to apply:** do not expect pcg_solvent (any lambda) to separate micelle from protein on detergent/nanodisc data; a
protein/micelle region split needs another statistic (local resolution / high-band e-o correlation, model map, radial cut).
The iteration-count question is moot for the micelle given this discrimination. Related:
[[cnga1-continuum-is-the-micelle-halo]], [[flex-pca-solvent-mask-dominates-rna-basis]].
