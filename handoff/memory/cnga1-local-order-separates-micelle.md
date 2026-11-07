---
name: cnga1-local-order-separates-micelle
description: "2026-09-24 MEASURED: local even/odd correlation in the 6-4.5 A band + Otsu separates the cnga1 detergent belt from the protein cleanly (0-0.4 % of r<40 A called disordered, 75 % at 50-60, 100 % beyond 60) where the solvent-prior density weight could not; candidate knob-free region labels for the complement/region-block fit"
metadata:
  node_type: memory
  type: project
  originSessionId: 2590d10f-7aff-4aff-b2cc-d9d85f962ad3
  modified: 2026-09-24T14:27:12.198Z
---

Script `locres.py` in `/mnt/beegfs/elmlund/afan/cnga1/solvent_w/` (job 53337325, ~minutes, 8 cpu): band-pass the
prior-free `_unfil` halves, local cc in a Gaussian window sigma 6 A, Otsu on it inside the automask; maps in `locres/`
(`lcc_<lo>_<hi>.mrc`, `locres_lcc05.mrc`). FSC0.143 = 3.30 A.

Fraction called disordered (below Otsu) by cylindrical shell, 0-10/10-20/20-30/30-40/40-50/50-60/60-70/70-80 A:
- 8-6 A (thr 0.74): 0/0/0/0/0.01/0.27/0.78/1.00
- **6-4.5 A (thr 0.51): 0/0/0.003/0.004/0.19/0.75/1.00/1.00** <- the clean split
- 4.5-3.5 A (thr 0.30): 0.11/0.21/0.27/0.66/0.92/1/1/1 (starts eating protein)
Local resolution (finest band with lcc>=0.5): 3.8-4.1 A at r<30, 4.4 at 30-40, 4.8 at 40-50, 5.7 / 6.6 / 7.3 A in the belt.

Contrast [[cnga1-solvent-prior-weight-cannot-separate-micelle]] (density-level w called 51-67 % of the protein core
solvent). **Caveat:** order also flags genuinely flexible protein periphery -- at r 40-50, z -60..-40 (outside the belt's
z range -40..+20) lcc is 0.46 < thr; a disordered-label block would absorb real peripheral protein motion too.
Knob-free band choice to decide: tie it to FSC0.143 (6-4.5 A is ~1.4-1.8x) or pick the band with the most bimodal lcc.
Related: [[cnga1-continuum-is-the-micelle-halo]], [[multibody-flex-means-region-pca]].
