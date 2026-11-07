---
name: cnga1-angular-harmonic-label-weak
description: "2026-09-24 MEASURED: the angular-harmonic (m>=3 about the membrane normal) micelle label does NOT separate cnga1 protein from belt (71% of protein called smooth); blind axis TILT is recovered (0 deg, clean landscape) but the CENTRE is not (-9,-5 A, flat +-10 A) because removing m=1 absorbs shifts"
metadata:
  node_type: memory
  type: project
  originSessionId: 2590d10f-7aff-4aff-b2cc-d9d85f962ad3
  modified: 2026-09-24T14:52:13.264Z
---

Script `harmonic.py` in `/mnt/beegfs/elmlund/afan/cnga1/solvent_w/` (job 53342206, 4 cpu): Cartesian ring bins (1 px in
r and along the axis), subtract each ring's m<=2 expansion, local m>=3 fraction of ring power (8 A consensus, Gaussian
4 A), Otsu. Map `solvent_w/harmonic/highharm_frac_knownaxis.mrc`.

Known axis: mean m>=3 fraction by shell 0.22/0.24/0.22/0.26 (protein 0-40 A), 0.17 (40-50), **0.066 (50-60)**,
0.14/0.19 (60-80, weak density -> noise). Otsu 0.266 calls 77.6% of the automask smooth, incl. **70.9% of protein-labelled
voxels** (belt 90.2%). **Why:** at 8 A a C4 tetramer's rings are dominated by m=0 (the helices merge into rings), so the
fraction of power above m=2 is small in the protein too; the outer belt's ratio rises from noise over weak density.
Order ([[cnga1-local-order-separates-micelle]]) stays the only working label.

Blind recovery (box 96, 15 A): tilt landscape 0:0.069 5:0.074 10:0.077 20:0.083 45:0.135 90:0.174 -> tilt found.
Centre landscape along x is flat/symmetric (+-10 A: 0.0752 vs 0: 0.0757), refined offset (-9,-5) A = wrong. **Fix if
revisited:** fit the centre by minimising m=1 power (remove only m=0 for the centre stage). Caveat: the map is C4-imposed,
which makes the axis easier than a true C1 case (PfCRT). A toy is needed before trusting any micelle label/prior.
