---
name: micelle-linear-leak-toy
description: "2026-09-24: the cheap micelle toy is a LINEAR noise-free least-squares test (lin_leak.py, K=<CTF^2>/|k| metric, box 96, minutes on 8 cores) not simulated particles; real cnga1 belt variance is 63% coarser than 30 A and only 6% at 15-8 A, so belt 'texture' must be ~0 (0.6 was 5x the smooth energy, unrealistic)"
metadata:
  node_type: memory
  type: project
  originSessionId: 2590d10f-7aff-4aff-b2cc-d9d85f962ad3
  modified: 2026-09-24T15:11:32.167Z
---

Dir `/mnt/beegfs/elmlund/afan/cnga1/toy_micelle/` (scripts also in the session scratchpad): `lin_leak.py <TEX>` builds 120
belt fluctuations on cnga1 geometry (breathing, oval, m3-6 ripples, thickness, shift, texture*TEX), confined to an eroded
belt (0.0000 % on protein voxels), solves the protein-window least-squares fit of each in the image metric
K = <CTF^2>/|k| (8 A band) jointly with a K-PCA belt block of r_b columns, and reports leaked/belt energy, TM-edge share
and whether the top-10 protein columns keep injected modes A (TM-edge bulge) and B (cytoplasmic twist).
**The particle-simulation toy (build_toy.py / run_sim.sh) was abandoned:** belt 82.6 % of variance and 13.9 % of the null
variance sat on protein voxels through the soft split -- no clean truth. Use the linear test first; real-flex toy only after.

Calibration (`spec.py`): real sweep_n10 belt variance band split >60 / 60-30 / 30-15 / 15-8 A = 0.27/0.37/0.30/0.06
(protein 0.25/0.38/0.24/0.13; FSC-shaped basis, so the real fine share is a lower bound). Toy TEX=0.6: 0.04/0.11/0.36/0.49
(texture/smooth energy 4.98); smooth-only 0.17/0.28/0.38/0.17, eigen share 10 modes 0.88, ev1/ev10 22.7. Real whole-basis
ev1/ev10 1.35 -- its flatness is noise, not detergent texture. TEX=0.6 result (unrealistic, keep only as the extreme):
r_b=0 leaks 20 % of belt K-energy, 40 % of it at the TM edge (28 % of protein voxels), yet top-10 recovery of A/B 1.00
(B 0.96 at belt/protein 10). Reruns at TEX=0 and 0.134 (texture/smooth 0.25): jobs 53343856/7.
Related: [[cnga1-protein-cut-leaks-belt-variance]], [[cnga1-local-order-separates-micelle]].

**RESULTS (TEX 0, belt split 0.19/0.30/0.36/0.14, K-spectrum 10 modes 0.88 / 20 0.98):** leaked/belt K-energy by belt-block
rank r_b 0/2/5/10/20/40 = **0.204 / 0.123 / 0.078 / 0.042 / 0.011 / 0.001**; 37-40 % of every leak at the TM edge (28 % of
protein voxels, ~1.35x). About 20-50 % of whatever belt energy the block leaves out leaks (finer residual modes leak more).
TEX 0.134: 0.203 / 0.133 / 0.093 / 0.061 / 0.033 / 0.018. **Recovery of injected A/B in the top-10 protein columns = 1.00
in every case** (belt/protein variance 0.5-10, even r_b=0): in the ideal LS limit the leak is a spread of weak spurious
TM-edge directions (false positives), not a displacement of real modes. L2 leak energy is 3-10x the belt's (ill-conditioned
high-frequency protein-voxel junk the real Wiener/FSC-regularised M-step would damp) -- judge by K. Not tested here: noise,
finite data, E-step coupling, deflation -> needs a small real-flex toy. Belt block of ~20 columns is what TEX-0 needs for
leak ~1 %; the region branch gives neigs per region uniformly (10+10), so 20 belt columns needs a per-region rank.
