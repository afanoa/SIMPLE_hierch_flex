---
name: cls-expansion-solvent-prior-2026-10-09
description: "cls_expansion solvent prior (3D-PCG-style real-space ridge on the basis, inscribed-circle particle mask) - user judged RNase P averages best with it; half-agreement metric rewards reproducible background and misleads"
metadata:
  node_type: memory
  type: project
  originSessionId: b8bc4431-3743-4b17-aafd-6d31e354ea19
  modified: 2026-10-09T19:58:15.690Z
---

Branch cls-expansion-solvent (~/SIMPLE_clsclean, uncommitted, off cls-expansion-clean + its uncommitted cleanup): particles
masked by the box's inscribed circle (not mskdiam: keeps beyond-mask heterogeneity, removes rotated box corners); the PPCA
M-step solved by PCG with lambda*||(1-m) U_q||^2 in real space (separable partial DFTs, origin at box/2 - the planes'
convention; getting the origin wrong made the support 7-9%); m = soft support = Otsu core of |mean| plus spread above Otsu
within 3.3x0.6 band-resolutions of the core, dilated 0.6 band res, Gaussian edge. Env: SIMPLE_FLEXCLS_SOLVENT_LAMBDA
(0 = diagnostics only, bitwise == no prior), SIMPLE_FLEXCLS_TIGHTMASK, SIMPLE_FLEXCLS_UNCAPPED (experiments).
Realistic toy (gen_toy3 TOY_NEIGHBOR=0.3 TOY_ICE=0.3 TOY_ARM=1): tight mskdiam loses the beyond-mask arm (amp 0.54) and its
main direction is background; inscribed circle keeps the arm but neighbours enter; circle+prior keeps the arm, solvent share
of deviation variance 0.37 -> 0.07, deviation error 0.76 -> 0.63.
RNase P top-5 (flex3sets/rnasep/clsexp_{mask,circle,circle_prior}): without the prior the main direction is a background
level (69% / 50% of variation, bright/dark discs); with it the main direction is structural. **The user judged the prior
best (2026-10-09).** The even/odd deviation agreement ranked it LAST (0.80-0.84 vs tight 0.86-0.90): reproducible
background inflates that metric. **How to apply:** do not rank masking/solvent variants by deviation half-agreement; use
truth (toy), the strips, or downstream 3D. Related: [[cls-expansion-fortran-toy]], [[cavg3D-states-findings-2026-10-08]].
