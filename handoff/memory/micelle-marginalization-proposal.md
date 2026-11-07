---
name: micelle-marginalization-proposal
description: "2026-09-29 user-requested 'more elegant use of the mask': marginalise the micelle as structured noise with a GEOMETRIC fluctuation model built from the order-mask geometry (tangent templates -grad(rho_B).v in the membrane frame + amplitude modes, elastic spectrum prior, marginal-likelihood hyperparameters); first test = do learned belt blocks lie in span(templates)?"
metadata:
  type: project
---

User pushed back (2026-09-29): the hard protein support is too crude a use of all the mask work; wants something more
mathematically elegant. Proposal made: y_i = A_i(mu + U z_i + b_i) + e_i with b_i ~ N(0, Sigma_B) MARGINALISED
(noise covariance sigma^2 I + A_i Sigma_B A_i^T, Woodbury in the E-step, no M-step for the micelle block), where
Sigma_B = T Psi(theta) T^T is built from the mask-derived geometry: T = tangent templates of the mean belt density
rho_B (dRho = -grad rho_B . v, v = radial/tangential/axial x cos/sin(m phi) x P_n(h/H)) + amplitude templates
rho_B x harmonics; Psi = elastic-sheet-like spectrum (~1/(kappa q^4 + gamma q^2 + mu)), 2-3 hyperparameters by marginal
likelihood. Handles the projection-overlap leak exactly, no hard edge for the protein, parsimonious for the flat
belt spectrum, physical readout (micelle softness). Implementation path: a FIXED generated block in the region branch's
joint E-step (prior variances instead of learned columns). Precedent: normal-mode bases (HEMNMA), structured noise
generally; micelle-specific geometric marginalisation not seen.
First validation (`flex_region/tangent_test.py`, 108 templates, m<=4, n<=2): fraction of the learned belt blocks
(v2_R2_block10 cols 11-20, s_soft1 cols 11-20) inside span(T) vs a random belt-confined null of equal size, with
s_hard protein columns as the control. Related: [[agnostic-order-protein-mask]], [[cnga1-micelle-leak-is-negligible]],
[[e1-rigid-tangent-subspace]], [[flex-soft-blocks-branch]].

**TANGENT TEST RESULT (2026-09-29 01:57, job 53422477):** 108 templates capture 0.218 of the R2 belt block's eigenvalue-
weighted energy (random belt-confined null 0.100), 0.177 of s_soft1's nuisance block (null 0.074), and 0.022 of the
s_hard PROTEIN components (null 0.001) -> specific (2x null, 10x the protein control) but only ~1/5 of the learned micelle
variance. Captured part: radial deformation m=1 (shift) 0.18, m=2 (oval) 0.14, amplitude m=1 0.12 / m=2 0.11, radial
m=3 0.07, m=0 (breathing) 0.07. Caveat: the learned blocks also hold solvent noise outside the belt (region 2 = sphere
minus protein), which belt-confined templates cannot represent -- re-measure inside the belt only and with higher orders
before judging. Implication: Sigma_B = geometric low-rank templates + a belt-confined smooth random field (texture).
**REFINED TANGENT TEST (tangent_test2.py, 2026-09-29 02:5x):** only 6-39 % of each learned micelle-block column's energy
lies INSIDE the belt -- the learned "micelle" blocks are mostly solvent/background variance (another nuisance). Belt-
restricted capture, templates vs random belt-confined null of equal size: m<=2,n<=2 (60): 0.45/0.19 (R2), 0.57/0.21
(s_soft1); m<=4 (108): 0.54/0.31, 0.64/0.33; m<=8,n<=4 (340): 0.64/0.66, 0.73/0.69 (no better than random). =>
~half the in-belt micelle variance is LOW-ORDER shape change (60 templates at 2.5x chance); the rest is generic texture.
Model: Sigma_B = T_60 Psi T_60^T + belt-confined smooth random field; plus a smooth solvent/background nuisance
(constant + gradients; DEFLATE_BG covers only the flat term).
