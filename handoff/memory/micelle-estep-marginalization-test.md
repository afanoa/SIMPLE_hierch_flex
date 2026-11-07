---
name: micelle-estep-marginalization-test
description: "2026-09-29 compose-route test on cnga1: embedding a 50-column geometric micelle/solvent block jointly with the protein basis cuts the automask basis's belt leak 14.0% -> 4.7% (order mask 0.93% -> 0.75%); the remaining third is belt texture beyond the low-order templates"
metadata:
  node_type: memory
  type: project
  originSessionId: 2590d10f-7aff-4aff-b2cc-d9d85f962ad3
  modified: 2026-09-29T14:20:12.477Z
---

Test (flex-softblock binary, compose route, no refit): `flex_region/nuis_cnga1` = 50 columns from `nuis_model.py`
(29 tangent/amplitude templates of the belt density in the order-mask membrane frame, 12 belt-confined smooth random
field, 9 solvent polynomials; prior variances calibrated on v2_R2_block10's rest block, whose variance the whole 90-column
span captures only 32 %). Arms `m_{hard,auto}_{base,nuis}` (SIMPLE_COV_COMPOSE=nuis,s_X; COMPOSE_CUT=0 DECONV=0
RELPRIOR=0; ~13-19 min each, 32 cpu, nparts=1); scorer `marg_eval.py` (CV R^2 of the 10 protein latents on the R2 belt-
block latents).

Leak: hard 0.93 % -> 0.75 % (protein latents unchanged, cancorr 0.98-1.0); automask 14.0 % -> 4.7 %, and the automask
latents move toward the hard-mask ones (cancorr 0.91/.87/.85/.73/.60/.48 -> .93/.90/.88/.81/.73/.56). The nuisance
latents explain 39 % of the belt reference. Compose also DEFLATES the basis against span(T): s_auto columns keep only
0.56-0.92 of their energy (s_hard 0.98-0.995), so this is deflation + marginalisation together, not the E-step alone.

**Why it matters:** the geometric block removes ~2/3 of the micelle leak without any hard cut, but not all -- the rest is
belt texture outside 50 low-order templates. With the order mask the block is nearly redundant for leak.
**How to apply:** don't expect the geometric block alone to replace the order mask; the full EM version (Schur-fold
E-step, branch flex-nuis) is the real test. Related: [[micelle-marginalization-proposal]],
[[flex-nuisance-block-code-map]], [[agnostic-order-protein-mask]], [[cnga1-micelle-leak-is-negligible]].
