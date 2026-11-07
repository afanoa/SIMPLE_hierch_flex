---
name: relion-multibody-and-symexpand-facts
description: "RELION 5 multi-body as read from /Users/atifao/relion (f6d3ad20): bodies refined SEQUENTIALLY per particle, other bodies subtracted on the fly in Fourier space (no stacks), soft masks must be pre-made, overlaps allowed via extra projectors, per-body maps from unmasked (optionally subtracted) images, flex_analyse = mass-weighted PCA of 6 DOF/body with maps by rigid COMPOSITION (no backprojection); sigma_offset-in-pixels never stored (bug); DynaMight model is not in the clone"
metadata:
  type: reference
---

Local clone `/Users/atifao/relion`. Key file:lines: `src/ml_model.cpp:1168-1457` initialiseBodies (no softening, range
check only; COM rounded; rotate-relative-to = axis from body COM to reference body COM, `orient_bodies` aligns it to -Z;
`rlnBodySigmaAngles` sets tilt AND psi sigma, no rot sigma; offsets in A via `rlnBodySigmaOffsetAngst` -- the deprecated
pixel label is parsed but never stored -> sigma 0 -> no translational search, and the GUI example STAR uses it;
overlap projectors `:1408-1455`; `--multibody_norm_overlap` loop only visits the first X line `:1374-1377`).
`src/ml_optimiser.cpp:4433-4715` expectationOneParticle: `for ibody` sequential (bodies < i see this iteration's
estimates, > i last iteration's); `:6533-6742` Fsum_obody: other bodies projected at
`Aori * orient^T * A_rot90 * Aresi * orient`, summed, CTF on the sum, soft-masked and subtracted from the masked
alignment image ALWAYS, from the unmasked backprojection image only with `--reconstruct_subtracted_bodies` (GUI
default yes). Priors: local search around previous body angles nested with Gaussian on deviation from (0,90,0),
3 sigma cutoff (`healpix_sampling.cpp:882-924`); zero-centred offset prior (`ml_optimiser.cpp:7736-7749`); norm/scale
corrections off. Reconstruction: one BPref per body, body mask re-applied to the projector each iteration
(`ml_model.cpp:1523-1553`), never to the stored map; per-body phase-randomised FSC with the body mask
(`ml_optimiser_mpi.cpp:3221-3400`); converged when all bodies stall. Per-body orientations live in `images_body_N`
data blocks with ordinary labels, identity residual = (0,90,0) (`exp_model.cpp:463-499`).
`src/flex_analyser.cpp`: PCA over 6 numbers/body (Euler degrees + 3D-lifted offset) weighted by the L2 density change
of a 1-degree rotation / 1-px shift (`:99-168`); maps along eigenvectors = rigid composition of body maps averaged
where masks overlap (`:554-654`); `--select_eigenvalue_min/max` particle subsets. `src/apps/particle_symmetry_expand.cpp:171-189`:
one row per operator, Euler angles recomposed, everything else (incl. rlnRandomSubset) copied verbatim; not a pipeline
job; no documented multibody interaction. DynaMight: only `scripts/python_dynamight.in` + job options
(`pipeline_jobs.cpp:5538-5682`: 10,000 Gaussians default = 1-2 per residue, regularization factor, half-set 0 =
validation set for deformation error, inverse deformations 200 epochs, deformed backprojection, mask option commented
out, no latent-dim option); the model is the external `dynamight` package.
