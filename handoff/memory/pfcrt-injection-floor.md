---
name: pfcrt-injection-floor
description: "2026-09-11: known-motion injection on the real PfCRT images (half-bundle rocking about an in-membrane axis, run-56 pipeline, box 100, protein support + complement 4, shell/BG deflation, budget 0): 16 deg recovered (pc1 cos 0.85, R2 0.72, AUC 0.97), 8 deg weak (R2 0.13, AUC 0.68), 4 deg nothing -- although label-split 2D class averages show the 8 deg signal in the images (corr 0.5-0.8 with expectation). The detection floor is the pipeline's, between 2.8 and 5.6 A tip displacement; tools + calibrated forward-model conventions"
metadata: 
  node_type: memory
  type: project
  originSessionId: e5141833-c145-408a-ba48-7c3a3ffe0fbd
  modified: 2026-09-11T19:52:16.206Z
---

**Setup (all under /mnt/beegfs/elmlund/pfcrt/inj, tools ~/pcg_tools/inj on batch2):** `make_delta.py theta out.mrc`
rotates the half of the consensus on one side of a soft plane through the protein-mask centroid (normal in the
membrane plane) by +-theta/2 about an in-plane axis -> delta_theta.mrc (box 300; only ~30% of its energy lies inside
the protein support because the half-space includes micelle+Fab). `reproject` at every ptcl3D pose (needs EVEN
nspace + pgrp; applies NO shifts). `inject.py calib` fits the forward model against the real particles: image ~=
s * CTF * shift(P(V)) with **ctf sign +, shift sign - (exp(+2 pi i k.x)), astig sign irrelevant, s = 2.10** (mean
per-particle corr 0.082 full band). `inject.py inject` adds label_i * 0.5 * s * CTF_i * shift_i(P_i(delta)),
label = hash(index) sign, chunked memmap write of a 6 GB stack; then `reimport_particles ctf=yes` into a project copy
and the run-56 flex script (flex_pca_inj_{4,8,16}_script -> 73/74/75_flex_pca). `inj_eval.py theta run` scores
|cos|(pc_q, windowed delta at box 100, 20 A), R2 in span(pc1..4) and (pc1..8), AUC(label | z_q) and held-out LDA AUC.
`cmp_inj8_cavgs.py`: label-split make_cavgs of the injected stack vs make_cavgs of the pdelta stack (ctf=no) --
the image-level check that the injection is real-sized (corr 0.53-0.82 at 10 A in the core, even/odd 0.2-0.6).

**Result:** theta 16 (5.6 A at helix tips): pc1 |cos| 0.85, R2 0.72, AUC(z1) 0.97. theta 8 (2.8 A): pc2 |cos| 0.35,
R2 0.13, LDA AUC 0.68. theta 4: R2 0.03, AUC 0.55. Injected state difference rms in-mask ~25/35/48% of the consensus
in-mask rms at 4/8/16 deg -- these are LARGE motions; a few-A real helix motion is far below this pipeline's floor
at 16.9k particles. het-resolution stays 19-20 A in every arm.

**Ablations at 8 deg (full-half rotation incl. micelle):** box 140 (79): R2 0.14 / AUC 0.69 (no change); Wiener off (76):
0.09 / 0.79; deflation off (77): 0.14 / 0.54; both off (78): 0.09 / 0.62; **15 probe iterations + Wiener off (80): R2 0.15 /
AUC(z1) 0.90** -- the E-step sorting improves a lot with iterations, the delivered basis shape does not.

**PROTEIN-ONLY injection (delta x protein mask, 33% of the full delta energy; run 82, soft mask, 16 deg):** R2 0.21, pc1
|cos| 0.26, AUC 0.66 -- i.e. the full-16-deg success (75) was carried by the micelle+Fab half rotating in the images
driving the per-particle sorting; the same in-protein motion alone is NOT recovered. Mechanism: the M-step is a
z-weighted residual average, so a protein-only mode that cannot drive the sorting gets a noise-dominated basis.

**Mask finding:** flex_focus_protein_native.mrc is soft: >0.95 only 14 kA^3 of the 223 kA^3 0.5-region (helices are thin,
median depth 3.7 A), >0.05 footprint 511 kA^3 (the hard PCG solve domain, 2.3x the protein, 6-8 A into the micelle).
Tight mask flex_focus_protein_tight.mrc (binary at 0.5 + 3 A cosine edge). Real data with the tight mask (run 83):
SAME pc1 as 56/72 (cos 0.86/0.94), same localisation, het-res 20.5 A -> the mask treatment is not what shapes the real pc1.

**Next (launched 2026-09-11 evening):** 8-deg ablations to attribute the floor: SIMPLE_COV_MSTEP_WIENER=0
16p with the tight mask (flex_pca_inj16p_tight_script; first attempt raced the chain and died on a missing project) and
16p with n_probe_iters=15 + Wiener off (flex_pca_inj16p_it15_script). Trap: sed "$s_OUTPUT" expands to an empty var.

**Protein-only 16 deg, closed (runs 84/85):** tight mask (84): AUC 0.67, pc1 |cos| 0.31 -- the mask does not change
recovery. **15 probe iterations + Wiener off (85): AUC(z1) 0.915, pc1 |cos| 0.40, R2 0.16, het-res 18 A** -- the
per-particle SORTING recovers a protein-only motion with more EM iterations, while the delivered eigenvolume stays a poor
picture of it. Practical recipe for real data: n_probe_iters=15, SIMPLE_COV_MSTEP_WIENER=0, read the latent and
reconstruct on-axis states (state_axis=1) instead of judging eigenvolumes. Launched: real data with that recipe
(flex_pca_real_it15_script) and on-axis states from 85's embedding to check the state difference reproduces delta_16p.

**END-TO-END PASS (2026-09-11 night):** on-axis states (state_axis=1, nkern=4, gridding at box 300) from run 85's
embedding (run 87): the end-to-end state difference inside the protein correlates with the injected protein-only
16-deg delta at +0.98 / +0.94 / +0.90 (20 / 12 / 8 A), even/odd reproducibility 0.95 / 0.92 / 0.89, fitted amplitude
1.1-1.2x the injected. So with 15 probe iterations + Wiener off the pipeline recovers a protein-only helix-scale motion
through the STATES even though the 20-A eigenvolume (cos 0.40) never looks like it. Judge by states, not eigenvolumes.

**Real data with the recipe (run 86):** basis completely different from the 5-iteration runs (cos 0.08 to 56/83), flat
spectrum (70-90), balanced pos/neg, less concentrated (top-10% 0.65); on-axis states (STATES_86) and pcafit=1|2 halves
(real_it15_h1/h2) launched to test reproducibility; score with ~/pcg_tools/inj/states_real.py <run>.

**VERDICT on the real data with the recipe (runs 86/88/89/90):** on-axis states from 86 (run 88) show an end-to-end
difference that reproduces even/odd inside the protein (0.90/0.87/0.83 at 20/12/8 A, micelle 0.58/0.31/0.14, rms 25% of
consensus at 8 A, dispersed: top-10% 0.47) -- but that test is CIRCULAR (placement by z on the same images makes any fixed
direction reproduce). The non-circular test fails: independent pcafit=1|2 halves with the same recipe (89, 90) have
best cross-half |cos| 0.09 (max 0.09), half-2 het-res 61 A with 0 of 8 components in band. So PfCRT has no covariance
mode above the floor; the injected control says a protein-only 16-deg half-bundle rocking (5.6 A tips) WOULD have been
recovered end to end. Remaining honest routes: more particles, or discrete subsets from the well-posed 53%.

**How to apply:** any claim that flex "misses" a mode on a dataset must be tested this way (known mode injected into
the real images, rotation-invariant span/AUC readout), not by eye on eigenvolumes. See [[pfcrt-poses-are-the-blocker]],
[[flex-pca-10028-finds-real-mode-not-ratcheting]] (the earlier spike harness on the box).
