---
name: subclass-node-embedding-2026-10-21
description: "per-node sub-class embedding (cls_expansion sub-classes in direction bins, node flex latent pooled over members at their own poses; no sync): pooling DESTROYS the node's heterogeneity (swing R2 0.41 -> 0.09; IgM 10 arm angles 0.24-0.49 -> 0.01-0.06) because 8 sub-classes/bin are mixtures; per-view 2D latent sees no swing (0.06-0.10 vs 3D 0.43-0.45) but ~matches 3D on IgM arms; sub-class split-half test gives 3 spurious dims on ONE conformer; junk isolation works"
metadata:
  node_type: memory
  type: project
  originSessionId: 6d5934fe-4a66-44cc-9cf7-17b64c67240a
  modified: 2026-10-21T17:38:53.003Z
---

2026-10-21 (task: make the recursive tree's per-node embedding meaningful with cls_expansion sub-classes; user also suggested IgM).
Scripts routes/toy55: subcls_node3/4.sbatch (node members at the parent's poses -> dir_bins -> fresh coupled-snapshot project; members
set by `prg=selection oritype=ptcl2D infile=rows` because import IGNORES the oritab state (without it cls_expansion bins held all 400k
rows); make_cavgs; cls_expansion pose mean vol1 = node mean; DUMP=1 = SIMPLE_FLEXCLS_DUMP per-bin 2D latents; node4 COUPLED=1 = infile
flex_pca_embedding.bin -> CCA placement), subcls_embed3.py (pooled latent, split-half reliability vs within-bin shuffle, XD groups,
children, truth), subcls_reg.py (cross-fitted per-bin regression 2D latent -> 3D latent = zhat; k-means windows; within-bin truth R2),
subcls_diag.py, subcls_nullcheck.py, subcls_post.py (children view balance), subcls_rootprep.sbatch (root of a synthetic set in rflex4
layout), subcls_oracle.sbatch, subcls_maps.sbatch + subcls_mapscore.py. Outputs rnasep_klem/subcls_*, rflex_root_<set>.

Truth readouts (within-bin held-out R2 per particle: 2D latent vs node 3D flex latent; sub-class partition eta2 = ceiling of pooling):
- r22 swing node (rflex_s30a): 2D 0.105 vs 3D 0.451; sub-class eta2 0.096; pooled R2 0.086 vs per particle 0.405; children 0.06-0.09
  (rflex refined 0.406, coord leaves 0.531). ORACLE swing-only members: 2D 0.059 vs 3D 0.429 (composition is not what hides it: the
  per-view PPCA of ~1300 particles cannot resolve the swing at this SNR). k-means windows on z2 no better; zhat adds nothing.
- IgM armsigm1 root (flex neigs 16, 150 bins x 8): 2D 0.23-0.58 vs 3D 0.38-0.59 per arm angle (2D better on a3_ph, a4_ph, a5_ph);
  sub-class eta2 0.02-0.07; pooled R2 0.01-0.06 vs per particle 0.24-0.49 (ALL 10 angles, both per arm: a plain flex latent encodes
  them linearly; boards had the second angle <= 0.22); z3+zhat +0.01-0.03 over z3.
- Composition is what sub-class MEANS carry: cross-half R2 of the pooled latent explained by conformer shares 0.50 (r22) / 0.25 (r21),
  nuisances ~0; conformer NMI of the partition only 0.05 (r21).
- Validity: null_50k with truth-junk removed (ONE conformer): bin-centred sub-class means show 3 "reproducible" dims (0.62/0.48/0.40 vs
  null 95th 0.31) = artefact (shared image noise between 2D assignment and 3D latent, or in-bin nuisance); zhat 8 dims. The within-bin
  shuffle null is invalid; a stop test needs disjoint information (e.g. assignment <= 20 A, measurement from a disjoint band).
- two_50k (49 + closed 18 + junk): closed-vs-open within-bin R2 2D 0.198 vs 3D 0.466; sub-class eta2 0.06 (k-means windows on z2
  0.155); pooled 0.04 vs per particle 0.436; partition conformer NMI 0.019; 8 groups, main ones 62/29 and 75/21 % (49/18).
- Maps r22 (subcls_maps, reconstruct3D at the parent's poses, end-children difference per band 117-30/30-20/20-12/12-8 A): sub-class
  children repro 0.11/0.02/0.02/0.00, vs truth 0.21/0.09/0.01/0.00 (ends 19.5/29.3 A, view-starved); rflex coord leaves repro
  0.93/0.86/0.35/0.05, vs truth 0.87/0.85/0.51/0.15.
- Cheap PER-PARTICLE comparator (no sub-classes needed): PC1 of the bin-centred node flex latent |corr| swing 0.636 (= the sub-class
  reliability direction 0.620); view-balanced within-bin quantile bins eta2 0.371 with view KL 0.000 (width bins 0.382), ends 7.4-22.1 A
  = close to the focused refinement (0.406) for the cost of the node flex alone. MAPS of those 10 view-balanced children (12.8k each,
  6.7-7.3 A, reconstruct3D at the parent's poses, maps_pp_binq): end difference repro 0.98/0.94/0.38/0.06, vs truth 0.80/0.88/0.60/0.23
  = as good as or better than the coord leaves below 30 A (0.87/0.85/0.51/0.15) without focused refine (18 min) or score loop (32 min):
  the useful 2D ingredient is VIEW STRATIFICATION (direction bins for latent centring + balanced children), not sub-class pooling.
- CHEAP CHILDREN ON EVERY s30a NODE (cheapkids.sbatch/.py -> rnasep_klem/cheapkids_s30a_<node>; flex latent centred per direction
  bin, PC1, 10 within-bin quantile bins, reconstruct3D at the parent's poses, ~7-10 min @24 cores): labels match the node's flex states
  (NMI r 0.19/flex 0.22, r1 0.21/0.21, r2 0.22/0.25, r12 0.16/0.16, r21 0.09/0.07, r22 swing eta2 0.39/0.31) but the FOCUSED REFINEMENT
  separates more at compositional/mixed nodes (NMI r 0.45, r1 0.30, r2 0.32 + swing 0.38, r12 0.34, r21 0.15); only at the pure swing
  node r22 do they match it (0.39 vs 0.41). End-children difference maps vs truth (117-30/30-20/20-12/12-8 A) are as good or better than
  the coord leaves (r12 0.92/0.91/0.66/0.31 vs 0.86/0.76/0.43/0.15; r21 0.75/0.86/0.51/0.21 vs 0.72/0.72/0.36/0.12; r22 0.86/0.89/0.61/0.24
  vs 0.87/0.85/0.51/0.15) EXCEPT r11 (no resolvable heterogeneity): ends reproduce 0.97 but match truth 0.12 = a nuisance PC1, flat PC
  spectrum (PC1 share 0.17 vs 0.30-0.48 elsewhere) = a candidate truth-free stop flag.
- FOCUS MASKS (user: PC-derived masks noisy, use real-space variation? nodemask_check.py): focus_mask.py already uses the real-space
  cross-half variance of the flex STATE maps; corr with the node's true variance: states 0.67-0.86, cheap children 0.74-0.86, PC maps
  (sum lambda PC^2) only 0.28-0.59 (r11: 0.21/0.00/0.28). Masks: current 8-50 % of the sphere (26-50 % at r1/r12/r21/r22, recall 95 %);
  focus_mask.py on the cheap children's half maps 4-15 % (r22 4.1 % vs truth support 3.8 %, recall 80-89 %, enrichment 6-20x vs 2-3.6x).
- MASK A/B (maskab.sbatch -> rnasep_klem/maskab_s30a_<node>_kids_mask; rflex3's refinement repeated exactly, only pcg_mskfile = the
  cheap-children mask; subcls_statescore.py = per-state deviation fidelity vs Fourier-cropped truth): r22 (mask 41 -> 4 %) swing eta2
  0.41 -> 0.44, extent 8.8-22.4 -> 8.0-22.8 A, state deviation fidelity 0.63/0.76/0.49/0.21 -> 0.73/0.78/0.65/0.45, deviation
  split-half 0.98/0.93/0.59/0.11 -> 1.00/0.99/0.93/0.53, median FSC 5.7 -> 5.4 A; r12 (32 -> 8 %) NMI 0.34 -> 0.36, fidelity
  0.87/0.85/0.65/0.36 -> 0.87/0.83/0.72/0.49, split-half 0.99/0.95/0.67/0.18 -> 1.00/0.99/0.90/0.47. Labellings move a lot between
  the runs (agreement NMI 0.28 r22, 0.54 r12). REPLICATE of the current-mask run (maskab_s30a_<node>_mask, maskab_cmp.sbatch): 97.1 %
  (r22) / 97.8 % (r12) of particles in the same state, NMI 0.91 / 0.93, truth metrics equal to 0.001, map fidelity equal to 2 decimals
  -> the focused refinement is nearly deterministic and ALL mask differences (labels and maps) are real effects of the mask.
- IgM ICA (igm_ica.py): FastICA of the bin-centred root latent = PCA (matched |corr| per arm angle mean 0.31 vs PCA 0.32, supervised
  ceiling ~0.59); components are super-Gaussian (+0.5..+13 kurtosis: outliers/nuisance) while uniform angles diluted by noise sit at
  ~-0.2; truth-free junk removal (largest sub-class group) changes nothing. Variants (igm_ica3.py): SUB-GAUSSIAN PURSUIT (deflationary
  kurtosis minimisation) best, mean 0.33 truth-free / 0.37 oracle-junk vs PCA 0.29-0.33, FastICA cube 0.30-0.33; a few angles near
  ceiling (a4_th 0.64-0.65, a1_ph 0.54-0.58), others lost (a3_ph 0.08-0.14): no clean linear unmixing; the factorised boards
  (dominant angle 0.60-0.66 for all five arms) remain better. Python stdout through `| grep` is block-buffered: use -u + grep
  --line-buffered or results appear only at the end (a cancelled job loses them).
- COUPLED placement (cls_expansion infile = node flex_pca_embedding.bin, CCA of class vs 3D latent, canonical corr 0.6-0.87): r22
  sub-class swing eta2 0.096 -> 0.153, pooled R2 0.086 -> 0.142, children 0.085 -> 0.14; cls_expansion 14 min instead of 7.
- Junk works at the sub-class level: null_50k 6 groups 93-100 % junk, main group 93 % pure (80 % overall); IgM 86 % of junk in
  junk-majority sub-classes, junk groups 100 %.
- Children from divisive sub-classes: view KL 0.5-2.2 vs 0.004 random (one leaf holds ~half a bin: 772/360/221/112/95/15/7/2).
- flex_pca's own particle XD picks K=2 (r22, r21), 12 (null_50k), 15 (IgM): not a group count either.
- Cost: cls_expansion 7 min @32 (128k), 11.4 min (200k); import+make_cavgs 1.2-1.5 min; embedding 1.4 min; node flex 5.7 min @64
  (r22), 2.5 min (50k), 12.7 min @48 (200k, neigs 16); particle route r22 focused refine 18.4 min @64 + coordinate stage 32 min @64.

**Why:** sub-class AVERAGES are high-SNR but membership is a per-particle decision by a per-view model; at realistic SNR the 3D latent
is the stronger (or equal) per-particle feature, and 8 leaves per bin cannot represent a continuum or a product space.
**How to apply:** use sub-classes for junk / composition and for 2D inspection, not as pooling units for continuous or multi-domain
states; the 2D -> 3D translation that carries information is the per-bin regression of the 2D latent onto the 3D latent (cross-fitted);
any reproducibility test on sub-class means needs a calibrated null. Next for IgM: truth-free unmixing (ICA) of the root flex latent.
Related: [[routes-2d-to-3d-states-tests-2026-10-14]], [[direction-bin-route-synthetic-2026-10-18]], [[igm-like-synthetic-armsigm1-2026-10-19]],
[[klem-integrated-hierarchical-model-vision]].
