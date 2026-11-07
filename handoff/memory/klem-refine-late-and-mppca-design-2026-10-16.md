---
name: klem-refine-late-and-mppca-design-2026-10-16
description: toy refine3D_states matrix (2026-10-16): every local refinement against the 10 k-means mixture references worsened poses (1.41 -> 2.8 deg median, p90 5 deg, baseline included), FSC unchanged, frequency march kept identification (8/18 vs 5/18); policy = grow the tree at consensus poses, refine only supported leaf states (pure references); real A/B f692c (with) vs f692d (without). Agreed design for the one-paper flex_pca integration: mixture of PPCAs on the lattice (K=1 -> flex_pca, q=0 -> klem)
metadata:
  type: project
---

Refinement matrix (toyah fam1, 87k, 10 seeds, true poses known; job 53911367): start 1.41/2.06 deg; baseline local
20->8 maxits 6 (65 min, 38 its) 2.79/5.00 deg, 8/18 identified at 0.50; no march 10->8 (32 min, 22 its) 2.81/5.03,
5/18; single band 8 A (29 min) 2.82/5.05, 5/18; pose_policy=fixed unsupported (local|global only); coarse sampling 1000/100 (43 min, 22 its) 2.81/5.04, 5/18; small sample nsample 20000 (58 min, 35 its) 3.39/6.59, 8/18 at 0.42: matrix complete 21:25, every setting worsens the poses, none changes the FSC. FSC 10.18/9.76 A in every row. Matches [[toyl-direction-bin-scheme-results-2026-10-14]] (mixture
references worsen poses, pure references help). refine3D_states also segfaults on toyah fam2 (8 min in, twice).
Policy to test: no family-level refinement; refine supported leaves against their own maps; rebuild their level 2.
Real-data A/B launched 19:5x: f692c (chain8, with refinement, job 53911867) vs f692d (chain8_sr, no refinement,
NPARTS 20, PARMAX 4, job 53912605), both CTFABS=1 PSCALE=1.

MPPCA design (user wants flex_pca's per-image pose/CTF model inside the method, one paper): y_i = C_i P_i (V_k + B z_i)
+ noise; cells V_k free (klem), one small basis B per node fitted by flex_pca's PPCA solver on the responsibility-
weighted per-image residuals; E-step responsibilities = evidence marginalised over z (Woodbury term = FLEX2D_CLS_MARG
formula in 3D); readout = posterior mean of z per image. Guards against the 2D failure (basis absorbing the state):
fit on residuals after the cells, deflate B of the cell-deviation span, rank <= 2, keep only cross-half reproducible
directions; pose tangents in B = linearised pose marginalisation. Limits: K=1 -> flex_pca, q=0 -> klem. Cost: E-step
x(1+q). Prototype in klem.py (pose_marg machinery), toy10 must keep 10/10. Cell-PPCA (simple_ppca on the cells, BIC
rank) implemented in klem_tree as rank ppca (test job 53912578) = the q-from-cells special case and fallback.
Related: [[klem-f692-fortran-chain]], [[class-marginal-evidence-to-3d-not-worth-it]], [[klem-priority-every-state-before-refinement]].

Prototype 20:00-20:15: klem.py --node_basis q [--node_basis_after 3 --node_basis_marg 1 --node_basis_rep 0.5
--node_basis_deflate 1]: E-step (numpy path only) projects the other half's basis once per image (q slices), per cell
beta = 2 Re(Q^H W r), M = I + 2 Re(Q^H W Q) (klem's ll has no 1/2), evidence += 1/2 beta' M^-1 beta - 1/2 log det M,
ZPOST (n, K, q); zb = responsibility-weighted posterior mean written as zb columns of klem_latent.txt; M-step =
weighted backprojections with weights r_ik z_ikj (K q columns) and r_ik z_ikj^2 (q columns) through _mstep_chunk /
fort_mstep, B_j = sum_k (num_kj - den_kj V_k) / (den_j + eps) per half, deflated of the half's cell-deviation span,
Gram-Schmidt, sign-aligned across halves, zeroed if cross-half corr < 0.5. First attempt died: init from the cell
PCA lies in the cell span -> deflation annihilated it (amplitude 2e-16), halves' signs opposite (-0.58), and the
block used G['REF'] before assignment (IndexError). Fixed: local REF, init iteration without deflation and without
the evidence term, sign alignment. Toy10 A/B jobs 53912835 (q1), 53912836 (q1 readout only), 53912837 (q2) vs base2
(10/10, 15-NN 0.387, closed AUC 0.870, purity mean 0.51). QOS csbdevel caps the user at 642 cores (afm jobs count).

Cell-PPCA in the Fortran tree (rank ppca, simple_ppca on the 49 mass-weighted cell deviations, job 53912578, 20:28): axis |corr| 0.98 with the Gram-matrix axis (tree_rank_pc1), i.e. the same chart; BIC rank saturates at the largest candidate (N = 49 samples, D ~ 4e5, sigma2 6e-6) so rank selection must use cross-half reproducibility of components, not BIC. Sequential board path validated with the new build (toy 20k, 46 s).

MPPCA toy10 result (jobs 53912835-37, 20:40): the init basis (cell PC1, cross-half 0.98, amplitude 0.64 of the RMS cell deviation) gives a per-image coordinate with Spearman 0.74 vs the truth axis (lattice z2: 0.81) and closed AUC-dist 0.20-0.26; the residual-fitted, deflated basis does NOT reproduce between halves (0.18, 0.08) and is zeroed -> no evidence term ever acts -> results identical to the baseline (10/10, 15-NN 0.387, AUC 0.870, purity 0.51). Correct behaviour: toy10 has no within-conformer motion (10 discrete maps), so the residual basis is noise and the cross-half guard rejects it. The benefit test needs a toy with continuous within-node motion (dense intermediates along a deformation path).

Real data, 100k random subset of the 692k, Python parent board (|CTF|, pscale; jobs 53913272/73, 21:00): the init basis (cell PC1) reproduces 0.81 across halves, but the residual-fitted deflated component reproduces only 0.24 at 10 A and is zeroed -> no reproducible within-node residual mode at the parent level; results identical to the baseline (9/10 at 0.36, 18/55 at 0.44, cells 49/49 at median 0.74). Follow-ups submitted: no deflation at the parent (TAG nodefl) and the family-1 region board at refined poses, 12->8 A, focus region_1 (TAG fam1reg).

Follow-ups 21:35: parent WITHOUT deflation: the residual basis reproduces 0.57-0.64 and is kept (amplitude 0.73-0.79), correlates with the lattice coordinates (0.37-0.40 with z1) = largely the between-cell direction; with the marginalised evidence active the results are unchanged (15-NN 0.442 vs 0.443, purity 0.469 vs 0.472, 9/10 at 0.38, 18/55 at 0.45): no harm, no gain. Family-1 region board (refined poses, 12->8 A, focus regions_child/region_1): deflated residual basis 0 at every fit; the board itself identifies 7/10 at 0.45 and 20/55 at 0.55 (cells 49/49 at median 0.95).
Caveat: the fam1reg run used meta_all (consensus poses), not meta_fb: the script's own META assignment overrode the env (fixed in run_nb_real.sbatch afterwards). Its init basis reproduced 0.72/0.58, the deflated residual fits 0.05/0 -> same verdict at the region level: no reproducible residual mode after the cells. VERDICT: the per-image PPCA has nothing reproducible to carry on RNase P at any level tested; marginalised evidence over a kept (undeflated) basis is harmless (no purity change); flex_pca's per-image part = the chart readout (per-image coordinate along the cell chart).
Refined-pose region test (meta_fb, regions/region_1, 12->8 A, job 53913835, 23:00): deflated residual basis 0 at every fit -> verdict final at every level; the region board at refined poses identifies 7/10 at 0.48 and 19/55 at 0.57 (consensus poses: 7/10 at 0.45, 20/55 at 0.55), cells 49/49 at median 0.98.

LANDSCAPE-LEVEL A/B ON THE TOY (both chains finished, 23:05): refinement at the family node HELPS despite the worse global pose metric: refined chain 22/48/90/173 states, 6-7 of the 10 apo conformers identified, family-1 leaf purity 0.21-0.29 (6 states, 5-6 conformers covered each); no-refinement chain 14/32/79/142 states, 3 identified, family-1 leaf purity 0.12-0.16. Family 2 (unrefined in both after the segfault) identical. So the pose-error metric vs the simulation truth is not the relevant quantity; keep family-level refinement in production until the real A/B (f692c vs f692d) says otherwise; the 'leaf-only' policy is NOT supported by the toy.

**Confirmed 2026-10-17 (late evening):** after the re-explanation (flex_pca = linear decoder, projection/CTF-aware
embedding; klem board = same E-step with a kernel decoder on the latent side; flex latent only initialises the parent
board), the user said "I would put it in one paper": one model, two limits (lattice of one cell -> flex_pca; basis rank
zero -> klem), with the per-node PPCA basis (mixture of PPCAs on the lattice) as the integration. Still Python-only
(node_basis); the Fortran chain runs boards only.
