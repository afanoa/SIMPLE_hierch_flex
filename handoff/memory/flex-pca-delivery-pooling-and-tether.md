---
name: flex-pca-delivery-pooling-and-tether
description: "2026-09-09: state maps were 'breathing' because of per-state low-pass + noise, not a mode (dilation <= 1.6% of any state's deviation); shipped a RECOVAR-style per-shell adaptive kernel at delivery (own / macro-cluster pool / global pool by even-odd error, one common low-pass, diff maps), a focus mask volume (SIMPLE_COV_FOCUS_MSK), and the continuation tether (SIMPLE_COV_TETHER) for marching across boxes; work lives on branch flex-delivery in worktree ~/SIMPLE_wt_delivery because a parallel session is doing the file cleanup in the main checkout"
metadata: 
  node_type: memory
  type: project
  originSessionId: 7e382bce-970c-4cdb-bc32-4a1a1c1801a9
  modified: 2026-09-09T15:23:33.576Z
---

**Diagnosis (v5 focused zoom maps, breath.py on the box):** fitting each state's deviation from the
consensus with scale / dilation / Laplacian templates explains 0-19% (dilation itself <= 1.6%); Rg within
2%. The visible "breathing" was (a) each state low-passed at its OWN eo-FSC (4.4 vs 6.4 A) and (b) an
independent noise realisation per map (deviations 1.1-1.35 noise units). RECOVAR's answer is the per-shell
adaptive kernel (adaptive_kernel_discretization.py) + one local-res filter for all states; its per-image
"contrast" is grid-searched jointly with the latent and re-run once ("do-over"), and its continuum handling
is deconvolved density on <=4 PCs with a FOCUS MASK on the moving domain, never clustering.

**Shipped (branch flex-delivery, commits 08c8b1af3 + d3506191e on top of 0fcd70404; pushed to the box
as flex-pca-resident-in and built into ~/wt_res):**
- rec3D delivery: per state three nested estimates at no extra reconstruction cost (own kernel, the
  macro-cluster pool = accumulators summed, the global pool); per shell keep the one with the smallest
  error = noise (|e-o|^2/4) + bias (own-vs-pool power above its noise, own share f of the pool accounted
  for), smoothed +/-3 shells, NO monotone rule (the first version locked on one noisy crossing at ~23 A);
  one Butterworth at the worst blended res for every state; `flex_pca_state_diff_XXX.mrc` vs the global
  pool; `flex_pca_state_kernel_XXX.txt` per-shell table. `SIMPLE_COV_STATE_POOL=0` = old per-state
  low-pass. `state_group` (macro of each state by kernel mass, `state_groups_by_mass`) passed from the model.
  First result (v14, monotone rule still on): own detail only to 17-25 A, blended res 3.40 A for all.
- Focus: `SIMPLE_COV_FOCUS_MSK=<mrc>` ([0,1] on the box_crop grid or Fourier-cropped), one `focus_img`
  per fit multiplied into every basis volume (annulus builds the same image). 40S candidate mask =
  density-weighted 2-way spatial split of consensus_cdrgn_box64.mrc (`~/seg40s.py`, numpy only; the box
  has no scipy): 47/53 mass, centres 106 A apart -- a positional half, may be the 60S.
- Tether: `SIMPLE_COV_TETHER=<stage-1 flex_pca_embedding_deconv.bin>`; per fit row the stage-1 posterior
  mean/precision; iteration 1 bypasses the solve (z = m, E[zz'] = mm' + P^-1, so the first M-step re-solves
  the basis on the new grid from the old latent); later E-steps add P_i / P_i m_i to both solvers
  (optional Pt/Ptm on probe_solve_ecm / probe_solve_mix); eigvals seeded from the tether spread. Low band
  is double-counted (likelihood + tether) in this v1; the exact form restricts the E-step to the new shells.

**Why:** the old lp-schedule marching (within one box) was measured neutral; the fine-band run from
scratch (v13) dissolved the populations. Continuation across boxes is the marching that can help.
abinitio3D for reference: 8 staged refine3D runs, MAXITS=[20,20,17,15,12,12,12,25], lp ladder from class
FRCs (20 A floor -> lpstop), early stop on orientation overlap 0.9/0.95 except the NU stages, FSC=0.5
promotion at stage boundaries only.

**How to apply:** runs queued as ~/queue33.sh (v14 pooled delivery + zoom_continuum_v6_focus40s) and
~/queue34.sh (v14b with tables; v15 = paired_rankcut_br360_v15_tether96: box_crop 96 lp 12 neigs=NT
n_probe_iters=4, new cache ~/cache_ribo10028_plane96). Low CPU during a states run = the master's
single-process box-360 delivery (FFT/filter/mask per state, ~6 cores), not I/O. Do NOT edit
/Users/atifao/SIMPLE while the "Flex comment cleanup" session runs there; the box branch does not have
that cleanup, so the two lines must be rebased onto each other later.
Related: [[flex-pca-deconvolution-v8-result]], [[flex-pca-cv-ica-recipe-2026-09-01]], [[flex-pca-speed-levers-2026-09-08]].
