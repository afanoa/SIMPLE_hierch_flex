---
name: flex-pca-complement-block-and-step0
description: "2026-09-11: on branch flex-pca-pcg, step 0 (smpd_target working box, default 2.2 A) and the complement block (SIMPLE_COV_COMPLEMENT=k_out, two supports in one fit) are implemented and tested on PfCRT; the leading PC of every unmasked/envelope run is a one-signed MASS mode (belt thickness outside, density scale inside), not a difference volume; the 4-shell + flat-in-window deflation removes it; a PCG-alignment plan doc exists; production PCG budget is 2 iterations from ZERO (warm starts ditched 2026-09-10)"
metadata: 
  node_type: memory
  type: project
  originSessionId: 254a27e2-0d20-447a-8523-54dc0a6f3dfb
  modified: 2026-09-11T15:59:49.463Z
---

**Branch:** `flex-pca-pcg` (Mac worktree `/Users/atifao/SIMPLE_pcg`, batch2 `~/SIMPLE_pcg` via `flex-pca-pcg-in`, snapshots
`~/SIMPLE_pcg/snap_<rev>` for running jobs). The main checkout `/Users/atifao/SIMPLE` (flex-pca-resident, dirty) and the
worktree `~/SIMPLE_wt_focus2d` (branch flex-focus-2d, one snapshot commit, 2D residual masking ABANDONED after Hans's note) are
separate. `~/SIMPLE_flex` on batch2 = the old focus-mask build (env SIMPLE_COV_FOCUS_MSKFILE), superseded by `pcg_mskfile`.

**Step 0 (258556429, 90ddb71fc):** `smpd_target` primary (default 2.2 A), explicit `lp` -> smpd_target=lp/2.5, explicit
`box_crop` = override; magic-box autoscale, floor 64. Embedding cache v4 stamps box_crop/smpd_crop (module state
`cov_box_crop_glob`) and refuses a resume at another lattice. PfCRT (box 300 @ 1.035 A): default -> box 140 (5.5 A band,
14 min per probe iteration, ~80+ min per paired run, master RSS only 10.6 GB); smpd_target=3.2 -> box 100 (7.8 A band,
6-18 min per single-fit run); 2.4 -> box 128. Box 140 memory is NOT a problem; time is.

**Complement block (114a2eca4, 801728976):** `SIMPLE_COV_COMPLEMENT=k_out` with `rec_backend=pcg pcg_mskfile=<env>`; block 1 =
envelope window, block 2 = sphere minus envelope (`flex_cenv_img`); block-aware `flex_window_apply_blk` at every tail site,
block-wise deflation and orthonormalisation (`orthonormalize_blocks`), two-block data-free init, `fit%ncomp_in`,
`cov_ncomp_in_init` carried in the probe meta file. Single-fit and `maxits_pcg=0` only (paired merge and the PCG operator are not
block-aware yet). PfCRT runs 51 (ctrl) / 52 (16+4) / 54 (4+4, EM_DEFLATE=4 + DEFLATE_BG=1): supports separate exactly
(inside 93-99.9% energy in envelope, outside 0.1-1.4%); the belt mass mode lands in the outside block (pc17: one-signed, Gamma
21245, latent std 162 vs 6-9); the INSIDE pc1 of ctrl and comp is the same in-envelope mass mode (cos 0.993, pos/neg 5, R^2 0.75
on the outside block = the same per-particle density-scale nuisance) and the 4-shell + BG deflation removes it (54: all
components balanced, block coupling max corr 0.04, flat Gamma 222-522, flat deconvolved spectrum, one UMAP cloud).

**Standing diagnosis on PfCRT:** every unmasked run's pc1 (27 box 64, 50 box 140) is a one-signed belt blob (pos/neg energy
95-207), "looks like the molecule" = mass mode; check pos/neg energy and |mean/std| inside the support before believing a PC.
After envelope + complement + shell/BG deflation there is still no dominant reproducible inside axis at 7.8 A with 16.8k
particles; the next decisive test is the pcafit=1|2 half harness on the run-54 config (single-fit runs measure nothing about
reproducibility).

**PCG defaults / alignment:** the overnight session's flex M-step warm-starts from the gridding solution and defaulted to 20
iterations; production (`c9808eff9`, policy 2026-09-10) starts from ZERO with `maxits_pcg=2, rtol=0` (warn >8, refuse >100,
cold finals max(maxits,5)). Defaults changed to 2/0 on 9e9c4b091 (unbuilt); the full alignment plan (start from zero, shared
ridge rule, P_tau replay inside the operator, shared stop rules, per-component support) is
`doc/implementation_notes/flex_pca_pcg_alignment_plan.md` (6550dc57f). The user wants parity with production before trusting
any PCG M-step result. See [[flex-pca-pcg-branch-state]], [[flex-pca-real-space-priors-proposal]].

**FIRST POSITIVE ON PfCRT (2026-09-11 afternoon, runs 56/58/59/63):** protein-only support (density-threshold mask, the
100 kDa choice is arbitrary -- use the construct mass or an atomic-model mask; the density histogram inside the envelope is
NOT bimodal, Otsu lands at 195 kDa-equiv) + 4 inside + 4 complement components + shell/BG deflation, box 100, single fit:
inside components reproduce across the pcafit=1|2 halves at best-|cos| 0.94/0.74/0.75/0.78 (null 0.005); gridding states
resumed from the embedding with `nkern=4` (placement on the inside block only; log line "state stage restricted to the
leading 4 of 8 latent components") give 4 states x ~4200 particles at 6.3-6.9 A whose deviations from the ensemble mean
reproduce even-vs-odd INSIDE THE PROTEIN (same-state 0.40-0.47 vs different-state -0.14) and NOT in the micelle (0.03-0.06).
State maps are in /mnt/beegfs/elmlund/pfcrt/63_flex_pca and ~/Downloads/pfcrt_63_states. The automask support (runs 54/55)
is protein+micelle (33% of the mskdiam sphere; inside components had 51-89% of their energy in micelle) -- never use it as
the transporter support. States at 20 PCG iterations were the overnight default; production floor is 5 (maxits_pcg=0 on
the fit gives 5 on a states-only resume without a rebuild).

**CORRECTIONS + LOCALISED MODE (2026-09-11 evening):** (1) the state even/odd "gap test" is CIRCULAR when placement uses the
same images: states from the noise-dominated no-Wiener basis (run 66/68) reproduce at 0.18-0.49 vs -0.08, the same pattern as
the run-63 "positive" -- never cite it; the micelle "control" is null by the support. (2) The half-set basis match SURVIVES its
init control: init-deflated 0.91/0.71/0.74/0.74 (init basis from a 0-iteration run, 67_flex_pca). (3) The leading inside
component is 0-2% rigid-body/affine/quadratic/mass/dilation, and NOT a CTF frequency modulation (DEFLATE_CTF template,
aa9eaa609, left it at cos 0.95; per-shell corr with consensus only +-0.3). (4) It is LOCALISED: top 10% of protein voxels
hold 77% of its energy, one 22 kA^3 patch 29 A off centre in the 4-5 A core, the same patch in both halves (Jaccard 0.8);
"dispersed" was a display-level artefact. pc1_energy_native.mrc (box 300) in 56_flex_pca / ~/Downloads/pfcrt_56_pc1.
(5) SIMPLE_COV_MSTEP_WIENER=0 (9870c022f) adds noise (66-85% power beyond 19 A, het-res unchanged at 19 A), not resolution.
(6) Box 140 vs 100: het-resolution 22 vs 20 A -- the band is not the limit. (7) On-axis states (state_axis=1, run 70 from
56's embedding): 4 equal-occupancy slices along z1, eo-FSC 4.25-4.31 A, end-to-end difference rms 20% of the consensus at
10 A with 24% of its energy in the 1.6%-of-mask pc1 patch; maps in ~/Downloads/pfcrt_56_ax1_states.
