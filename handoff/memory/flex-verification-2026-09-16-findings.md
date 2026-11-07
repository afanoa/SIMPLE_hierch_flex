---
name: flex-verification-2026-09-16-findings
description: what the 09-16 branch verification on batch2 established (delivery reproduces its reference exactly; pcg+mask doubles reproducible components on matched inputs; PCG state maps lose PD; the spliceosome "regression" was the 3_selection pose mismatch; PC1 slices beat the mixture on the funnel continuum)
metadata:
  type: project
---

Verification matrix (batch2 only, `~/flex_verify/`, run roots `flex_verify_20260915/` under 10028_simple
and 10180; doc: `doc/refactoring_notes/flex_pca_branch_reconciliation_2026_09_15.md` section 3):

1. **flex-delivery a94d6c790 works as is** on the cluster without the plane cache: ribosome loglik/particle
   -91.778339, K=16, kNN 0.804 vs the `bare_b64_lp16` labels, floor drops the 232/285 states, 6 maps.
2. **flex-pca-pcg is worse than the delivery on the ribosome at every setting, gridding included**
   (6 vs 9 reproducible components; kNN vs reference labels at chance 0.29, K=5): the unported state and
   deconvolution defaults, not the backend.
3. **PCG state maps fail**: two of two ribosome runs died with "flex PCG state solve lost
   positive-definiteness (cold start)" -> gridding for the state maps.
4. **The spliceosome "regression" was an input error**: every 10180 run since 09-15 used 3_selection
   (imported star poses) with the 4_abinitio3D volume -> mean-scale ratio 6.4/4.0/-1.2 (must be ~1).
   See [[spliceosome-project-is-4-abinitio3d]]. On the right project: pcg+mask 10 of 20 reproducible
   (>=0.7 both halves) vs 3 for delivery gridding and 5 for the 09-08 run; band boxing numerically inert.
5. **Placement on the corrected latent**: funnel-shaped continuum along z1, no islands. Mixture
   (state_axis=0) = 5 lobes overlapping in z1/z2, maps 9.5-10 A (FSC0.5 19-22 A); PC1 equal-occupancy
   slices (state_axis=1, nkern=3) = 8.6-8.9 A. Volumes in `S4_Ppcg_gridst_ax0/` and `_ax1/`; figure
   `doc/implementation_notes/flex_pca_10180_umap_ax0_vs_ax1_2026_09_16.png`.
6. Harness traps: delivery/merge masters write nparts*nthr into the part scripts (cpus-per-task=160)
   -> cap with SIMPLE_COV_MASTER_NTHR on pcg or keep nparts*nthr<=96; `pcg_mskfile` refused without
   rec_backend=pcg; SOLVENT_MSK not read on pcg; a states-only resume from flex_pca_embedding.bin
   re-runs sigma2 + the 20-min deconvolution ladder (resume from the deconv cache instead);
   `update_project qsys_name=local` needs SIMPLE_QSYS_PARTITION defined.

**How to apply:** trust the mean-scale table before any basis number; use PC1 slices on 10180; score
bases by the "merge matched" rows of flex_pca_paired.txt (both halves >= 0.7). See
[[flex-pca-fork-map-2026-09-15]], [[flex-tests-run-on-cluster-only]].
Addendum: user judged the corrected `state_axis=0` spliceosome volumes (S4_Ppcg_gridst_ax0) as looking perfect,
so on matched inputs both placements deliver; PC1 slices are sharper by eo-FSC. The UMAP figure is now produced
IN-ENGINE (simple_flex_pca_plot.f90 + write_rgb_jpeg in simple_jpg.f90 + hooks in model.f90, on the Mac
SIMPLE_pcg worktree, uncommitted; new file -> cmake re-run). No python/external plotting in the pipeline.
**Placement conclusion (user, 2026-09-16):** PC1 equal-mass slicing was a workaround for a weak basis; with the
masked PCG basis (10 reproducible comps) the mixture placement (state_axis=0) finds the spliceosome states
itself. Default = mixture in the nkern reproducible components; one-axis slices only as the weak-basis fallback.
**Resume trap (pcg line, found 2026-09-16):** a states-only resume in a fresh directory re-deconvolves with noise
scale 1.0 (flex_pca_noise_scale.txt not found in cwd) instead of the fit's calibrated a (0.39-0.55) -> different
latent/K than the fit (ribosome: K=16 vs 5). All `*_gridst_*`/`_ax1` resumes today and the user's sp_*.sh resumes
were affected. Fixed in the SIMPLE_pcg worktree (one-file cache block adopted; legacy path reads the file next to
infile). Until built: run resumes IN the fit directory or copy flex_pca_noise_scale.txt + the deconv cache there.
**User (2026-09-16):** the ribosome pcg+mask + gridding-state maps (R_Ppcg_gridst_ax0, 8 states at 4.1-6.2 A,
mixture placement) "look great" too. So the recipe of record for BOTH benchmarks is: pcg + solvent mask basis,
mixture placement (state_axis=0), gridding state maps with the eo-FSC filter.
