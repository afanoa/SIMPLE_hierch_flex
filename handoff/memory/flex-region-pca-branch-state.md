---
name: flex-region-pca-branch-state
description: "2026-09-20: branch flex-region-pca (Mac /Users/atifao/SIMPLE) = master 34e69f0c4 + the joint region-block fit for flex_pca (pcg_mskfile as an integer-label volume, neigs per region, per-column block windows in M-step/deflation/orthonormalisation/merge/PCG operator); BUILT on batch2 (~/SIMPLE_region, gcc 15.2, same flags as the user build) and VALIDATED on the toy: the joint block fit removes the leak (theta2/theta3 from region-1 columns 0.61/0.69 -> 0.014/0.033) at equal attribution/prediction; toy harness + RNA launchers committed under scripts/flex_pca_verify/region_toy"
metadata:
  type: project
---

Commits: 705047b96 (note + harness), then the implementation commit (8 files, +609/-51). Design: section 7.3 of
`doc/implementation_notes/flex_rna_symmetry_multibody_investigation_2026_09_19.md`.

**How it works:** `flex_pcg_support_volume` detects an integer-label pcg_mskfile (all values within 1e-3 of integers,
max >= 2) -> `region_windows` (indicator, 20 A low-pass, Fourier crop, clamp/floor) -> `flex_blk_img(1:B)`, union in
`flex_env_img`; ranges `flex_blk_lo/hi` set by `flex_set_block_ranges` at init (`init_basis_datafree`: impulses cut to
each block, orthonormalised per block, `neigs` columns PER REGION) and re-registered by `orthonormalize_blocks` (per-block
Gram eigendecomposition, then `basis_recs_from_images`). Per-column windows: `flex_window_apply_col/_rec_col`, per-region
`flex_window_apply_region`; mstep realisation + deflation loop (`nblk_dfl = flex_nblocks`), pairmerge (window per column,
`merge_deflate_mean_shaped` per block, Mba and gauge Mra block-diagonalised), `flex_pcg_t%set_window_blocks` +
`mask_work/mask_comp` per component in apply_operator/finalize_rhs/apply_precond/solve. Workers (no init) install the
union window only. Gridding backend ignores labels (needs rec_backend=pcg, the default).

**Validation recipe (no build):** `scripts/flex_pca_verify/region_toy/` is on batch2 as `~/toy_*`, `~/launch_toy.sh`;
the syntax checker is in the session scratchpad (`flex_syntax_check.sh`: devtools-gfortran -fsyntax-only in dependency
order, baseline master = 0 errors, branch = 0 errors). **Next:** build on batch2 (user's call), then the toy with
`pcg_mskfile=<labels at box 96>` (write labels from the reg96 masks: argmax label, 0 outside) vs the composed union
(leak: region-1 block must stop predicting theta2/theta3), then the RNA with `G10cov_labels_B3.mrc`.
See [[multibody-flex-means-region-pca]], [[rnasep-flex-seeds-vs-refine3d-finals]].

**2026-09-20 validation:** batch2 build `~/SIMPLE_region/build` (worktree of staging ref `flex-region-pca-in`, pushed
with `git push atifao@batch2.ncifcrf.gov:SIMPLE flex-region-pca:refs/heads/flex-region-pca-in`); launchers take
`EXE=` and `SIMPLE_PATH_OVERRIDE=` from the environment. Toy block fit `toy/T_B3`: leak gone, column energies 0.96-0.99 in
their regions, 12/12 pairs >= 0.7. Regression: branch vs user binary differ at 1e-4 in pre-EM numbers (compilation FP
drift), 1-4% after 4 EM iterations; exact gate = `T_G4master` (master built the same way in `~/SIMPLE_masterchk`).
RNA block fit `B_G10cov_B3` (labels `G10cov_labels_B3.mrc`, neigs=10 per region) launched with the branch binary.

**2026-09-20 later:** `rec_states=regions` (36bbe210a) = built-in per-region on-axis delivery (3 equal-mass states per
block, halves, `flex_pca_region<b>_state_00N.mrc`); resume recovers ranges from the manifest or an equal split
(616d2b012); `SIMPLE_COV_DECONV=0` now really skips the deconvolution (ad80cf3cd; the K ladder at rank 30 on 932k
particles takes HOURS, G30 sat in it > 3 h -- always set it for states-only reruns). Second batch2 build
`~/SIMPLE_region2` (staging ref `flex-region-pca-stage`, never checked out: push there, then `git reset --hard` in the
worktree; pushing to a checked-out ref is refused). RNA block fit `B_G10cov_B3`: 30/30 columns >= 0.85 (mean 0.995),
91-98% in-region, vs G10 10/10 at 0.96. Toy `T_B3_regions`: built-in delivery = manual reruns (0.974/0.974/0.979).
flex_pca is NOT run-to-run deterministic (same binary: eigenvalues 1141 vs 1129 on the toy) -- see
[[flex-pca-distribution-changes-results]].

**2026-09-20 afternoon:** branch head 00f20b84c (`flex-region-pca`), batch2 build `~/SIMPLE_region2/build/bin/simple_exec`
14:08 (includes the rotated-resume adoption guard 90c90a30b) (staging ref `flex-region-pca-stage`; a NEW source file needs `cmake .` first, see
[[batch2-new-source-needs-cmake-rerun]]). New on the branch: 9027314f3 sparse L1 M-step (`SIMPLE_COV_L1`, Lowdin frame,
[[flex-pca-sparse-l1-mstep-verdict]]), e5dbeadb1 varimax stage (`SIMPLE_COV_VARIMAX=1`, `simple_flex_pca_varimax.f90`,
[[flex-pca-varimax-rotation-recovers-region-axes]]); toy runs `toy/T_L1_*`, `T_LOW_r4/r12`, `T_G4vmx` (Fortran varimax
resume = prototype); RNase P `V10_ax1..10` = G10 cache resumed in the varimax frame, on-axis deliveries, scored by
`rna_axis_score.py <dir> 1.464 G10cov_B3_{1,2,3}.mrc` against C_regions/B_regions (78-80 % / 68-81 % in-region, e/o 0.98).
Scorers: `toy_score.py` (label-free attribution by energy + localisation index, prefers flex_pca_varimax_pc*),
`varimax_score.py` (post-hoc rotation prototype).
