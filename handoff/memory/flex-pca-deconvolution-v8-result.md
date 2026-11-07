---
name: flex-pca-deconvolution-v8-result
description: "2026-09-08 ribosome (10028, full N) latent deconvolution outcome: v8 turned the UMAP blob into a shaped manifold (d' 3.09->3.24, purity 0.54->0.68); v9's 'seed on worst-explained particles' caught noise outliers and LOST the compact populations (revert); v9b (quantile seeding, mixture-defined macro-clusters, min_neff floor) = 5 macro-clusters, island d' 3.71 / purity 0.70 / detached 0.57, best so far; the held-out K hit the ceiling 8 -> raised to 16; deconvolution took 900-1600 s at full N until the K ladder moved to a 20k subsample"
metadata:
  node_type: memory
  type: project
  originSessionId: 7e382bce-970c-4cdb-bc32-4a1a1c1801a9
  modified: 2026-09-08T23:13:32.051Z
---

Runs on the box (branch flex-pca-resident), delivered under `~/Downloads/flex_paired_20260908/`:
`ribo10028_v8_deconv` (maps+halves, embedding, coordinates, prior, UMAPs, figures), `ribo10028_v9`,
`ribo10028_v9b_states` (states-only rerun of the v9 embedding: figures + maps). Figures per run:
umap_raw_vs_deconv.png, umap_deconv_by_state.png, umap_deconv_by_component.png (scripts ~/umap_*.py on
the box; axes clipped to the 0.5-99.5 percentiles because wide-component particles sit far out).

**Pipeline as of the evening** (`simple_flex_pca_deconv.f90`, hooked in `apply_latent_deconvolution`
at BOTH delivery points of the model -- the production recipe returns at the probe-only stop): noise
scale a from the even/odd half solutions (2.6 at full N), extreme-deconvolution mixture with
projection R_i = A^-1 D and noise a A^-1 D A^-1, K by particle-half held-out likelihood on a 20k
strided subsample (final fit on all), quantile init at every rung, K ceiling 16, 150 iterations, tol
1e-6; posterior means replace z, posterior precisions replace latent_second; components -> macro-
clusters of GMM AUTO (`macro_in`), floor min_neff, over-budget clusters fold smallest-first; the paired
stage caches `flex_pca_embedding_deconv.bin` + `flex_pca_deconv_labels.txt` and a resume adopts them
(no re-deconvolution, no UMAP); `SIMPLE_COV_STATE_REC=0` skips map reconstruction for tests; UMAP over
every particle (was a 20k subsample; ~55 s per call at 105k).

**What was learned.** (1) Seeding a new component on the 5% worst-explained particles catches noise
outliers (five of eight components tiny and 20x wider than the body) and destroys the compact
populations: v9 d'/purity fine but the island became an outlier spray. Quantile init is right. (2) The
5000-particle floor (GMM_MIN_OCC) folds a 3-4k island clump; mixture-defined clusters use min_neff.
(3) The two island clumps are 4-6 sigma apart in the 17-D latent (clump_sep.py) yet stayed one
component at K=8 = the ceiling -> ceiling 16 (queue18 tests it on the v10 embedding, no maps,
npreimages=8). (4) v9b state table: 65k continuum in 2 sections (3.68 A), 20.6k lobe (4.05 A), 6.1k
(5.4 A), 5.7k (7.1 A), 7.5k island (5.1 A). (5) The principal-angle rank gate agrees across
distribution (19 vs 20) but prunes nothing at 0.5; needs a null. (6) v9 anatomy at full N with the
plane cache: paired EM 6 x 45-58 s, joint 34 s, embed 24 s, UMAP ~55 s, deconvolution 1178 s (fixed);
states round 156 s + reduce 52 s; a full run should land near 17 min after the fixes (queue17 = v10
measures it).

Related: [[flex-pca-continuum-kernel-fix]], [[flex-pca-distribution-changes-results]],
[[flex-pca-cache-taper-mismatch]], [[flex-pca-speed-levers-2026-09-08]].

**Update (late evening 2026-09-08).** v10 = fresh full run on the final build: 1042 s end to end (paired 653
incl. deconvolution 157 s on the 20k-subsample ladder; states 389 s adopting `flex_pca_embedding_deconv.bin`),
vs 61 min for v9; science as v9b (4 macro-clusters, d' 3.67 / purity 0.70). v10b = v10 embedding with K
ceiling 16 (K=16 chosen -- the held-out curve keeps rising, so K is set by the ceiling and the min_neff
fold decides the states), npreimages=8, maps off: the island's two clumps are separate states (3.8k and
7.3k), purity 0.71, detached 0.58; the deconvolution at K=16 took 454 s. Delivered under
`~/Downloads/flex_paired_20260908/ribo10028_v10` (maps) and `ribo10028_v10b_k16` (coordinates, figures).
Next natural run: v10b settings with maps on (~25 min). Box branch flex-pca-resident tip bf02ca00a.

**Sectioning a continuum is geometric, not probabilistic (tested and reverted 2026-09-08 late).** Fitting a
bm-component mixture INSIDE a continuum macro-cluster collapses to one component whichever noise model is
used: with the posterior covariances as noise (v10c) because the shrunk coordinates plus their posterior
covariance double-count the spread (sections of width 3e-3); with the raw coordinates through the
calibrated measurement model (v10d) because the top-level mixture already found that cluster to be ONE
Gaussian, and a Gaussian has no internal modes (98.4% in one section). Both commits reverted
(2cf673b5a); the arc-length Epanechnikov sections along the cluster's dominant axis (v10b) stay. If
finer sectioning is wanted, it has to be equal-mass regions on the denoised coordinates (k-means /
Voronoi with posterior-precision distances), never a likelihood fit.

**Continuum sectioning, final form (2026-09-09 00:17, branch tip 7f817efd2).** `continuum_region_weights`
in simple_flex_pca_gmm: k-means targets (targets module, deterministic farthest-point seeding) on the
denoised coordinates in the STANDARDISED metric over all nk dimensions (so a breathing axis does not
dictate the cut -- the user's point), each particle assigned to its nearest target with unit weight
(Voronoi cells partition the cluster: v12 regions 16.5k/14.7k/19.6k, neff = count). A compact kernel of
half the target spacing was tried first and held 10% of the cluster in 17 dimensions (v10f). v11 (one-axis
arc-length sections) gave three identical continuum maps at 3.86 A; v12 gives three disjoint-population
maps at 4.05-4.19 A. Delivered: `~/Downloads/flex_paired_20260908/ribo10028_v11_k16_maps` and
`ribo10028_v12_voronoi_maps` (8 states + e/o halves each). Reminder that bit me twice today: never
`pkill -f`/`pgrep -f` a pattern that the ssh command line itself contains (it kills the ssh shell); kill by
PID / process group from a script.

**The continuum is one structure at the current basis resolution (2026-09-09 ~01:00).** At a common 8 A over
the envelope (`~/state_diff8.py`), the v12 continuum regions differ from each other by 1.8-2.0 half-map noise
units while discrete populations differ from them by 3.6-4.1. A ZOOM (flex on the 50.8k continuum particles
alone via `pindfile=` on the ORIGINAL project, `zoom_continuum_v1/v2`) did not help: its own deconvolution put
99.8% in one Gaussian, its six k-means regions differ by 1.6-2.5 noise units, and its UMAP is a round cloud.
So a within-population basis pass is NOT justified. The user's hypothesis, which fits: ratcheting is a fine
motion (~1 px at box_crop 64 / lp 16), below what the basis grid resolves -- v13 tests box_crop=96 lp=12
(`queue28`, dir paired_rankcut_br360_v13_crop96); SIMPLE_COV_SNR_SHELLNORM no longer exists (dead-code cut).
Also: a lone mixture-defined continuum never falls back to the tied GMM (decb0b6b5); the multi-state input
check is now worded correctly (a rerun must start from the ORIGINAL project; select with pindfile=).
Cleanup 2026-09-09: superseded runs/gates removed on the box, plane caches 92 -> 27 GB (one per kept run),
local Downloads 26 -> 12 GB (maps kept only for v11 and v12; splice10180_v5 kept).

**Continuum experiments, closed 2026-09-09 morning.** Fine band on everything (v13, box_crop 96 lp 12):
calibration a=7.9, populations dissolve back into the body, island d' 3.66->2.64 -- dead end. Fine band on the
continuum alone (zoom v3): a=10.2, K=1, deconvolved coordinates collapse (regions at distance 1e-4) --
dead end. Focused periphery (zoom v4, SIMPLE_COV_FOCUS_R1=80 FOCUS_EDGE=15, coarse band): the latent finds 4
populations in the continuum (7.7k/5.6k/25k/12.5k, axes 4-35% signal) but their maps differ by 1.8-2.3
half-map noise units at 8 A even inside the shell beyond 80 A (`~/state_diff8_shell.py <dir> 6 1.31 80`),
same as the plain zoom (1.6-2.4). Noise is NOT solvent: eigenvolumes 84-94% inside the envelope, half-map
noise flat with radius (2.4), signal falls 12.9 core -> 4.4 at 120-140 A, so the periphery (40S) has the
worst SNR. Conclusion: within the continuum nothing separates beyond ~2 noise units per 5-17k-particle
state; the ratchet, if inside it, needs a DISPLACEMENT representation (per-domain rigid-body pose search of
the 40S per particle, the "flex fades, refine3D displaces" lesson), not a finer or focused density basis.
The resolved heterogeneity of 10028 is the discrete populations (v11/v12).

**Breathing, and the final word on the continuum (2026-09-09 10:20).** The focused-periphery zoom's four
"populations" were the dilation mode: with the consensus dilation (x-c).grad rho deflated
(`dilation_template` in util, appended to the deflation set in mstep + pairmerge, default on,
SIMPLE_COV_DEFLATE_DILATION=0 opts out; commit 0fcd70404), the states inside the continuum deviate from
their consensus by 1.1-1.4 noise units and from each other by 1.6-2.1 at 8 A (zoom_continuum_v5_focus80_dil,
delivered under ~/Downloads/flex_paired_20260908/). Every route -- coarse, fine band, zoom, focused,
focused+dilation-free -- says the 51k continuum of 10028 is one structure at this resolution. Keep the
dilation deflation on for every run (it is a nuisance in general). Local checkout switched to branch
flex-pca-resident (master WIP in stash@{0}); the scratch worktree is gone -- edit in /Users/atifao/SIMPLE.
Pending user decision: the flex cleanup (comments to intent lines, merge small submodules, planes+cache,
rename state merge, delete dead placement paths; ~22 files / ~17k lines, gated bit-identical).
