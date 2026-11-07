---
name: cnga1-has-reproducible-continuum
description: "2026-09-13: cnga1 (245,424 ptcls) DOES have a cross-half reproducible flex_pca basis (comps 1-4 at |cos| 0.87-0.997, unlike PfCRT's 0.09) but the latent is a pure CONTINUUM -- unimodal marginals, smooth z1/z2 gradient, no gaps -- so the discrete macro-cluster readout will always look like nothing"
metadata:
  type: project
---

**Run:** job 53201193 on batch2/csbdevel, `1_flex_pca` under
`/mnt/beegfs/elmlund/afan/cnga1/flex_pcg_env`. 18 min wall, 10 parts x 16 thr, embedding only
(`SIMPLE_COV_STATE_REC=0`). `rec_backend=pcg maxits_pcg=8 rtol=1e-3` with
`pcg_mskfile=10_refine3D_auto/automask3D_state01.mrc` (the refine3D automask, 8.0% of box 288),
`neigs=10 n_probe_iters=4 box_crop=64 lp=16 box_rec=288`, mskdiam 180 from the project.

**The basis is real.** Cross-half matched |cos| (`flex_pca_paired.txt`, fit A / fit B):
comp1 0.997/0.995, comp2 0.942/0.913, comp3 0.885/0.872, comp4 0.894/0.957, then
5 0.777/0.841, 6 0.907/0.690, 7 0.717/0.664, 8 0.869/0.644, 9 0.878/0.732, 10 0.863/0.615.
So ~4 solid axes. Contrast [[pfcrt-injection-floor]] where the real data gave cos 0.09.
Merged eigenvalues are FLAT: 25.5 23.1 20.2 21.0 19.8 19.1 19.5 18.1 18.3 15.8 (only 1.6x range,
and not monotonic) -> no dominant mode, variance spread over all 10.

**But the latent is a continuum, not states.** z1..z4 marginals are clean unimodal Gaussians
(median-mean within 0.017); the UMAP blob shows a SMOOTH z1 and z2 gradient across it with no
density minimum. GMM AUTO reported kfit=4 macro-clusters (51303 / 99149 / 86682 / 8290) and then
labelled 3 of the 4 as continuum needing k-means regions. Mean max-responsibility 0.589 with
**29.4% of particles below 0.5**: the three big clusters are an arbitrary partition of one blob.
Deconv CV gained only 0.083 nats/particle from K=1 (-18.271) to K=4 (-18.188).
Signal fraction (deconvolved/observed) came out 0.92-1.13, i.e. ~1.0, because prior variances
(19.7 ... 11.8) dwarf sig2_eff=0.189 -- the deconvolution had essentially no noise to remove.

**Two real satellites, both small and suspect.** Separated islands at UMAP (-17,30) n=2076 and
(-7,10.5) n=1735 (together 1.6% of data), both inside macro-cluster 4. Same `corr` as the bulk
(0.374 vs 0.375) so not junk by correlation, but defocus is shifted (1.96 and 1.88 vs 1.78 um)
and e1 is narrower (18.8 and 23.5 vs 25.3 sd), so view geometry / CTF is not excluded.
**The per-micrograph test is impossible on this project: every particle has `stkind=1`** (imported
from one concatenated `final_particles.mrcs`), so micrograph concentration cannot be checked.

**Next readout (NOT run, user deferred reconstruction):** continuum data needs equal-mass bins
along the reproducible axis, i.e. `state_axis=1` sextiles, exactly the fix from
[[afm-rna-flex-pca-needs-continuum-readout]] and [[flex-pca-gmm-undoes-equal-mass-placement]] --
not the discrete macro-cluster maps. Artifacts: `umap_deconv_by_component.png`,
`umap_deconv_core.png`, `cnga1_embedding_summary.png` in that run dir and on the Mac at
`~/Downloads/cnga1_flex_pcg_20260913/`.
Related: [[simple-script-yes-compenv-trap]], [[flex-pca-pcg-branch-state]].

---

**CORRECTION 2026-09-14 — the `state_axis=1` recommendation above is REVERSED for cnga1.** Both
placements were run at box_crop=128 / lp=8 / rank 10 and scored by the cross-halfset same-state vs
different-state gap on mean-subtracted maps (the test in
[[flex-pca-state-reproducibility-gap-test]]):

| placement | same-state | different-state | **gap** |
| `state_axis=0` (Voronoi in the full latent), all 8 | 0.588 | -0.035 | **0.622** |
| `state_axis=0`, 6 well-populated only | 0.482 | -0.079 | **0.561** |
| `state_axis=1` (equal-mass z1 sections), all 8 | 0.272 | -0.038 | 0.310 |

ax1's per-state values are 0.43 0.32 0.18 0.05 0.03 0.18 0.41 0.59 -- **three of its eight sections
are noise**, which is why its maps looked nearly identical. ax0 wins even after controlling for
occupancy. So "continuum data needs equal-mass bins" did not hold here; verify per dataset.

**But ax0 needs an occupancy guard.** Its 8 states had neff 96815 / 113916 / **7781** / 121261 /
115665 / **15101** / 119700 / 116406: states 3 and 6 are 8-15x under-populated, are the two visibly
blurrier maps, and pairwise map cosine tracks occupancy at **r = -0.660**, so every one of the top-5
"most different" pairs was a low-neff vs high-neff comparison. Exclude them before reading anything.
Among the six well-populated states the dominant difference axis carries **57 %** of the
between-state variance (second 27.6 %) and orders them **1, 8, 2, 7, 4, 5** -- a one-dimensional
ordering, consistent with the continuum verdict. Artifacts: `states_ax0/montage_ax0.png`,
`montage_ax0_axis.png`; scripts `cl_gap.py`, `cl_neff.py`, `cl_axis.py`, `cl_rad.py` in
`flex_pcg_env`. What the ordering actually is: [[cnga1-continuum-is-the-micelle-halo]].

**2026-09-17:** neigs=10 on cnga1 (box 128, pcg+mask, 4 EM iters) gives 9/10 components reproducible
at both-halves 0.7 (meanmin 0.808) vs 3/32 at neigs=32 (0.477), 6x faster per iteration; neigs=10 is
the recipe of record. Over-ranking dilutes reproducibility here as it did on TRPM4 (K=20 -> 2/20).
