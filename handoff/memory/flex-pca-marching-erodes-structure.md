---
name: flex-pca-marching-erodes-structure
description: "CODE REMOVED 2026-09-10. 2026-09-09 afternoon: box-to-box frequency marching (tethered continuation 64->80->96->112) makes the box-64 clusters LESS separable at every finer band, with or without the rank cut (fixed-label kNN purity 0.87 -> 0.71 -> 0.69 -> 0.65 with the cut, 0.63/0.62 without); the cut helps but does not stop it; the per-particle S/N threshold of 1 keeps only 2 axes (use 0.1 -> 8-9); box 128 collapses the working rank to 1 and crashes in the prune; the distributed reduce hung on a mid-stage rank change (fixed)"
metadata: 
  node_type: memory
  type: project
  originSessionId: 7e382bce-970c-4cdb-bc32-4a1a1c1801a9
  modified: 2026-09-09T18:06:30.885Z
---

**Setup (branch flex-delivery, box branch flex-pca-resident-in, tip 292911a42):** `SIMPLE_COV_TETHER=<prev
stage's flex_pca_embedding_deconv[_cut].bin>` re-solves the basis on the new grid from the previous
posterior (iteration 1 bypass), one tethered EM step, merge + joint step, then a PLAIN embedding
(`SIMPLE_COV_TETHER_EMBED=1` opt-in only: tethering the halves collapses the calibration, a=0.27, K=1).
Stage-boundary cut `signal_subspace` (noise-whitened population variance > `SIMPLE_COV_CUT_SNR` x mean
noise, default 1.0) writes `flex_pca_embedding_deconv_cut.bin`; `SIMPLE_COV_STAGE_CUT=0` disables.
Scripts: ~/queue36.sh (no cut), ~/queue37.sh (cut), ~/queue38.sh (fresh box-48 start, never run),
~/queue39.sh (extension-only test). Readouts: ~/montage_ladder.py (per-stage UMAP montage with a, K, kept
rank, island score), ~/montage_watch.sh, ~/fixedlabel_score.py (kNN purity + between/within of FIXED
box-64 labels in each stage's deconvolved latent -- the frame-free comparison; UMAP layouts are not).

**Results (10028, 105k):**
| stage | fixed-label kNN | B/W | island d'/purity |
| box 64 origin (self) | 0.87 | 1.19 | 3.66 / 0.71 |
| 80 no cut / 96 no cut | 0.63 / 0.62 | 0.99 / 0.37 | 3.18/0.67, 3.10/0.64 |
| 80 cut r8 / 96 cut r8 / 112 cut r7 | 0.71 / 0.69 / 0.65 | 0.60 / 0.58 / 0.42 | 3.02/0.71, 3.07/0.70, 2.96/0.60 |
Small clusters dissolve first; per-direction S/N RISES up the ladder (1.8 -> 3.7 -> 4.6 on axis 1) while
separability falls, so the finer band adds population variance that is not the coarse structure. The
S/N=1 threshold kept 2 of 20 at box 64 (spectrum 1.77 1.17 0.90 0.41 0.26 0.20 0.15 0.13 0.10 0.04 ...);
0.1 keeps 8-9 = the CV+ICA rank. Sensible rule: cut ONCE at the origin at the knee, hold the rank.
Box 128 (lp 9): d_tilde energy floor collapsed the fit to ncomp=1, merge |cos|=0, `double free` in the
prune (SIMPLE_COV_MERGE_PRUNE) -- degenerate-rank crash, unfixed.

**Open question the extension test answers (marchext_b80_lp14, SIMPLE_COV_TETHER_EXTEND=1, single fit,
n_probe_iters=1 = basis re-solve only, plain embedding):** if it scores ~0.87 the EM refit at the fine band
is the culprit (march by basis extension, never refit); if ~0.71 the plain fine-band embedding drowns the
coordinates and the embedding must be tethered with a split noise model (data part calibrated, inherited
part carried).

**Infrastructure fixed today:** master hang = stale mixture accumulators after a rank change misaligned
the part stream (1eee41d14; reader now sizes from the header); serial part reads -> concurrent
(128654bb1); threaded FFT plans on the big reconstructors + M-STEP SPLIT timing (ea1a75d57). Worker
E-step is 92% efficient; the coupled K x K insert is 70% of a worker pass; the exchange is ~2.5 GB/part
at box 112 on the padded grid (compress-before-exchange would cut it 8x; not done).
Related: [[flex-pca-delivery-pooling-and-tether]], [[flex-pca-polish-destroys-island]], [[flex-pca-cv-ica-recipe-2026-09-01]].

**Extension-only answer (14:19):** basis re-solved from the box-64 cut tether at box 80, NO E-step refit,
plain embedding: calibration a=5.46, fixed-label kNN 0.60 (worse than the refit's 0.71), small clusters
gone (7:0.22 12:0.12 9:0.04), island d' 2.78. So the culprit is the plain fine-band embedding itself:
the coordinates at the finer band are dominated by per-particle noise, and the newly opened shells of a
basis extended from a coarse latent are noise-fitted (the refit cleans them, a 1.7). Latent marching to
finer bands does not refine the coarse-band clusters on 10028 at this SNR; the fine band belongs to the
maps (delivery already reconstructs states at full box). Only a TETHERED embedding (split noise model:
data part calibrated on untethered halves, inherited part carried) could keep the coarse structure while
adding fine-band information, and its expected gain is small. Single-fit path: SIMPLE_COV_PROBE_ONLY=1
stops BEFORE the embedding (the paired path does not) -- drop it for single-fit runs.

**Tethered-embedding arms (15:20), final table (fixed box-64 labels, kNN purity; origin 0.87):**
plain 0.71/0.69 (80/96); posterior-tether 0.76/0.76 but box 96 == box 80 to 1-2% (inherited precision 5e7 vs
data 2-5e5: the mixture posterior compounds, chain freezes -- a shrinkage artefact); raw-measurement tether
(`flex_pca_embedding_raw_cut.bin`, data precision / a) 0.73/0.73 with B/W 0.30 and island purity 0.65
(worse than plain 0.71). 30 A origin (box 40): only 2 directions above 0.1 x noise, K=5 -- the coarse start
is weaker, not cleaner. VERDICT: latent marching is closed on 10028; keep the origin rank cut; fine band ->
maps only. Env switches that exist now: SIMPLE_COV_TETHER, _TETHER_EMBED, _TETHER_EXTEND, _STAGE_CUT,
_CUT_SNR; pose refinement (P1 in-plane/shift, marginal likelihood) lives at git 44ed91143^ em_pose.f90.
Next lever chosen with the user: agnostic zoom on the head-less state (two recipes) + mask-diameter arm.

**Floor measured (21:01):** an exact repeat of the sphere origin scores 0.72 against the original labels, so 0.71/0.69/0.65 were never below chance-of-rerun; the solvent-mask chain (own labels 0.83 -> 0.57 -> 0.52 -> 0.49) and the box-40/30 A origin (0.43) close marching for good. Successor: [[flex-pca-basis-composition]].

**REMOVED 2026-09-10 (flex-delivery, after the compose results):** SIMPLE_COV_TETHER / _TETHER_EMBED / _TETHER_EXTEND, the tether fields of probe_fit_t, tether_read/tether_load, the E-step bypass and Pt/Ptm solver terms, SIMPLE_COV_STAGE_CUT and the flex_pca_embedding_deconv_cut.bin / _raw_cut.bin files, and deconvolve_latent's wcut/kcut/tether_prec optionals. signal_subspace + SIMPLE_COV_CUT_SNR survive only for the composed-basis cut (compose_cut_reembed). The env switches named above in this memory no longer exist.
