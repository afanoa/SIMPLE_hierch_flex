---
name: flex-pca-basis-composition
description: "2026-09-10: SIMPLE_COV_COMPOSE (commits cbd5e58a4, fc0848596, a4660e89d on flex-delivery / box flex-pca-resident-in) composes an lp-agnostic basis from finished runs; the compose path is exact (control reproduces the origin to 3 decimals), half-fit shell weights HURT (opt-in now), the 24+20+16 A union gives 3 distinct islands (d' 3.21, K=11, 30 columns after the whitened cut) but keeps the coarse bands' view (0.63 vs box-64 labels); marching closed (floor 0.72)"
metadata:
  type: project
  originSessionId: 7e382bce-970c-4cdb-bc32-4a1a1c1801a9
  modified: 2026-09-10T01:57:13.586Z
---

**Why composition:** PCA at a band returns the leading eigenvectors of the covariance projected
onto that band, and adding shells re-orders them; marching one basis through the bands re-orders
it at every stage. RECOVAR never marches: one band (downsampled box, default 256), ~300 covariance
columns at |k| <= column_radius=5 each estimated to full band, per-column per-shell RELION-style
FSC regularisation (compute_covariance_regularization_relion_style), SVD of the columns, PCA of the
projected covariance in that span. The per-column FSC shrinkage is their "reduce noise in finer
shells"; composition ports it at component level.

**What ships (src/main/flex/simple_flex_pca_em_compose.f90, wired in the model's block_eigenbasis,
refuses SIMPLE_COV_PAIRED):** `SIMPLE_COV_COMPOSE=<dir>,<dir>,...` finished run dirs (any order,
sorted by crop box, coarse first; namespace polished > merged > plain with the matching
flex_pca_probe*.txt for rank/sig2/prior variances). Per component: optional per-shell weight
2c/(1+c) from the FSC of its matched fit-A (`flex_pca_pc*`) and fit-B (`flex_pca_fitB_pc*`) half
fits (c<0.143 -> 0, shells < 6 exempt, match |cos| >= 0.5 required; `SIMPLE_COV_COMPOSE_FSCW=0`
off); Fourier-pad to box_crop with coefficients copied verbatim (fft is 1/N-normalised, so a
coarse column's real-space values are the band-limited interpolant; prior variance x fine-grid
norm^2 after renormalising); coarse-first Gram-Schmidt in real space, columns with residual
fraction < 0.05 dropped, prior variance x residual fraction; sig2_eff = finest run's. Writes
flex_pca_pc*.mrc + flex_pca_probe.txt (distributed embed workers load them) + flex_pca_compose.txt.
No covariance moment, no EM: mean -> compose -> the plain single-fit embedding -> deconv/cut/UMAP.
Run with SIMPLE_COV_RELPRIOR=0 to mirror the paired path's plain-prior embedding.

**Evening controls (10028, fixed box-64 labels from marchcut_b64_lp16 unless noted):**
- repeat of the sphere origin (repeat_b64_lp16): 0.72 -> the reproducibility FLOOR of the score.
- sigma-shape recal origin + box 80 (sig2_*): box 80 0.63 vs plain box-80 cut 0.49 (same labels).
- solvent-mask marching (solv_b64 -> 80 -> 96 -> 112, its own labels): 0.83 -> 0.57 -> 0.52 -> 0.49.
- box 40 / 30 A origin: 2 directions above 0.1 x noise, island d' 1.32, 0.43 vs box-64 labels.
Marching is closed: with or without the proper solvent mask it erodes below the floor.

**Queued (queue54 -> queue55, dirs under flex_h0901):** coarse origins coarse_b40_lp24 and
coarse_b48_lp20 (6-iter sphere fits, scored both ways), then compose arms: compose_ctrl_b64
(repeat's own basis re-embedded, FSCW=0 = machinery control, expect ~self score vs repeat labels),
compose_ctrlw_b64 (weights only), compose_40_48_64 (box 64), compose_64_80 and
compose_40_48_64_80 (box 80, lp 14, marchcut_b80_lp14 as the fine source). Score against BOTH
marchcut (floor 0.72) and repeat labels (self chain).
Related: [[flex-pca-marching-erodes-structure]], [[flex-pca-delivery-pooling-and-tether]].

**Night of 2026-09-09/10 results and shipped fixes:**
- compose control (repeat origin re-embedded, no weights): 0.768 / 0.804 = identical to the origin -> the path is exact.
- half-fit FSC shell weights: control drops to 0.72 / 0.71, the two leading components' half fits
  disagree beyond 80 A (mean weights 0.12, 0.35) -> weights are OPT-IN (SIMPLE_COV_COMPOSE_FSCW=1).
- 58-column union without a cut: the deconvolution K ladder takes >1 h (cost ~ d^2) -> shipped
  compose_cut_reembed (fc0848596): embed the union once, calibrate, moment-based population
  covariance, signal_subspace at SIMPLE_COV_CUT_SNR, rotate the basis onto U W+ and re-embed
  (SIMPLE_COV_COMPOSE_CUT=0 skips). noise_and_projection + signal_subspace are now public in deconv.
- union 24+20+16 A at box 64 (compose_40_48_64_nw): 57 columns -> cut 30 -> re-embedded; K=11,
  island d' 3.21 / purity 0.61, THREE distinct islands + a blob in the UMAP; vs box-64 labels 0.63
  (= the coarse origins), vs repeat labels 0.58, B/W 0.16: the coarse-first union inherits the
  coarse bands' view (fine columns enter as residuals with r2-scaled priors).
- coarse origins: box 40/24 A d' 3.25 (0.63 vs box-64, box-64 keeps only 0.48 of ITS clusters);
  box 48/20 A d' 3.30 (0.64 / 0.46). The bands separate partly different populations.
- SIMPLE_COV_COMPOSE_FINE_FIRST=1 (a4660e89d) reverses the band order.
**Queued overnight (logs ~/queue56..59.log):** box-80 unions (64+80, 40+48+64+80); union
reproducibility (coarse refits *_r2 + marchcut origin -> compose_40_48_64_r2 scored vs union-1
labels) + union-1 STATE MAPS (rank30 control = rank30_b64_lp16, a single box-64 fit with neigs=30);
fine-first and CUT_SNR=0.3 variants. Next code step agreed with the user: DISJOINT band-pass fits
(lower shell from an env hp) so the union is orthogonal by construction and the per-particle Gram
is block-diagonal (Fourier slice theorem); touches projected_model_kfromto + covariance_kfromto,
the 4 bp(0.,lp) calls, the hybrid polar E-step ring range and the GPU prep -- do it attended.
RECOVAR for reference: one band, ~300 low-frequency-indexed full-band covariance columns with
per-column per-shell half-set FSC priors, PCA in their span; no marching.

**Reproducibility verdict (02:22):** union 2 (independent coarse refits + the OTHER box-64 origin)
scores 0.60 against union-1 labels (union-1 self 0.84, chance 0.20) and union 1 scores 0.53 against
union-2 labels (self 0.79); the three largest union clusters reproduce at 0.77/0.57/0.66. Any
single-band origin scores only 0.48 against the union labels, so the union carries structure no
single band has, but its own repeat agreement (0.53-0.60) is BELOW the single-band floor
(0.72-0.77): the islands are partly reproducible, not established. Nested unions with box 80:
64+80 0.44 (cut kept 7 of 28), 40+48+64+80 0.58 -- each nested fine band dilutes. The first
state-map resume failed (delivered project carries state labels; copy the original project again
before infile= resumes) -> requeued as queue60 behind queue59 (fine-first, CUT_SNR=0.3) and
queue58 (rank-30 single-band control rank30_b64_lp16).

**Morning table (all three-band unions at box 64; labels: union-1 / box-64 origin / repeat origin):**
| arm | kept dims | island d'/purity | vs union-1 (self 0.84) | vs box-64 (floor 0.72) | vs repeat (self 0.80) |
| union 1 coarse-first cut 0.1 | 30 | 3.21/0.61 | self | 0.63 | 0.58 |
| union 2 independent refits | 28 | 3.34/0.65 | 0.60 | - | - |
| fine-first | 16 | 3.11/0.64 | 0.63 | 0.63 | 0.58 |
| CUT_SNR=0.3 | 14 | 3.04/0.64 | 0.63 | 0.62 | 0.58 |
| rank-30 single band (control) | 30 | 3.63/0.66 | 0.48 | 0.68 | - |
| any 20-comp single band | 20 | | 0.48 | | |
VERDICT: the union's clusters are a BAND effect (rank-30 single band scores 0.48 vs union labels,
same as rank 20; independent union 0.60), partly reproducible (below the 0.72 single-band floor);
ordering and cut threshold are neutral on purity, CUT_SNR=0.3 is 3x cheaper (14 dims, 8 min) ->
use 0.3 for unions. Union 1 state maps (8 states, 13 min, compose_40_48_64_nw/flex_pca_state_00?.mrc
+ _diff_): occupancies 5.2k/8.9k/40.9k/6.0k/5.2k/22.6k/7.1k/9.3k, own-kernel detail to 25/27/15/25/
25/16/23/13 A, common blended res 3.28 A. NOT yet looked at. Next: disjoint band-pass fits (attended).

**2026-09-10 morning, shipped:** (1) 3D rank gate at compose (2bcdd0277): per band keep only components
whose two half-fit eigenvolumes agree to cross-FSC >= 0.5 over the band (read from flex_pca_paired.txt
"# merge pair" via "# merge matched"; SIMPLE_COV_COMPOSE_GATE=0 off). Measured: the reproducible set is
NOT a prefix of the variance order (box-64 origin: comps 1-8,12,17,20 pass, 9-11 fail; rank-30 fit: 17
pass; box-40/24: 15 pass); the latent 0.1 cut kept 8/16/10. The merge's principal-angle subspace gate
is saturated (all >= 0.7) and useless as a cut. (2) SIMPLE_COV_HP band-pass fits (e6e31cd00): util
cov_band_hp/cov_band_klo/cov_band_limit; planes zeroed below the edge after prep (highpass_fplane in
latent_ops, post-cache so the plane cache stays exact), mean_rec band-limited in init_mean_reconstructor,
the four basis band limits (mstep x2, pairmerge, basis) via cov_band_limit (one-shell edge), init
impulses inside the band, covariance_kfromto(1) = klo; device prep path refuses it. Run band-pass fits
with EM_DEFLATE=0 DEFLATE_BG=0 (no nuisance power in a high-passed band). Facts: there was NO high-pass
anywhere before (band = shell 1 at ~480 A .. lp; lp at the crop Nyquist is inert, the crop box IS the
band); ~93% of basis power > 60 A so every nested band is dominated by the same low shells.
Queued: queue61 = three-band union with gate (+/- latent cut); queue62 = band fits 24-20 A at box 48,
20-16 A at box 64 (band_b48_24_20, band_b64_20_16) + disjoint union compose_band_40_48_64.
Band schedule idea (not built): per-shell cross-half FSC of the coarse fit's kept components ->
sum FSC/(1-FSC) = heterogeneity signal per shell -> band edges at equal cumulative signal, stop where
the coarser fit's components stop reproducing (the abinitio3D class-FRC ladder analogue).
