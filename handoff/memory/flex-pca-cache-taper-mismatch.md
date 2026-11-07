---
name: flex-pca-cache-taper-mismatch
description: "Why cache=yes changed flex_pca numerics (rank 7 -> 13-15, purity 0.48 -> 0.40 at 20k): the shared particle cache re-tapers and pads a Fourier-cropped 64-px image, a different prep from the full-box path (taper width AND wrapped edge ringing); RESOLVED by a flex-owned PLANE CACHE storing the full-box prep's padded transform on the box_crop grid: gate s1c vs s1 = 1.5e-5 (rounding) vs 0.9 before"
metadata:
  node_type: memory
  type: project
  originSessionId: 7e382bce-970c-4cdb-bc32-4a1a1c1801a9
  modified: 2026-09-08T18:44:50.534Z
---

The particle cache entry is the noise-normalised, Fourier-cropped particle at box_crop (written by
`norm_noise_fft_clip_shift` + ifft in simple_ptcl_cache). flex's cached prep then calls
`norm_noise_taper_edge_pad_fft(lmsk_crop, ..., renorm=.false.)`, whose taper/edge-mean strips are
`winsz = nint(COSMSKHALFWIDTH)` = 6 PIXELS regardless of box: 6 px at box 64 (7.4 A/px) is a
44 A window, the raw path's 6 px at box 360 is 8 A. Per-particle residual energies agreed to 0.1%
while eigenvalues differed 43% on the gate; at 20k the rank cut went 7 -> 15 and purity 0.48 -> 0.38.

**Fix (box branch flex-pca-resident, 1e46ddb7e + d0feeaec8):** optional `taper_width` on
`norm_noise_taper_edge_pad_fft` (interface in simple_image.f90, body in simple_image_fft.f90; the
one necessary shared-code edit), passed as `nint(COSMSKHALFWIDTH*crop_factor)` in the cached branch
of `prep_imgs4projected_model`; and the mean + embedding passes now also read the cache
(`init_rec(cropped=)`, `prepimgbatch(box=box_crop, smpd=smpd_crop)`), because the resident store
asserts one read path per run. Verification = the 20k cache arm at 10 x 8 against the raw 10 x 8 arm
(rank 7, d' 2.21, purity 0.48); pending at the time of writing (queue9b). Remaining known
differences: real-space pad fill after a Fourier crop vs the full-box pad (aliasing-level).
Related: [[flex-pca-distribution-changes-results]], [[flex-pca-speed-levers-2026-09-08]].

**Resolution (2026-09-08 evening):** the taper width was only part of it (rank 15 -> 13). The exact fix
is `simple_flex_pca_plane_cache.f90` (box branch flex-pca-resident, 276e61f39 + 49c29e9f4): per particle
the (fdim(128) x 128) complex block of the 720-padded transform = the cmat a 128-grid image carries
(same physical frequencies, both scaled 1/n^2), injected into the heap image with set_cmat/set_ft and
run through gen_fplane4rec unchanged. Gate: s1c vs s1 z max|d|/max = 1.5e-05 (rounding) against 0.88-0.94
with the shared cache. Direct-access file under cache_dir, record p+1 = row p, header LAST (a build that
dies mid-way leaves a zero header) plus a size check on adoption -- a header-only leftover from a crashed
build was adopted once and killed two arms. The shared `simple_ptcl_cache` is no longer used by flex; the
image-module taper_width edit was reverted. Residency (`simple_flex_pca_planes`) is the in-memory face of
the same object and only engages with cache=yes.
