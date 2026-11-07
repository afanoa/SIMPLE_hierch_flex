---
name: empiar-10509-ribosome-inputs
description: "EMPIAR-10509 70S ribosome (2 A, Watson 2020): raw K3 super-res TIFFs only, both gains identical, one broken movie, EMPIAR dose/count typos, gain flip x; work dir /mnt/beegfs/elmlund/afan/ribosome_10509"
metadata:
  node_type: memory
  type: project
  originSessionId: a9a4000f-ee34-4500-a4e0-aed13b26e342
  modified: 2026-10-20T15:15:05.798Z
---

Inventoried 2026-10-19. Download: batch2 `/mnt/beegfs/elmlund/EMPIAR_downloads/10509/data/` (lowercase downloads,
fetched by meanapanedar2 2026-10-16). Work dir `/mnt/beegfs/elmlund/afan/ribosome_10509` (empty at that date).
- Movies only: 70S-Complex-1 3,880 TIFF (XML says 3,380 = typo; EMPIAR FTP also has 3,880), 70S-Complex-2 3,394.
  K3 super-res 11520x8184, 40 frames, 8-bit LZW. Gains: `Complex-[1/2]-SuperRef_*.mrc` (float32 super-res).
- Every file byte-matches the EMPIAR FTP listing except `70S-Complex-2/19Jun25_5f_162-27_0008.tif` (4 KB vs 323 MB):
  exclude or re-fetch. 65 movies ~57 MB and 19 ~510 MB are the same size in the archive (as deposited).
- TRAP: 9 movies are DEPOSITED TRUNCATED (39 of 40 frames; frame-40 IFD/strip table past EOF, sizes 4096-multiples
  identical to the archive), so a byte-size check passes them. One killed preprocess part 5 (TIFFFetchStripThing
  StripOffsets). The movies are BigTIFF; scripts/tiff_integrity.py walks every IFD. 7,264 good movies.
- TRAP 2: 12 more movies are structurally intact but have corrupt LZW strips (rTiffSlices 'Failed to read from
  strip'); only a full decode finds them (scripts/decode_check.sbatch: info_image per movie, 48-way xargs, ~2 min/movie/core).
  Usable: 7,252 movies. Square Complex-2 124 and 162 hold most of the bad ones.
- Optics (paper methods, PMC7550191): 300 kV, Cs 2.7, super-res 0.355 / physical 0.71 A, calibrated 0.7118
  (super-res 0.3559), total dose 39.89 e/A2 (EMPIAR's "0.09972 e/A2/frame" is mislabelled counts/px/frame), defocus -0.6..-1.5 um,
  authors binned 2x to 0.71; 874,943 initial / 307,495 final particles. EMD-22586: map 0.7118 A/px (496^3), 1.98 A,
  dose 40.0. MEASURED: raw movies average 0.100 counts/px/frame = super-res 0.3559 (0.79 of 40 e/A2; at 0.7118 it would
  be 0.20), and EMPIAR's 0.09972 is counts/px/frame. User pushed to import at 0.7118 (2026-10-19); evidence says 0.3559.
- Gain "flipped about vertical axis" = SIMPLE flipgain=x, verified on 20 movies (frame-average corr: flip_x -0.89,
  flip_xy -0.76, flip_y -0.34, unchanged -0.28). image%flip off-by-two (ii=ldim-i-1) was fixed upstream in e165e2aea
  (cyrilr 2026-10-02); branches older than that (flex2D, batch2 ~/SIMPLE) still carry it.
- Processing: batch2 tree ~/SIMPLE_proc = origin/master aea5b1846 (dedicated, do not rebuild under jobs).
  USER'S TARGET: everything downstream at 1.42 A = 0.7118 physical binned 2x ("import at .7118 ... binned 2x to 1.42";
  I wrongly ran 0.7118 micrographs first and argued instead of hearing the goal). Implemented correctly as import
  smpd=0.3559 + preprocess smpd_downscale=1.4236 total_dose=40.0 -> 2880x2046 mics. Current run: job 54004424
  (scripts/import_preprocess.sbatch, 7,252 readable movies). Earlier cancelled/failed runs are in
  /mnt/beegfs/elmlund/afan/ribosome_10509_discarded/ (user may delete). A single bad movie aborts the whole
  distributed preprocess. 3_selection = 6,883 mics (ctf_estimatecc >= 0.2). Pickrefs branch 4-8 (100 mics, segdiam,
  box 256, solve2D 50 cls; user kept 33 -> 8_selection/cls2D_selected.mrcs). 9_pick on 3_selection: 800,098 picks.
  Ribosome measured ~235-260 A; user chose box 264 -> 10_extract: 710,631 particles (2026-10-20). Paper: 874,943
  initial, 307,495 final.
  11_solve2D 300 cls (all 710k) -> user kept 231 (575,870) -> 12_selection -> 13_cls_expansion ncls=23100 ("100 per
  class"; ncls is the TOTAL): 22,919 subclasses, ~25 ptcls each, repro median 0.59. Honeycomb-sphere contaminant is common.
- The two SuperRef gains are the SAME gain: pixel data byte-identical (cmp past the 1024-byte header), only the
  EMAN label date differs (6/26 vs 6/28/2019). User 2026-10-19: use ONE gain for everything —
  `data/Complex-1-SuperRef_19Jun24_mmal_88-3_0001.mrc` — for both sessions, in one preprocess run.
Related: [[persist-project-lineage-no-new-projects]], [[head-node-runs-nothing-but-sbatch]].
