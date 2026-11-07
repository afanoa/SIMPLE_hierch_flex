---
name: simple-relion-eer-geometry-identical
description: SIMPLE and RELION render EER micrographs identically (no flip); SIMPLE xpos/ypos (0-based box corner) -> RELION centre = (pos + box/2) * bin
metadata:
  node_type: memory
  type: reference
  originSessionId: 851a7b77-9a45-48b1-b2b6-396c1a7b5829
  modified: 2026-11-03T03:26:58.260Z
---

Checked 2026-11-02 for RNase P r1312 polishing. The SIMPLE and RELION 4.0.1 source both place an EER electron at x = pos & 4095, y = pos >> 12. SIMPLE rmat(x,y) and RELION image(y,x) are written in the same MRC order.

SIMPLE window(): xpos/ypos are the 0-based lower-left corner of the box, so the centre is xpos + box/2. RELION rlnCoordinateX/Y are 0-based column/row centres. Conversion: centre = (xpos + box/2) * (SIMPLE smpd / RELION angpix). No y flip.

Measured on 4 micrographs: SIMPLE 2048-px intg against RELION 4096-px binned 2× correlate at 0.86-0.95 as-is, and about 0 flipped or transposed. SIMPLE micrographs live in /mnt/beegfs/elmlund/rnasep/preprocessing/motion_correct/. A RELION auto-refine that sits at one resolution for 4-5 iterations is normal: it holds while healpix order 2 converges, then refines the angular step (here 17.6 -> 16.0 Å at iteration 9). It is not a sign of bad coordinates.
