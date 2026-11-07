---
name: rnasep-periphery-ladder-2026-09-29
description: RNase P on-axis bin ladder (3/6/12/24 states per region axis from the BW_n1 cache) and region-masked e/o FSC; e/o FSC cannot see motion blur; ~12 reproducible positions per axis; region 1 softest in every state
metadata:
  node_type: memory
  type: project
  originSessionId: a0f08b9d-781d-471c-b0a6-698a8a8afe30
  modified: 2026-09-29T19:47:49.405Z
---

Runs `/mnt/beegfs/elmlund/afan/rnasep_regions/LAD_ax{1,2,3}_n{3,6,12,24}` (states-only resumes of the BW_n1 3-column
block embedding, `state_axis=q npreimages=n rec_states=yes SIMPLE_COV_DECONV=0`, ~3 min per 12-state arm with 10 workers);
scorer `~/rnasep_ladder/regfsc.py` (region-masked even/odd FSC + end-to-end and adjacent-bin difference reproducibility),
table `ladder_logs/regfsc_table.txt`.

**Findings (2026-09-29):**
- Region-masked EVEN/ODD FSC is the wrong instrument for "I cannot see base pairs in the periphery": the consensus and
  the refine3D_states finals hit Nyquist (2.9 A) inside every region, including the unreadable one. Motion blur is
  reproducible in both halves. Sharpness needs a sharp target (toy) or a map-feature metric, not e/o FSC.
- Each region axis is a clean single-region motion: in-region share of the end-to-end difference 0.67-0.71, other
  regions 0.03-0.08, in-region e/o cc 0.998. Adjacent-bin differences reproduce at ~0.87 with 12 bins but ~0.5 with 24
  -> about 12 distinguishable positions per axis; finer slicing adds noise.
- Region 1 (G10cov_B3_1) is the softest region in EVERY state (e/o 0.5 threshold ~6.0-6.6 A vs 4.2-4.8 for regions 2/3
  at ~130k particles per state), including the states sliced along its own axis: slicing one axis does not remove its
  residual softness -> multi-axis motion within the region or per-particle noise, not a missing slice.
- Whole-molecule GMM states (BW_n1) are worse than on-axis slices at matched neff.

**How to apply:** do not judge peripheral sharpness by local e/o resolution. The toy3 line ([[toy3-body-frame-test-system]])
is where sharpness vs the truth is measured; the RNA transfer needs a label-free sharpness readout. Related:
[[rnasep-combined-block-embedding-states]], [[flex-state-maps-selection-bias]].

**Body-frame transfer to the RNA (2026-09-29 night, `rnasep_regions/rnabody_s{+1,-1,0,ctrl}`, script `~/toy3/rna_body.sbatch`
-> `rna_body.py`):** region-1 generators from recvol_state01 x G10cov_B3_1, latent from BW_n1, autog from LAD_ax1_n3 (scale
0.067: the RNA latent coefficients are LARGER than the state-difference coefficients, opposite to the toy). Calibrated
motion is SMALL: |omega| median 0.69 deg (p90 2.5), |t| 0.4 px, rotation covariance rank 1 (1.49/0.07/0.03 deg). The
three sign variants give identical maps (4.4/6.8 A e/o whole, region 1 5.0/6.3) -- worse than the consensus pair, so a
like-for-like control (unchanged poses through the same reconstruct3D, `rnabody_ctrl`, CAL=0) was submitted. Pose
import route that keeps the stack segment: `oriops oritab=<txt> oritype=ptcl3D outfile=<ptcl3D-only .simple>` (no
projfile) then `replace_project_field projfile=<copy of full project> projfile_target=<that> oritype=ptcl3D` (mkdir=yes
forced: the result lands in 1_replace_project_field/). `oriops projfile=... oritab=...` IGNORES the table (builder reads
oritab only without projfile).

**CORRECTION + sub-region result (2026-09-30 early):** every RNA arm up to `rnabody_ctrl` had UNCHANGED poses:
`replace_project_field` throws "Inconsistent # of ptcls" because the ptcl3D-only source has no ptcl2D rows, and the
script fell through to the copied project. Working import = `~/toy3/patch_project.py <project> <new table> <orig
table>`: in-place numpy patch of the ptcl3D segment (binoris header 20x5 int64; segment 6; 50 float32 per record;
1-based e1=8 e2=9 e3=10 x=27 y=30), verified against the print_project_field export (max diff 0). With real poses
(region-masked e/o at 0.143/0.5, control 4.4/6.8 whole, region 1 5.0/6.3): region 1 body-frame 4.7/7.3, 5.1/6.4;
sub-regions of region 1 (`sub1/R1sub_B3_{1,2,3}`, co-variation split of the 10 region-1 block columns, within/between
ratio only 1.6; 3-column BW_n1 latent, autog) motions 0.8-1.9 deg median: maps unchanged or worse (sub-region 2, the
largest fitted motion, 7.0/9.4 in its own mask). VERDICT: the latent-visible rigid motion of region 1 or its
sub-regions is <= ~2 deg and body-frame reconstruction does not sharpen the RNase P periphery; the softness is not a
latent-visible rigid body motion (non-rigid, or per-particle disorder below the 16 A band). The toy line stays valid.

**Two decisive tests, 2026-09-30 (both negative):**
1. Focused local refinement of region 1 (`FOC_a{5,10}_lp7`: refine3D maxits=1 refine=neigh athres 5/10 trs 2/3 lp=7
   vol1 = consensus x G10cov_B3_1, automsk=no; score the *_unfil halves, the filtered ones read ~5.4 A everywhere):
   region 1 e/o 5.5-5.6/6.1 vs control 5.0/6.3 (no gain), core 3.9 -> 5.0-5.1 (pose drift from the region-only
   reference). No alignable coherent per-particle motion in region 1 at the 7 A band.
2. Finer-band region-1 latent (`HI_r1`: pcg_mskfile=G10cov_B3_1 neigs=4 lp=10; NOTE smpd_target is overridden by
   lp/2.5 -> box_crop stayed 64, kstop=23): 4/4 columns reproducible, but the ladder is the SAME picture as the 16 A
   fit: end-to-end 0.75 in-region, ~12 positions along the axis (adjacent-bin 0.86 at 12), region 1 the softest in
   every state (4.6-5.0/6.1-6.3). Nothing new above 16 A.
**Verdict:** the RNase P periphery blur is not a latent-visible or locally-alignable rigid motion; it is either
non-rigid within the region or per-particle disorder. Stop pursuing rigid per-body poses on this dataset; the toy3
machinery ([[toy3-body-frame-test-system]]) is validated for cases with real multi-axis rigid motion.
