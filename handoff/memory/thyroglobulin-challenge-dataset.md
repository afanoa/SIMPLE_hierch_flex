---
name: thyroglobulin-challenge-dataset
description: "Flatiron 2023 heterogeneity challenge thyroglobulin set (bioRxiv 2025.07.18.665582, PMID 41280101): where the data is, what is in the tar, and the bar other methods set (trimodal GT, ~9 A per-volume FSC)"
metadata:
  node_type: memory
  type: project
  originSessionId: 230e38d3-a42d-4cd7-ab95-1ddcf5c084f4
  modified: 2026-10-06T18:02:29.244Z
---

Downloaded by the user 2026-10-06 to the box's LOCAL disk (`fwl-c142339:~/Downloads/`), copied the same
day to `/mnt/beegfs/elmlund/afan/thyroglobulin/cryoem_heterogeneity_challenge_2023_first_dataset.tar`
(27,460,925,440 B, size verified), extracted there. SIMPLE project (2026-10-06, box master
1f890abc5): `thyroglobulin/tg_ds1/tg_ds1.simple` -- 33,742 ptcls, cryoSPARC poses + CTF verified
against the star; recipe `thyroglobulin/tg_import.sh` + `tg_star_convert.py` (A->px, e1/e3 wrapped
to [0,360), per-particle smpd, stack-dir prefix). `tg_ds1_failed_import/` = empty first attempt.
When the user says "it's on the box", look in `~/Downloads` sorted by mtime first.

Tar = dataset 1 = EXPERIMENTAL only (no ground truth): one `.mrcs` 33,742 x 448 px float32 at
1.073 A/px (K2, 300 kV, Cs 2.7, Q 0.1); RELION-3.1 `.star` with cryoSPARC poses, defocus, and
`rlnOriginXAngst` shifts in ANGSTROM (convert, see [[simple-star-shifts-are-pixels]],
[[simple-star-pose-import-traps]]); unsharpened consensus `thyroglobulin_volume.mrc` 448^3; aligned
`thyroglobulin_model.pdb` (6SCJ-based). 660 kDa (pseudo-)C2 homodimer; consensus 3.1-3.3 A.
**The provided consensus map is EXACTLY C2 about z** (180-deg self-corr 1.0000) and SIMPLE's Cn axis is z:
run consensus reconstruction/refinement with `pgrp=c2` (the user expects it). My C1 runs (2026-10-06):
bootstrap_rec3D 4.29 A sphere-FSC / 3.53 A density-FSC; C1 refine3D_auto got WORSE (3.82, 3.85 A, poses
flipping, mean 9-12 deg) and was cancelled. C2 bootstrap+refine3D_auto = job 53712270 (`tg_ds1/c2.sh`): 3.97 -> 3.53 A sphere-FSC0.143, 4.01 A FSC0.5,
3.48 A density-masked; final map `tg_ds1/4_refine3D_auto/recvol_state01*.mrc`. 2026-10-07: flex_pca npreimages=3 pgrp=c2
(PCG, pcg_mskfile=4_refine3D_auto/nu_envmask3D_state01.mrc) job 53720431 in `tg_ds1/flex_pca_c2/`, and refine3D_states
nstates=3 flex=yes pgrp=c2 job 53720432 in `tg_ds1/r3Dstates_c2/`; each dir holds a COPY of the project + sigma2_state.bin
(flex reads the project's relative sigma2 path; check `sigma_whitened=T`, no white-noise fallback line).

Dataset 2 (cisTEM-simulated, same poses/defocus, GT = 3 Gaussians on MD-PC1 at -150/0/+150 with
populations ~41.5/17/41.5 %). OSF 8h6fz's Dropbox add-on is DEAD (every file 404s, even the repo's
demo link; the repo deleted data/fetch_data.sh). Everything lives at
`https://users.flatironinstitute.org/~ccb/cryoEM_heterogeneity_challenge_2023/`:
`..._second_dataset.tar` (37.8 GB) and `submissions.tar` (262 GB, 2026-09-07: 22 synthetic + 16
experimental submission .pt at 7.2 GB each, GT_submissions/{averaged,sampled}_gt.pt, no Salted Caramel,
no `maps_gt_flat.pt`/`metadata.csv`). The dataset-2 tar (38,764,789,760 B, downloaded intact) has
particles at 448 AND 288 px (+ a star and consensus volume for each) and NO GT latents either.
Both wget -c'd into the same beegfs dir on 2026-10-06
(logs `dl_*.log`, ~30 MB/s). Tar members can be listed remotely by 512-B Range reads of each header.
Scoring code github flatironinstitute/Cryo-EM-Heterogeneity-Challenge-1
(`run_svd_pipeline`, `run_map_to_map_pipeline`, `run_distribution_to_distribution_pipeline`).
Submission = 80 vols, 2.146 A/vox, box 224, aligned to the provided reference, + populations.

Bar to beat: every method recovers the PC1 motion; no blind round-1 submission recovered the
three modes (only non-blind Cherry 2 / Chocolate 2 / Vanilla 2 did); best per-volume FSC0.5 vs GT
~9 A (Cookie Dough 1), worst ~18.5 A. Fig S9: ab-initio poses shrink the motion and fake a bimodal
ensemble because of the C2 pose/state degeneracy, which is our [[pose-state-coupling-toy-verdict]] problem.

C2 flex_pca (pgrp=c2, which the cnga1 rule says never to pass) gave a weak latent (d_signal=1, het-res 16 A,
states = k-means halves of one blob + an outlier shell). C1 SYMMETRY-EXPANDED flex_pca (2026-10-07, job 53720944,
`thyroglobulin/tg_symexp/flex_pca_c1/`, 67,480 = 2 x 33,740 hard-linked copies, e1+180): latent z1 is ANTISYMMETRIC
(corr(copy1,copy2) = -0.968), z2/z3 symmetric (+0.78); states 2/3 (27.9k/31.7k, 4.37/4.33 A) are C2 MATES
(corr(s2, rot180z(s3)) 0.956 vs corr(s2,s3) 0.852; copies cross-tab 10.2k/10.4k off-diagonal) = one asymmetric
structure in its two registers; state 1 (7.9k, 7.6 A) = outlier shell. Open question: real one-protomer motion vs
pose nuisance (rotation about an axis perpendicular to z is also C2-antisymmetric) -- test the s2-s3 difference
against rigid-rotation derivatives before claiming biology. Per-state C2 refine3D_auto of refine3D_states (C2,
flex init; 11.4k/11.8k/10.6k) = jobs 53721232-34 in `tg_ds1/r3Dstates_c2_auto_state0{1,2,3}/`.

