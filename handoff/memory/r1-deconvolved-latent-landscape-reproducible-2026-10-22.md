---
name: r1-deconvolved-latent-landscape-reproducible-2026-10-22
description: "f692st r1, 3 independent runs (script, port1, port2): the deconvolved flex-latent density and its poles REPRODUCE (port1-port2 BC 0.9988 > half ceiling 0.995; script 0.97; poles same maps within noise, 74-77 % same particles) though partitions do not; XD latent basin rule finds 2 kT/4 kT wells but misses 5 % states and unencoded states; masked per-run pole paths: 2 basins 0.38/0.62 in all runs"
metadata:
  node_type: memory
  type: project
  originSessionId: f0f8499a-d0a6-450b-9f93-483f98f27006
  modified: 2026-10-22T05:57:42.941Z
---

2026-10-22 overnight for the 2D-to-3D coordinating session (user bar: reproducible landscape and extrema, labels may be stochastic).
Scripts `routes/toy55/land_xdcompare.py`, `land_xdgauss.py`, `land_xdbasins2.py`, `land_latsep.py`, `land_runpath.sbatch`;
outputs `rnasep_klem/landscape/{latb_r1,runpath_*}`; write-up section 6d of `doc/implementation_notes/arm_motion_landscape_2026_10_21.md`.
- **Deconvolved latent density across runs.** Extreme deconvolution through SIMPLE's own noise model, read from
  `flex_pca_embedding.bin`; runs mapped by an affine fit, R^2 0.96-0.99.
  - port1 vs port2: Bhattacharyya 0.9988, above each run's even/odd ceiling of 0.995.
  - script vs ports: 0.97, under the script's ceiling of 0.98 (script is the odd run out).
  - Each density is far from Gaussian (0.73-0.76 against its moment-matched Gaussian), so the agreement goes beyond the covariance.
  - The 1 % and 99 % points along 4 axes agree within 0.03 sd (port1-port2) and 0.10 sd (script vs ports). SIMPLE's mixture K differs between runs (5 / 2 / 2) all the same.
- **Poles** (3 % tails of each run's own latent PC1):
  - always closed-site 0.90-0.91 vs open 0.92-0.99;
  - particle Jaccard 0.60-0.79 (74-77 % shared by all three);
  - pole maps agree within half-map noise (landscape_cmp 0.97-1.05 of the ceiling).
- **XD latent basin rule** (threshold 0.5 kT + half check):
  - swing30: continuum / flat / broad 1, 2 kT and 4 kT wells 2 (barriers read 0.84 / 0.92 kT), 5 % state missed;
  - simsites missed, because the focused latent encodes the arm at d' 0.63 only;
  - r1: 1 basin in every run. One latent basin is NOT evidence of a continuum.
- **Path rule between each run's own poles:**
  - **Masked** (10 % tails, morph only inside `rsmask_variance.mrc`, `land_morph2.py --mask`, `land_runpath2.sbatch`): 2 basins in all
    three runs.
    - Closed-pole side 0.37-0.39 / open-pole side 0.62-0.63 at lam 0.001, barrier 2.4 kT.
    - U-shaped p(s); the sparsest point is s 0.43 in all three.
    - p(s) BC 0.998-0.999; per-particle s Spearman 0.95-0.98 across runs.
    - Basin maps (`land_basinmaps.sbatch`, split at s 0.43) agree within half-map noise across runs (landscape_cmp 0.97-1.02). The
      closed-pole basin (30-34 %) keeps the OUTER arm (60-105 A) at the closed site; the inner arm is closed in both basins.
  - **TRAP, unmasked** (3 % tails): 25-49 % of each map "moves" (noise), so 60-70 % piles into one end bin. The path follows the latent
    only at |rho| 0.3. Always mask a morph between data-derived poles to the focus region.
- **Extreme LEAVES of 4 runs** (stage-B leaves, pooled stable PC1/PC2, `land_leafext.py`; lp 10 node mask):
  - **The 3 port runs reproduce their extremes:** port1-port2 within noise at all 4 poles (delta 0.76-0.90, rho >= 0.995); port3
    within 0.6-2.7x the noise.
  - **Script matches 1/4 poles** (axis-2 high). Its axis-1 extremes sit one leaf inward, and it has a fully closed leaf (closed-site
    1.00) that the ports lack (their most closed read 0.81-0.89).
  - **TRAP:** with the focus (variance) mask PC1/PC2 nearly tie (25 / 23 %), so single-axis "most extreme leaf" flips between runs.
    Use the node mask or the 2-D outer corners.
- **Paths between extreme leaves** (`land_leafpath.sbatch`, focused morph, shared node_r1/n.simple poses):
  - ports reproduce along axis 2: 3 basins 0.56-0.59 / 0.13-0.16 / 0.28 at both lam, p(s) BC 0.996-0.9995;
  - axis 1: 4 basins in port1/port2, 3 in port3, 2-3 at lam 0.01;
  - end bins hold 12-54 % (end pile-up), so read the interior wells, not the end masses.
**How to apply:** test landscape reproducibility on deconvolved densities and poles, not on partitions. Count states with
image-domain likelihoods (path rule, [[path-model-pose-marginalised-reweighting-2026-10-22]]), not latent densities
([[flex-latent-density-basins-no-power-2026-10-22]]).
