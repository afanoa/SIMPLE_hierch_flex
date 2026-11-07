---
name: averages-flex-pv-toy-validated-2026-10-14
description: "flex_pca on kernel sub-class averages gives 3D states on a preferred-view toy (model, neigs 10, lp 12); kernel enrichment needs view-tight 2D classes; mkcavgproj.sh smpd trap; one-transform map docking"
metadata:
  node_type: memory
  type: project
  originSessionId: 34b66a02-03d2-42df-91f8-915afbf80120
  modified: 2026-10-14T16:36:57.319Z
---

2026-10-14, batch2 /mnt/beegfs/elmlund/afan/toyexp/pv (scripts there: pvtoy.sbatch, flexmaps_pv.sbatch, mapcmp_pv.sbatch,
avgflex.sbatch for real data, mkcavgproj2.sh, pv_readout.py; local copies in the session scratchpad pv/).

- The user's route: kernel cls_expansion -> solve3D_cavgs (1 state) -> even+odd averages as ctf=no rows with the
  sub-class pose, class = parent (mkcavgproj2.sh, SMPD env) -> flex_pca flex_cavgs=model neigs=10 lp=12
  flex_view_centre=yes state_placement=gmm npreimages=K min_neff=5 rec_states=yes.
- pv3 toy (oc20 volumes 0/+20/-20 deg at 50/30/20 %, 80 view clusters of 4 deg): particle accuracy 0.897 vs ceiling
  0.902; each state map closest to its own truth; difference maps r 0.68-0.85 at 20 A. neigs=5 (default) is too few
  (latent R2 for composition 0.26 -> 0.54 on pv2); flex_cavgs=ls with neigs 10 worse (0.79).
- 13 % minority (pv2): reaches ceiling accuracy but only 3.3k of 15.6k particles; limit is the expansion.
- Kernel sub-classes are state-enriched only when 2D classes are view-tight: uniform-view toys gave 3.4x, the
  preferred-view toy 29-79x (RNase P 39-76x).
- TRAP: rnasep_coupled/flexcavgs/mkcavgproj.sh hard-codes smpd=1.464 (RNase P); every toy flex-on-averages run before
  2026-10-14 ran at the wrong pixel size (band ~39 A, mask wider than the box). Use mkcavgproj2.sh.
- TRAP: docking each state map to the truth separately puts alignment noise in difference maps (r 0.45); dock the
  consensus once and apply that transform to all maps with volops (needs mskdiam; convention e1 e2 e3 + shift as
  reported by dock_volpair).

- Model variants tested 2026-10-14 (snap reseed29; SIMPLE_FLEXAVG_PAIR=1, SIMPLE_FLEXAVG_BFAC=1, noise precision
  in ptcl3D corr_cart, sub-class in cluster -- import_particles DROPS unknown keys like subcls/wnoise): noise per
  average, even/odd pairing, per-row B-factor, per-view offset removal (avg_vcentre.py) -> no gain on the toys (flex
  sits at the expansion ceiling 0.90). Untested on real data with a wide neff range.
- What works: solve3D_cavgs already maps the averages' poses to the particles (spproj%map2ptcls), so the flex
  labels go to particles directly -> reconstruct3D per state (pv3 diff maps 0.68-0.85 -> 0.89-0.92) -> score_states
  E-step x3 (pv3 0.897 -> 0.950, pv2 13 % minority 0.897 -> 0.978, minority recall 0.95; pv2 diff map 0.44 -> 0.81).
  Scripts: ptclstates_pv.sbatch, estep_pv.sbatch.

Related: [[averages-unmixing-overnight-2026-10-14]], [[flex-on-expanded-cavgs-plan-2026-10-13]].
