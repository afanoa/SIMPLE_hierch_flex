---
name: refine3d-lp-set-scores-merged-map-trap
description: "TRAP (found 2026-10-21, again 2026-10-29): refine3D with lp= scores every particle against the MERGED map in both halves (read_mask_filter_refvols); hetero_tree node refinement is self-referential; held-out = 2 passes with volN = other half; solve3D too (every stage sets lp, 2026-11-06)"
metadata:
  node_type: memory
  type: project
  originSessionId: 9643ce3d-fc7a-4dc1-860a-59db10867e7c
  modified: 2026-11-06T13:29:16.630Z
---

Found 2026-10-29 by the micelle session (hetero_tree micelle fixture, soft multi-state prototype).
- Code: simple_matcher_refvol_utils.f90 read_mask_filter_refvols ~433-446, branch `params%l_lpset .and. .not. l_nonuniform_mode`
  reads the merged vols(s) and copies it to vol_odd; l_lpset = cline%defined('lp') (simple_parameters_phases.f90:591).
- Who passes lp=: hetero_tree stA0 / stage A / stage B (lp 15/15/10), the stage-B hard classifier, any refine3D run given lp.
  Without lp (automatic lp) the halves stay separate and blend_lowres_eo_for_registration averages < max(4 shells, 20 A).
- CORRECTION 2026-11-06 (IgM DJ session): solve3D is NOT an exception. Its stages set lp (stage manifest l_lpset=yes in all 5 stages
  of dj_ws2 r300k), and with filt_mode=nonuniform_lpset (lp promoted) l_use_merged_nu_ref copies the merged recvol_stateNN_nu_filt.mrc
  into both halves (master1030 build too). solve3D half maps are not independent: r300k s6 halves read identically to 3 decimals.
  Measured: the s6 halves agree .82-.90 at 20-9 A while s6 vs EMD-13921 and vs our consensus is ~0 beyond 25 A (16 A "FSC" is imprint).
- Consequence: classification inside a node can lock onto noise (8 soft or hard iterations from any start grew confident splits on
  a null), and groups built that way share noise across halves, so even the split test accepted them (micnull z 22, guard 1.025).
  The tree's protection is state_split's own cross-half EM re-test on stage-B groups.

- First found 2026-10-21 (tree-as-prior task, lrec/rflex_coord XHALF=1 swaps were no-ops): no FSC filter on the lp path (cosine
  low-pass at lp only). Held-out scoring recipe from then: keep lp set, pass volN = the OTHER half's map times its FSC, run refine3D
  maxits=1 twice over ALL particles (pass E keeps the even rows, pass O the odd rows; a per-half selection dies in the sigma2
  bootstrap, [[sigma2-bootstrap-group-empty-half]]). Measured then: own-particle contribution small for score-PCA coordinates (0.725 vs 0.734).
- Even with true halves (soft3 build), refine3D's own start maps reach the matcher merged at iteration 1: pass vol1..K from a separate
  reconstruct3D of the start labels for a clean cross-half run.

**Why:** a "cross-half" scoring claim is false whenever lp is set; the shrink arm logged an e/o difference FSC of exactly 1.000.
This trap was already in memory (2026-10-21) and was rediscovered on 2026-10-29: check this note before building any cross-half scheme.
**How to apply:** any test of gold-standard scoring inside refine3D must force true halves (the micelle session's soft3 build does,
under SIMPLE_SOFT_STATES only); never read node-level half maps as independent. Related: [[micelle-fixture-pose-scatter-verdict-2026-10-29]],
[[per-window-refine3d-shared-reference-inflates-reproducibility]].
