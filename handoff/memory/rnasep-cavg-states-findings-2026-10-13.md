---
name: rnasep-cavg-states-findings-2026-10-13
description: "RNase P states from 2D-expanded averages: solve3D_cavgs single-state works (orig 112 and expanded 1,252), multi-state de novo gives balls; refine3D_states_cavgs flex seeding is 65% random and the march collapses to ~2 states; empty-state reseed fix uncommitted in ~/SIMPLE_master1013; coupled view calibration breaks the weight publish after merge"
metadata:
  node_type: memory
  type: project
  originSessionId: 34b66a02-03d2-42df-91f8-915afbf80120
  modified: 2026-10-13T19:53:23.739Z
---

2026-10-13, batch2 runs under `/mnt/beegfs/elmlund/afan/rnasep_coupled/r3s_cavgs/` (newest master c6e344c6b build
`~/SIMPLE_master1013`, snapshots snap/c6e344c6b and snap/c6e344c6b_reseed).

- solve3D_cavgs nstates=1 reproduces the RNA from the 112 original 4_abinitio2D averages (today's code, Oct-8 build,
  and today's code with the 2025 a.txt params all fine) AND from the 1,252 expanded sub-class averages (CTRL_EXP).
  nstates=4/8 de novo on the expanded averages gives discs / noisy balls (S3D_K4/K8): the multi-state path, not the
  inputs. The user's 2025 reference map: `afan/reextract/7_abinitio3D_cavgs` (abinitio3D_cavgs, single state; same
  cavgs_iter044 + frcs + 112-class selection, byte-identical to atifao/reextract/4_abinitio2D).
- Audit 2025 abinitio3D_cavgs vs solve3D_cavgs: stage 1-2 search changed from per-average shc_smpl (deleted
  2026-10-06) to prob_neigh with prob_athres 90 from iteration 1; gauref Gaussian refs always on (params%gauref
  ignored: bug); distributed path new. Shared-memory solve3D_cavgs on 112 cavgs takes 2.5-5 min vs 34 min with nparts.
- The docked late split (multivol_mode docked / split_stage) was removed from master before Oct 8.
- refine3D_states_cavgs needs nparts (shared memory rejected). On the expanded averages its flex seeding cannot meet
  the population floor (2 qualifying clusters, 1,635 of 2,504 rows randomised); the march then collapses states:
  K4 segfaulted on a 24-row state; K8 consolidated to ~3 then 2 groups (1,253 / 1,188 rows) every block. Fix added
  (uncommitted, `~/SIMPLE_master1013`): policy l_reseed_empty_states refills EMPTY states at block boundaries
  (worst-fitting averages of the largest state, halves together) and rebuilds volumes; near-empty states still pass.
- Coupled placement with SIMPLE_FLEX_VIEW_CALIBRATE=yes failed at "state weight outside [0,1]; publish" after the
  two-gate merge folded a view-clustered state (calibrated weights summed > 1); rerun without calibration.

- refine3D_states_cavgs ALREADY docks: every row is posed against the project's consensus map (here
  rnasep_states/24_refine3D_auto/recvol_state01.mrc, lp 10) before any split; the user asked "can't you just go in
  with a consensus volume" - the answer is that it does; the split is the weak part.
- Label-seeded runs (infile, pooled 3D-latent k-means, snap c6e344c6b_reseed3): K3L finished (rows 125/1203/1176,
  FSC0.143 13.8/5.3/7.1 A). K4L/K8L died mid-stage when one state drained to a single row (K4 state 2 725->87->57->10->1;
  K8 state 7 reseeded 34->314 then 207->15->1): rec3D_pcg_strategy line 653 skips a state with NO rows but aborts on a
  state with rows in only one half. My first explanation ("drift to one halfset") was wrong; user called it out.

- Drain mechanism (logs): every iteration the state whose map has the worst FSC loses its averages next
  iteration (K8: four states at 26-117 A emptied in ONE iteration). User framed it as "the noise model doesn't
  support conditional cavgs": global sigma ignores 1/members per average, and the euclid score treats the
  state map as exact. K4 at lp 20 (lpstart=lpstop=20, job 53828294, r3s_cavgs/K4L20): stage 1 still drained
  2 states (maps 21-23 A); after the block refill all 4 maps were 5.6-8.7 A and held through stages 2-4
  (734/645/395/730 rows) -> equal map quality, no drain.
- Fix written 2026-10-13 (~/SIMPLE_master1013, uncommitted, snap c6e344c6b_reseed5): line-653 one-halfset skip;
  sigma_pop=yes (row sigma2 x (1/pop)/mean(1/pop), euclid_sigma2%read_groups, pop carried on cavgs work-project
  rows); state_refs=deviation (PCG distributed post-loop: refs = row-weighted consensus + deviation x FSC of
  even/odd base-pair deviations; final rec untouched). refine3D_states_cavgs defaults both on.
- Per-state PCG conditioning = distinct VIEWS, not averages (sub-classes share their parent's view): maxits_pcg=30
  converged for 1,252 avgs/112 views, diverged for 317/85 ("lost positive-definiteness", pcg:2910). Crowther pi*D/d
  (D=200): ~30 views at 20 A, ~105 at 6 A. The expansion IS view-balanced (5-16 kept sub-classes per parent, user
  pointed this out); the pooled-3D-latent seeds are not (K4: 28/112 parents span all 4 states; K8: 0/112). Fix idea:
  view-balanced state assignment (every parent contributes to every state; possible to K=5).
- K8 segfaults (K8R?, K8D2, K8L20): FSC never >= 0.143 -> get_resolution fsc0143 = 0 -> postprocess calls
  guinier_bfac(hp, 0.) (lplim < 5) -> get_find(0) overflows nk -> fit_straight_line reads past the slice. Guarded in
  simple_commanders_volops.f90 (lplim > 0) and guinier_bfac (lp <= 0 returns 0); snap reseed8.
- Coupled 692k placement v2 (15 states, 260k/250k + 13 small): user judged the maps "junk".

**How to apply:** seed cavg states from meaningful labels (coupled per-sub-class states or pooled 3D latent), K ~3,
not flex on averages; run solve3D_cavgs single-state first (dock), then split. Related:
[[rnasep-2d-subclasses-pooled-vs-particle-3d-latent-2026-10-13]], [[coupled-2d3d-reconstructed-tree]],
[[rnasep-24-refine3d-auto-quarter-unrefined]].
