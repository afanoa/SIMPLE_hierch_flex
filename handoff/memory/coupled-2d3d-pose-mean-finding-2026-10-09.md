---
name: coupled-2d3d-pose-mean-finding-2026-10-09
description: "toySE local toy: cls_expansion subclasses split on VIEW (12-deg cone > 15-deg swing); tangent nuisance columns do not fix it; subtracting the consensus projected at each member's OWN pose (pose mean) makes subclasses split on state (NMI 0.01-0.04 -> 0.41-0.70, purity 0.4 -> 0.6-0.8); SIMPLE slice convention euler2m([-e1,e2,-e3]), shift phase +s"
metadata:
  node_type: memory
  type: project
  originSessionId: a6c2687a-b3cd-4070-9eaa-6163aa554f8d
  modified: 2026-10-10T03:29:39.928Z
---

2026-10-09 night, autonomous local work (no VPN): scripts/cls_expansion_proto/subclass_pool/{run_toyse_local.sh,
pool_place.py, view_tangent_toy.py}, data ~/toyse_local (SNR 0.1, true poses) and ~/toyse_snr03 (SNR 0.03, refined
poses), note doc/implementation_notes/coupled_2d3d_flex_toy_findings_2026_10_09.md.

- At true poses / SNR 0.1 flex_pca's 3D latent is NOT a blob (15-NN state purity 0.99, k-means 0.96); the toy is the
  3D-favourable regime, not RNase P's. The user: "for RNase P the 2D embedding is far richer" - the structured test needs
  SNR ~0.02, refined poses, tight classes, 30 A motions.
- cls_expansion subclasses on this toy split on view (NMI 0.27 vs state 0.07; 48 classes = 12-deg half-angle cones).
  In-plane, first- and second-order view-tangent nuisance columns do not remove it.
- **Pose mean works:** per member subtract s * P(R_i) mu3D (consensus at the member's own pose, carried into the class
  frame with the same in-plane alignment), fit the class PPCA on the residual -> subclasses split on state. This is the
  coupled 2D-3D model's first concrete form (3D supplies the view term, 2D spans conformation). In-plane tangents on top
  hurt. Fortran port = per-member projected mean in flex_class_planes (needs ptcl3D poses + vol1).
- **SIMPLE projection convention** for Python: image FT = slice at loc = [h,k,0] @ euler2m([-e1, e2, -e3]) (euler2m =
  r3 r2 r1 with SIMPLE's rotmat), shift s -> exp(+2 pi i k.s/box); centred FFTs; volume [z,y,x]. Calibrate with
  simulate_particles snr=1000 ctf=no oritab=<known> (snr=0 gives NaN). The in-plane class alignment that sharpens the
  class mean was rot=-e3_2D, shift=+(x,y)_2D.
- solve2D (this build's abinitio2D) writes the classified project in 1_solve2D/toy.simple (master copy untouched);
  cls_expansion writes subclass ids into ptcl2D 'class' and the parent into 'cluster'; refine3D needs sigma2 in its own
  dir (use objfun=cc for toys); reconstruct3D needs pgrp.
Related: [[toyse-subclass-pool-test-setup-2026-10-09]], [[cls-expansion-vs-flex-pca-latent-diagnosis-2026-10-09]],
[[pose-state-coupling-toy-verdict]].
