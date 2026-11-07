---
name: averages-route-real-data-2026-10-14
description: "heterogeneity from kernel-expanded class averages without particle 3D alignment (user: \"pretty significant\"); spliceosome motion visible from motion-bin maps; NOT broken by orientation bias; balance=cavg in solve3D_cavgs (reseed30)"
metadata:
  node_type: memory
  type: project
  originSessionId: 34b66a02-03d2-42df-91f8-915afbf80120
  modified: 2026-10-14T19:54:05.522Z
---

2026-10-14. Route: 2D classes -> cls_expansion subcls_avg=kernel -> solve3D_cavgs on the sub-class averages (its
map2ptcls gives every particle a pose) -> flex_pca on the even+odd averages (flex_cavgs=model neigs 10 lp 12
flex_view_centre) -> readout. Scripts on batch2 /mnt/beegfs/elmlund/afan/toyexp/pv (fastchain = box-160 multi-part,
fullchain1n = full box one node, trajbins = equal-population bins along latent components, reconstructed from averages
and particles; flextests = tied parent poses / equal_occ path states).

- Spliceosome 10180 (46,783 selected, 34 classes, ~1,400/class; 30 sub-classes per class): flex on the averages has a
  10-dim reproducible motion subspace with a flat spectrum; GMM/kcenter give one state + view-coupled fringes, so read
  it as a continuum. Motion-bin particle maps (fast160/trajbins_*) show the foot missing and SF3b moving, smeared. User:
  "we can get heterogeneity out of expanded class averages without even aligning particles". Particle-flex 5 states
  are only x1.2-5 enriched in the kernel sub-classes (RNase P 39-76x): small classes limit the expansion.
- NOT (CCR4-NOT/ribosome, box 420, 1.3 A): strong orientation bias; 12,600 expanded averages made solve3D_cavgs
  unusable. Built solve3D_cavgs balance=cavg (+ nclust, nsample; reseed30): units = parent classes, view groups by
  cluster_cavgs_into_groups, class-balanced fractional sampling per iteration, a final all-average alignment pass, final
  map from a balanced subset. Box-160 test ran but was cancelled by the user; one-node lp 30-12 runs (1 and 3 states)
  launched. User prefers ONE node (nparts unset, nthr 32): the 8-part runs lose most time to per-iteration sbatch.
- Untested so far: head-to-head of averages-flex vs particle-flex bases on the same data; sharpening by E-step.

Related: [[averages-flex-pv-toy-validated-2026-10-14]], [[averages-unmixing-overnight-2026-10-14]].
