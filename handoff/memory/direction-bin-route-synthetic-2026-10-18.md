---
name: direction-bin-route-synthetic-2026-10-18
description: direction-bin 2D route on the 4 synthetic sets (synth_bins_route.sbatch): per-bin sub-classes carry the heterogeneity (tethered2 pure per site, arms 0.40/0.52) but 4 variants of a K-state EM over the averages with known directions (avgs_states3d.py) stitch nothing across views; failure modes and what is left
metadata:
  type: project
---

2026-10-18 00:30-01:30. Pipeline (toy55/synth_bins_route.sbatch, coupled snapshot 9eaef62df for cls_expansion pose-mean
mode): dir_bins.py on lat_$MODE/ptcl3D_r0.txt (100 bins, 90th-pct radius ~11 deg) -> import class=bin -> make_cavgs ->
cls_expansion ncls=6 lp=20 vol1=consensus pcg_mskfile=support sphere -> avg_stack.py -> avgs_states3d.py (K-state EM over
the 600 sub-class averages with KNOWN directions: Fourier-slice insertion/extraction on a c=80 grid, cross-half E-step)
-> particles inherit the state -> label_truth_readout.py.
Raw per-bin sub-classes (label_truth_readout, 1 = no info / tethered 0 = pure): arms 0.40 (2D classes gave 0.46),
arms2 0.52, tethered 0.21, tethered2 0.00 (median 0.56!). The 2D half works.
Stitching variants, all = no information (arms 0.97-1.0, tethered 0.73-0.84): (1) plain EM on the averages -> clusters
by view coverage; (2) deviations from the consensus slice + per-bin balanced (Hungarian) assignment -> states' deviation
volumes reproduce across halves (0.92-0.97) but are systematic per bin (CTF/scale/direction error), not conformational;
(3) bin-centred deviations (within-bin contrasts only) -> still reproducible (0.8-0.9) non-conformational structure,
probably residual view splitting inside the 11-deg bins. Geometry traps fixed on the way: FFTs must be origin-centred
(ifftshift) or tilted slices carry (-1)^l phases (corr 0); image arrays are [y,x], volumes [z,y,x]; the matched
convention is R^T, e1 +, no y flip (consensus slice vs average 0.84-0.91 from shell 2). Real-space check: projections
by scipy affine_transform reproduce the bin averages at 0.93.
Conclusion (same as [[routes-2d-to-3d-states-tests-2026-10-14]]): the within-view signal is clean; cross-view
assignment by a volume likelihood fails three ways; next would be proportions from the averages' own latent (unmix +
sync, 0.74) or tighter bins (<5 deg) with sub-classes forced to share the bin's pose mean. Stopped here for the night.

**Variant 5, flex_pca on the bin sub-class averages (02:00-02:20; synth_avgflex.sbatch, build
SIMPLE_master1013/snap/c6e344c6b_reseed27 = the only build with flex_cavgs + flex_view_centre; import needs kv/cs/fraca
even with ctf=no; coordinates land in the run dir, no *_flex_pca subdir):** flex_cavgs=model flex_view_centre=yes lp 12
neigs 10, full/even/odd fits. GMM placement collapses (95 % of sub-classes in one state, as on RNase P). The latent
itself (avgflex_readout.py): cross-validated R2 of the per-sub-class truth (arm angles / site occupancies) from
z1..z6 = 0.00 on all four sets (arms2 arm4 0.25 in-sample only, negative cross-half), while components reproduce
even/odd at up to 0.84: reproducible, non-conformational structure again. k-means placement: no information.
Standing conclusion: five stitching variants fail; the heterogeneity is in the sub-classes, the cross-view link needs
within-view contrasts propagated between NEIGHBOURING bins (ManifoldEM-style), not a global fit over all views.

**Variant 6, neighbouring-bin propagation (ManifoldEM-style; neighbor_propagate.py, 2026-10-18 11:xx):** bin-centred,
band-passed (60-12 A) sub-class contrasts; neighbour pairs (< 30 deg, 6 nn); in-plane frame calibration on the bin means
(corr 0.83-0.88); cross-half Hungarian matching of contrasts (matched 0.18-0.26 vs unmatched -0.03). Local truth check
(propagate_check.py: centred truth profiles of matched pairs): agreement mean 0.09 (arms) / 0.03 (tethered) / 0.19
(arms2) / 0.23 (tethered2) vs ~-0.03 chance, 19-35 % of matched pairs > 0.5 -> real but weak, improves with domain
mass. Tree propagation of K=6 labels -> no truth (product-space toys have no consistent global 6-labelling). Next:
permutation synchronisation over all edges (SYNC=1, added), and the Fab-mass arms3 set where the local agreement
should be stronger.
Sync result (20-deg neighbourhoods, 140-148 edges): local truth agreement of matched pairs 0.30 (41-43 % > 0.5) on
arms2/tethered2, but permutation synchronisation (top-6 eigenvalues 2.6-4.3 only) gives global labels with no truth
(arms2 0.90-0.98, tethered2 0.71-0.84). Verdict: K global labels is the wrong object for a product space (each bin splits
along a different arm/site); what can be propagated is a per-domain contrast coordinate. Parked here.
