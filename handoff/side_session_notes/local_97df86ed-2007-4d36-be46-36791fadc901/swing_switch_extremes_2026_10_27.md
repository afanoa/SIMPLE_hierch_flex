# Swing vs switch + extremes (session swing-switch-extremes-58027e, 2026-10-27), interim

Cluster dir: /mnt/beegfs/elmlund/afan/swing_switch_20261027 (tools/, diag/). No native code changed; nothing committed.

What matters for other sessions:
- The node_motion map call CANNOT fire "swing" on real data. Its non-blend share nb scales with the score's precision rho:
  synthetic RNase P-scale rigid swing (maps160/swing30, tip 30 A, noise-free) nb 0.23 / 0.07 / 0.036 / 0.020 / 0.010 at
  rho 1 / 0.8 / 0.7 / 0.6 / 0.5; swing60 0.03 at rho 0.6. Real scores track motions at rho 0.5-0.7, so real swings sit in the
  switch range (0.004-0.025). A noisy switch and a noisy swing give the SAME blend-like bin series (and the same step share).
- Separability of the split axis = share of node rows with |post_hi - post_lo| < 0.5 ("mid"). qtoy2 discrete splits 0.12-0.27,
  OLE root 0.36 (cut reproduces across seeds, ARI 0.95-0.96); OLE depth-1 0.53-0.78 (cuts do not reproduce, ARI -0.02..0.25 on
  base trees). IgM continuous arms are ALSO 0.13-0.29 (big motions give confident posteriors), so mid measures per-particle
  separability, not discreteness; IgM is caught by the map call (nb >= 0.03) instead.
- On graded OLE axes the minpop score-tail poles reproduce across seeds (Jaccard 0.46-0.70 vs chance 0.03-0.05) where the
  cuts don't. Their cross-read extremes (+-1.35-1.87 sd vs leaves +-0.06-1.33) are inflated by shared image noise. On IgM truth,
  a node's 2k pole pair is about as extreme as the best leaf below it (median 0.57-0.59 vs 0.56 of the oracle).
- Refined-evidence tails (stage-B prob_state_evidence ev1-ev2, fixed N) are purer than raw tails on qtoy2 classes (closed
  0.96-1.00, conf43 0.99) and about equal on IgM. Write-up: this worktree's doc/implementation_notes/planned/
  hetero_tree_swing_switch_extremes_2026_10_27.md.
- Stage-A/B refinement drifts a continuum pole inward (OLE r2: seeds at -1.17 sd -> refined pole -0.75; seed 2's pole grew
  2.6k -> 10.4k, +1.48 -> +0.62) but purifies discrete minorities. Raw tails at RNase P r1211 overlap the refined closed leaf
  r12112 only 0.63-0.78 -> raw-tail poles would break the closed-state gate; pole membership needs refinement-based evidence.

## Native A/B running (2026-10-27 ~16:00)
- Build: swing_switch_20261027/install (login node fsitgl-xfer03p = Skylake, -march=native; 8c1b66f31 + swing_delivery_on_8c1b66f31.patch);
  fast gate 0. Key swing_delivery=yes (default no), needs site_layer=no. Every arm uses this one build.
- 12 trees in swing_switch_20261027/trees/<ds>_<off|on>_s<seed>: ole s1/s2 x off/on, qtoy, qtoynull, igm, olenull x off/on (s1).
  nodepar 16, steps on csbdevel,norm, SBATCH_CONSTRAINT x6342|x6740P. RNase P 692k on/off follows once a swing node has run cleanly.
- 2026-10-27 ~17:30, per the parent's request:
  - The seed-2 OLE pair is paused. Controllers 54413229 / 54413231 and pending step 54419547 were cancelled; the trees resume with
    tools/tree2.sbatch.
  - The RNase P pair (54422093 / 54422094 and their mean steps 54422136 / 54422137) was cancelled 35 s in. Its folders are in
    aside_20261027/. It relaunches fresh only after the seed-1 truth trees and OLE show a benefit.
- First live swing: olenull_on_s1 root, mid 0.999 -> poles r_pA / r_pB of 3,000 + middle r1 of 74,599 (92.6 %).

## Wound down (2026-10-27, all trees cancelled by the coordinator; no jobs of mine left)
- The rule acted at OLE r1 (mid 0.56) and r2 (0.679), qtoy2 r211 (0.745) and the OLE null root (0.999). Every IgM node was a switch
  (0.15-0.35). Fallback nodes have no site axis: qtoy null r/r2, OLE null r1, OLE r21, IgM r2112.
- qtoy2 r211 poles are not pure. The node is 88 % conf28 / 6 % holo (557) / 5 % conf49. Pole A (900) has holo 0.20 and conf49 0.12;
  pole B is the quality end of conf28 (snrq 1.51 vs 1.14). The key-off refined cut isolated holo better (leaf 1,203 at 33 % holo)
  but then quality-sliced the 92 % conf28 rest three times.
- OLE r21 (middle) fallback split is unrelated to the consumed r2 axis (AUC 0.46).
- Cost about 1,670 core-h, about 450 of it in idle 8-CPU controllers.
- Write-up section 8 has the details. The code is uncommitted in worktree swing-switch-extremes-58027e.

## Two approved cheap tests (2026-10-27 evening, about 38 core-h, no jobs left)
- OLE r2 replay (tree 6 r2 rows, key on, seeds 1/2):
  - Both seeds were graded (mid 0.82) on the same axis. Poles of 3,000 + a 23,992 middle; pole Jaccard 0.91 / 0.92 (chance 0.05).
  - This is an upper bound: same rows, so the seeds differ less than the saved trees did (0.60-0.62).
  - Along the motion series r2_s4 axis (own-pixel): pA -0.67 vs tree-6/7 extreme leaves -0.60 / -0.53; pB +0.50 vs +0.49 / +0.57.
- RNase P r22 band-split held-out: the open pole's extra reach is imprint. Held out, pole -2.91 = leaf r221 -2.87 vs random
  -2.60 +- 0.22. Same band: pole -7.1, leaf -5.
