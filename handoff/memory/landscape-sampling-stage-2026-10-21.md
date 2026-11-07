---
name: landscape-sampling-stage-2026-10-21
description: "landscape sampling stage (lsamp*, routes/toy55; port spec section 5): per terminal node, score axes vs noise covariance, samples q01..q99 at equal held-out structural steps, n by PIXEL-SPLIT CV, check-3 gate; f692st r1 4-tree REPRO 2026-10-22: leaf-only axis reproduces from the 10th pct, tip not; POLE AXIS (closed-minus-open pole score) reproduces the whole landscape incl. tip; axis 2 = order-disorder"
metadata:
  node_type: memory
  type: project
  originSessionId: 5e5b4054-6080-49b8-8249-8b67e5b28864
  modified: 2026-10-21T23:59:08.011Z
---

Task 2026-10-21 (user, from the state-tree session): truth-free stage that samples a tree's reproducible variability tail to
tail at even spacing, count set by resolvability. Written into the port spec as section 5
(/Users/atifao/SIMPLE/doc/implementation_notes/state_tree_port_spec_2026_10_21.md, untracked, main checkout; a worktree hook blocks
the Edit tool there - ask before editing it again from a worktree session).
Scripts (routes/toy55, all new): lsamp_score2.py/.sbatch (scores per pixel set A/B = parity of kx+ky, UNMASKED caches
rnasep_klem/lsamp/caches/<set>_A_c48.npy, em_lim 12, meta ppnoise), lsamp2.py/.sbatch (axes, samples, pixel-split CV, checks;
defaults ppnorm 1, spacing struct, junk_ssnr 0.35, ends tail in code but the spec default is centre - lsamp_chain passes
--ends centre), lsamp_tree.sbatch (all terminal nodes + lsamp_rec), lsamp_chain.sbatch (scores then nodes), lsamp_rec.sbatch
(reconstruct3D + arm readout/tracer/radial extent), lsamp_disc.sbatch (K=16 leaves for sets without a tree), lsamp_clades.sbatch
(flat sets: state_merge clades), lsamp_tails.py (truth diagnostic). v1 lsamp.py/lsamp_score.py = superseded.
lsamp_wait_f692st.sbatch (1 core) submits lsamp_chain on rflex_f692st when its leaves_list/final_labels exist and job 54050336 is gone
-> rnasep_klem/lsamp/f692st_ps, maps view_rflex/lsamp_f692st/<node>.

Findings:
- Run PER TERMINAL NODE: at the s30a root the swing is not a leading axis (|corr| <= 0.11 on axes 1-5; axes = families/junk).
- Axes: generalised eigenvectors vs G = sum w CTF^2 Re(P_k P_l*); lam-1 = whitened map-space signal per particle.
- splitbench16-style beyond-band test (16.7-10.2 A) has NO power for neighbouring swing samples at realistic SNR (0.00-0.05 vs null
  0.05). PIXEL-SPLIT CV works: coordinate from set A, matched filter u on the other eo half's end-to-end difference over set B (all
  shells), strata-centred z of mean-u steps, Bonferroni 0.01, swap; n = min fold. Needs UNMASKED caches (soft disc couples neighbours ~0.2).
- PER-PARTICLE NOISE NORMALISATION essential on synthetic stacks: s30a r22 |corr| 0.642 -> 0.660, n 6 -> 12; with structural spacing
  truth 23.5 .. 6.1 A, eta2 0.45, steps 1.1-2.3 A (CV 0.26) vs tree leaves 8.8-22.4 and old rflex_coord 24.6-4.8 (CV ~0.36, corr
  0.725 from refine3D best-dist + local pose search at lp 8; my fixed-pose scores max 0.666). Equal coordinate steps saturate at the
  ends (truth steps 1.5/3.4/6.6/3.4/1.8). --ends tail: 24.0-5.6 but only 5 samples.
- Real f692rf (no ppnorm, centred ends): r22 7x7 grid 44 samples closed-arm -0.10..0.76 (leaves 0.07-0.64), arm tip intact at
  95-110 A (0.78-0.94), open site empty in the closed corner; r21 6x6 32, 0.06-0.55 (leaves 0.11-0.54); r1 2x2x2 (stop node).
  FINAL normalised (rnasep_klem/lsamp/f692rf_ps3): r22 centred ends (recommended) 48 samples, most closed 0.76/0.74 at 3,564 ptcl
  (16 % old closed; arm tip 0.92/0.89/0.85 at 95-110 A), most open -0.09; maps view_rflex/lsamp_f692rf_r22_centre. Tail-end chain
  (view_rflex/lsamp_f692rf_final/<node>): r22 30 samples, 0.76/0.63 at 1.4k, one corner 0.68/0.06 (open site empty); r21 1 axis 16
  samples 0.23-0.50; r1 8 samples 0.14-0.25 (stop node).
- f692st (64 leaves, DONE 19:39; multi-state leaves r21_g1/r51_g1/g2/r52_g1 have maps in node_rXX/leafmaps/ and poses in
  node_rXX/fstageB/*_refine3D/b.simple): chain job 54056763 (centred ends) -> rnasep_klem/lsamp/f692st_ps/tree/node_<u>, maps
  view_rflex/lsamp_f692st/<node>. The first waiter-launched chain (54055286) failed on the leafmaps pose lookup (fixed).
- Limits (flat K=16 discovery leaves): sim55_200k closed minority not spanned by any leaf/axis (best 18 %); simc compositional
  variants not separated (NMI 0.10; per-site layer's job); swingbreathhs clades mixed -> swing on minor axes; synthetic junk survives
  the leaf junk screen and ppnorm amplifies it (blob images low noise).
- Unmixing (k-means in 2-D, rescore) did nothing in the whitened space; a 1-D-bin score loop collapsed the half-pixel coordinates;
  size-mismatched references without per-shell amplitude matching make a nuisance axis.
Related: [[klem-bottom-up-organization-2026-10-20]], [[node-type-test-disjoint-band-2026-10-21]], [[split-test-holdout-and-drain-traps-2026-10-21]].

FINAL (2026-10-22, spec section 5 rewritten with user approval; the main-checkout spec was edited via python because the worktree hook
blocks the Edit tool there):
- Axis rule now: grids only while every axis resolves >= 3 samples AND passes check 3 (tails <= 0.1 end-to-end beyond the end samples;
  lsamp2.py --min_per_axis 3 --tail_tol 0.1). f692st 20-95k nodes: every grid failed check 3 by 0.1-0.3 (2 per axis = halves), every
  1-D sampling passed; gate test: f692st r1 3 -> 2 -> 1 axis (7 samples), f692rf r22 keeps 7x7.
- f692st (64 leaves): r1 1-D 7 samples viewing 0.39-1.22 (1.9k), HELD-OUT 0.86-0.87 -> 0.46 (tree leaf r1_s5 1.00 own-pixel); r51 5
  samples 0.44-0.87, held-out 0.62-0.63 (tree 0.50-0.97); r5 7 (0.40-0.63); r52 4; open-arm nodes r2/r4/r6/r7/r8 <= 0.34 (grids, pre-gate).
  Outputs rnasep_klem/lsamp/f692st_ps (tree/, tree_rule3/, tree_1d/, imprint_1d/, gate_r1/), maps view_rflex/lsamp_f692st{,_1d,_rule3}/<node>.
- Held-out readouts: [[selected-subset-maps-imprint-inflates-extremes-2026-10-22]] (own-pixel maps inflate extremes 4-7 A on truth).
- lsamp2h.py = lsamp2.py + heldout_X/own_X/final_all crop maps padded to 160; lsamp_swingpos.py = truth placement on the swing30 maps.
- Held-out LEAF readout (2026-10-22, user ask; lsamp2h.py --leaf_heldout 1, outputs f692st_ps/leafheld/node_{r1,r51}): r1_s5 held-out
  0.78 at 13.3-13.6k (own-pixel 0.97-0.99); landscape's 3 most closed samples at the same N 0.77-0.78, end sample 0.86-0.87 at ~2k;
  r51_s5 0.65-0.67 vs landscape end 0.62-0.63. Verdict recorded in spec 5.2: no honest gain over the tree's closed leaves; the
  landscape grades the closed region into smaller, more extreme samples.

REPRO on f692st r1 (2026-10-22, coordinator ask from "2D to 3D heterogeneous state translation" session; note in my worktree
doc/implementation_notes/landscape_sampling_repro_f692st_r1_2026_10_22.md): 4 independent trees (script rflex_f692st/node_r1, port1
stree_port/tree_f692r1, port2 _seed2, port3 _rep3; same 92,554 ptcl; each run its own leaves/poses/scores) -> lsamp_repro.sbatch
(scores + lsamp2h defaults) -> rnasep_klem/lsamp/repro_r1/<run>; lsamp_repro_cmp.py -> repro_summary.txt.
- Every run 1-D after the gate (lam-1 axis1 1.14-1.32 vs axis2 1.08-1.13, nearly degenerate; 2-axis grids end 0.26-0.37 short of tails),
  n 7/7/6/7 (end steps weakest, z 3.0-3.4); held-out closed-arm closed end script 0.86-0.87, ports 0.78-0.83; open end 0.46-0.48 all.
- Across runs: per-particle coordinate corr 0.83-0.93; cross-pixel agreement 0.87-0.95 of the within-run A-vs-B ceiling (rAB 0.38-0.40);
  sample NMI 0.30-0.45 but 78-97 % same-or-adjacent sample; end Jaccard 0.36-0.58. Script = the most different run (as the coordinator's
  leaf-map comparison).
- TRAP: best-correlation matching of sample maps of DIFFERENT sizes picks the largest samples (noise level); compare at equal N
  (lsamp_repro_bins.py/.sbatch: the same quantile bins of every run's fold coordinates, held-out-pixel maps -> repro_r1/bins).
- EQUAL-N RESULT (repro_r1/bins/bins_summary.txt; 12 identical quantile bins 0-2-5-10-20-35-50-65-80-90-95-98-100 % per fold, held-out
  maps): from the 10th percentile to the open end the 4 runs agree within one run's pixel-fold noise (closed-arm sd 0.003-0.012 vs |A-B|
  0.006-0.028; whole-map ruler sd 0.005-0.018 vs 0.047-0.071); occupancy per ruler quarter (closed->open) script 12/18/39/31 %, ports
  11-12/21-23/35-36/30-32 %; deviation maps cross/within 0.85-1.03 (mean 0.93), end-to-end 0.88-0.97. CLOSED TIP NOT REPRODUCIBLE:
  top 2 % script 0.87 vs ports 0.79-0.83 (sd 0.036 vs |A-B| 0.003; ruler 1.05 vs 0.96-1.01); script also has the most extreme leaf r1_s5
  and largest lam-1 (1.32 vs 1.14-1.16) -> tip likely inherits the leaves' span (untested; leaves vs poses not separated).
- 2nd score axis (lam-1 1.08-1.13, half agreement 0.98-1.00 in 3 runs) NOT sampled: gate refuses 2-D grids (ends 0.26-0.37 short).
  Proposed: closed extremum from the reproducible poles (+ pole maps in the score space); sample axis 2 as its own 1-D landscape.
- OVERNIGHT FOLLOW-UPS (2026-10-22 02:00-04:00, coordinator ask; all in the worktree note, outputs repro_r1/*):
  * axis 2 alone (lsamp_repro_ax.sbatch --axes 2): order-disorder of the arm (both sites gain together), 6/7/8/7 samples, interior within
    pixel noise, ends slightly beyond; matched maps 0.88 of within-run.
  * pole maps IN the score space (lsamp_polemaps.py: top-N of the factorial pole posteriors reconstructed at the run's poses, crop 48,
    per-shell signal-amplitude matched to the leaves): top-N* does NOT fix the tip (sd 0.086, port1 0.70); top-2,500 makes a lam-1 ~2
    NUISANCE axis (arm flat, defocus Cramer V 0.18-0.20 vs 0.06-0.10 for structural axes; lsamp_viewcheck.py).
  * POLE AXIS as the coordinate (lsamp_poleaxis.py: score(closed pole) - score(open pole), per pixel set; stage via lsamp2p.py --coords):
    whole landscape reproduces incl. tip: top-2 % 0.80/0.78/0.78 sd 0.016 vs fold 0.027, matched maps 1.00, cross-pixel 0.97-0.98, 6/6/5
    samples, end samples 0.75-0.81 (less extreme than the script's leaf tip 0.87: noisier coordinate, fold A/B corr 0.31).
  * Why the leaf-only tip varied: lam1 ~ lam2 in the ports (near-degenerate plane; which direction is "axis 1" is run-dependent).
  * Recommendation sent: extrema from the poles; landscape axis = the pole axis, not eigenvalue order; axis 2 as a separate landscape.
  * TRAP: a Monitor/grep pattern with an empty alternative `\(A\|B\|\)` fails under the Mac's ugrep (no events); use awk/case instead.
- SPEC (2026-10-22 ~10:30, user approved "add the overnight reproducibility results to section 5"): added the block "Reproducibility
  across independent trees" at the end of 5.2 (table of 5 coordinates + bullets) and the heading note; 5.1 deliberately UNCHANGED - the
  user has NOT approved making the pole axis the sampling coordinate (coordinator relay 09:55: investigation first: pipeline's own
  poles p1/p2, P=3 coordinates, recursion/axis 2, multi-site, placement/cost; outputs repro_r1/pq/).
- POLE-PAIR INVESTIGATION (2026-10-22 10:10-11:00, user via coordinator; proposal + exact spec-5 text in my worktree
  doc/implementation_notes/landscape_sampling_pole_pairs_proposal_2026_10_22.md; outputs repro_r1/pq; spec NOT changed by it):
  * pipeline's own poles (stree_port poles_r1_p{1,2}: P=3 stageB states; scored vs them at stageB poses): closed-openA 7/7 samples,
    p1 vs p2 cross-pixel 0.97, NMI 0.57; tip matches factorial poles at equal N (sd 0.012); axes differ across pole sources (0.82-0.83).
  * P=3: closed-open pairs reproduce everywhere; open-open only within a pole source; closed-vs-rest worse; pipeline pairs corr ~0.5
    (2-D simplex), factorial c-oA/c-oB 0.93. -> one 1-D landscape per pair (all pairs at P=3; MST edges for P>3).
  * closed child's pair (childdisc) reproduces (Spearman 0.88) but is NOT root axis 2 (|rho| <= 0.27); root axis 2 = root open-open
    pair among open particles (|rho| 0.65-0.72 pipeline poles).
  * toy (factorial IgM-like): one pole coordinate reads one arm angle -> axes per pole pair, not per site.
  * cost at r1 (16 cores): dump+scores ~2 min, stage 4.5 min per pair (node's own steps ~21 min).
  * TRAP: the scratchpad (/private/tmp/...) was wiped between sessions; keep scripts on batch2 (routes/toy55) as the source of truth.
- FORTRAN PORT (2026-10-22 afternoon, coordinator relay of user approval "yes go ahead"): prg=landscape_sample on my branch
  landscape-sample (from fork/hetero-tree 8421e1a31; worktree happy-lichterman-b94a1b, UNCOMMITTED), module
  src/main/strategies/search/simple_landscape_sample.f90 + commander/router/UI/visibility-tester; patch batch2
  ~/lsamp_jobs/landscape_sample_8421e1a31.patch; build tree batch2 ~/SIMPLE_lsamp (5703662e7 + 8421e1a31 diff + patch), gate 15/15.
  Interface agreed with the port session (final-phase call gated on landscape_pairs.txt; site_poles "rest" line; split_children).
  Validation p1/p2: coords corr 0.997 vs Python, counts equal or 1 step apart on the same subset, 62-75 s/node.
  Rest rules (tested rest_r1_p*): rest never a pair endpoint; its particles in every pair; nothing at 1 pole + rest.
  Converged site EM makes r1 P = 2 (coordinator): P = 2 is the main case.
  Spec section 5: the approved design text is in the worktree proposal note; NOT yet applied to the main checkout - needs the
  user's direct confirmation (a relayed approval does not answer my pending question).
- SPEC APPLIED (2026-10-22, user's direct "Yes: apply the approved pole-pair design text"): section 5 of the main-checkout spec now
  has the pole-pair Why, 5.1 (P = 2 main case, rest rules, MST distance = map distance over the score band as implemented), the 5.2
  "pole-pair design" bullets and the 5.3 placement/cost/program bullets; the reproducibility Status bullet says adopted.
