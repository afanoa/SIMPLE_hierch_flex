# flex -> refinement (f2r) replication: status and traps

VERDICT (16:15): replication holds on the task's criteria, but NOT worth building. Calibrated, the prior switches itself off on
real data (r131 full alpha 0 from it 3; qtoy2 gain halves 15k -> 60k: it-1 NMI +0.040 -> +0.018, held-out +0.077 -> +0.032).
Frozen at the between-stages calibration (chain P, the only form the hard gate allows) it holds r131 enrichment 1.42-1.43 but
by overriding held-out data (-0.004..-0.025 nats/particle by it 4) and anchors the olenull split (-0.015). Every frozen-pose
chain ends below the labels it started from (qtoy2 0.599-0.615 vs 0.655-0.658), so the proxy cannot show a gain over
production stage B. Nulls: qtoynull sites t6 alpha 0 all its; olenull alpha 0 after a +0.0000 it-1. No Fortran written;
the fflex call site in simple_hetero_tree.f90 is free for the refine -> flex session.
Per-run tables: runs/<key>/f2r_<tag>/f2r.log (P chain in *_p tags), flex weights saved as f2r_weights.npz in the p runs.

LATER THE SAME DAY (all under runs/<key>/, scripts f2r_fuse.py, f2r_regref.py, f2r_regref2.py, run via f2r_run2.sbatch):
- Two-distribution posterior (f2r_fuse_*; per-state learned Gaussians over [latent log-ratios, objective gaps], cross-half, shared
  covariance F / block-independent FI / objective-only OQ): worse than the calibrated product on qtoy2 15k (F -0.044 vs C +0.077
  nats/ptcl held-out; OQ -0.11..-0.25: the objective's gaps are not Gaussian within a state), equal on r131 at it 1 only, and
  overconfident + negative on both nulls (olenull F -0.035 by it 4, mean max post 0.72 vs objective 0.59). No off-switch on nulls.
  Its one virtue: label stability (qtoy2 full FI NMI 0.632 -> 0.631; r131 full F 1.396 at it 4) as an anchor to the input labels.
- Discrete-vs-continuum (inode fixtures, held-out latent gain at stage-B labels): f50 0.000 (alpha 0), cont +0.0006, null +0.0005,
  f15 (15 % minority) +0.028 stable. Does NOT detect a continuum; responds to a noisy minority-state reference. Line closed.
- Latent-span regularised references (f2r_regref_*: each state's deviation from the all-map projected onto flex modes + site
  factors; R = full projection, RH = projected inside flex's 15 A band, raw beyond; FSC weighting recomputed on the PROJECTED
  deviations): stops the drift on r131 full (objective-only enrichment 1.377/1.380/1.378/1.367, labels kept 0.76-0.77, where
  plain O goes 1.336/1.393/1.302/1.187 and C 1.379/1.423/1.422/1.287); at it 1 equals the calibrated prior with no prior.
  Costs signal on qtoy2 full (R 0.574, RH 0.588 vs O 0.615, C 0.633 at it 1; span holds 60-68 % of the deviation, raw deviation
  FSC 0.2-0.4 beyond the flex band). NULL HAZARD: projection makes a null's deviation look reproducible (olenull deviation FSC
  0.22 -> 0.95, qtoynull 0.36 -> 0.94; the site factors in the span come from the same arbitrary split), held-out best beta
  0.25-0.5 (+0.012-0.019 nats overconfident), posteriors harden 0.63 -> 0.84, drift not reduced (labels kept 0.67 -> 0.52).
  Calibration still gives alpha 0 on nulls (the latent adds nothing), but the references would pass a reproducibility gate.
- Raw-FSC weighting of the projected deviations (array 55003117, f2r_regref2_*, RF full / RHF hybrid): restores calibration
  (olenull overconfidence 0 at every it, peak +0.002; qtoynull +0.0001 at it 1 but +0.012 by it 4 under drift; r131 full best
  beta 0.75-1, posteriors 0.66-0.70 vs O 0.71 -> 0.91) at the cost of the stability: r131 full RF 1.382/1.377/1.369/1.356
  (projected-FSC variant 1.367 at it 4; O 1.187; C 1.287), r131 15k RF 1.391 -> 1.295 (BELOW O 1.316; C 1.468), qtoy2 full RF
  0.566 -> 0.557, RHF 0.581 -> 0.574 (O 0.615 -> 0.560, C 0.633 -> 0.599). The FSC weight dials shrinkage against stability;
  no setting beats the calibrated product on truth anywhere, and C itself fades on production references. LINE CLOSED.
  Untested and only testable in production: RELION-style sum over poses for the state evidence (logsumexp over the state's
  block at T_i in assign_rows / gather_state_row) instead of the (state, pose) argmin.
- DONE (prob_state_marg A/B, build merge_sm = claude/state-eval tip + the switch, fast gate 15/15; stage B rerun from recB/b.simple
  with the node's own keys, then state_split; runs math_reworks/latent/ab_sm/<node>_{argmin,marg}): r131 closed enrichment
  argmin 1.450 vs marg 1.466 (production, build_bc, 1.375), 94.7 % of labels shared, both split 2 ways like production; qtoynull
  and olenull one child in both arms (94.0 % / 98.3 % labels shared). The rule moves 2-6 % of the labels and the gain (+0.016) is
  far below the build/seed difference (production vs argmin 1.375 vs 1.450, 89 % shared). NOT worth keeping. The switch is
  uncommitted in the tree-master-merge worktree (seven files; branch claude/prob-state-marg if ever committed).
- Posterior temperature by held-out prediction (f2r_beta.py, runs/<key>/f2r_beta_{fs3,ffull}): at each particle's stage-B pose
  + 64 random rotations (median 2 deg, 90 % 6 deg; shifts fixed) every state's cross-half reference scored on the two
  checkerboard views; q_A ~ exp(beta l_A) over the 65 poses x K states predicts view B. Joint beta* (band A 15 A / B 10 A
  identical): qtoy2 0.30, r131 0.20 (15k and the full 25k node), qtoynull 0.05, olenull 0.10. Gain of beta* over 0.3: qtoy2 0,
  r131 +0.008..0.013, nulls +0.08..0.17 nats/particle; beta = 1 costs 0.34-1.46 nats/particle everywhere. Decomposition:
  pose-only beta* 0.10-0.20 (the pose posterior is the overconfident part: mean max over 65 poses 0.13-0.38 at beta 1 vs
  0.02-0.05 held-out, near uniform); state-only beta* 0.75-1.5 at the stored pose (calibrated or under-confident). One joint
  beta is a compromise -> 0.2-0.3. The full-band argmax leaves the stored pose for 99.5 % of particles: within +-6 deg the
  per-particle pose evidence is nearly flat at these SNRs.
  CORRECTED VERDICT (parent session, 2026-11-06 evening): master's prob_beta tempers ONLY the pose draw (assign_rows /
  gather_state_row take the state by argmin, then draw the pose within that state's block at prob_beta * T_i), so the
  optimum that matters is the POSE-ONLY one: qtoy2 0.10, r131 0.15-0.20, qtoynull 0.05, olenull 0.10 (gain over 0.3:
  qtoy2 +0.20-0.22, r131 +0.017-0.018 nats/particle). "0.3 is at the optimum" holds for the joint posterior only. The
  parent launched gate trees at prob_beta=0.15 (arm pb15, jobs 55003886-89) against mm2 (0.3) and mb1 (1.0): pb15 LOST
  the closed corner (0.91/0.53) and merged qtoy2's conformer 1/2 pair (NMI 0.748), so the local held-out optimum is not the
  tree's optimum; 0.3 stands. Caveat: the
  pose-only estimate here tempers a local +-6 deg neighbourhood of random rotations; master's draw is over the state's
  whole direction block, where far directions are strongly disfavoured, so the local estimate calibrates the spread near
  the mode. A per-iteration Fortran estimator still looks unnecessary if pb15 behaves (beta* 0.1-0.2 is stable).
- LATENT SEEDING of stage A (f2r_seed_labels.py: K-state mixture in the flex + site-factor latent at INPUT poses, K =
  production's seed K, cross-half assignment; f2r_replay*.sbatch: the node's own seeds -> stage A -> fold -> recB -> stage B
  -> split calls with the mm2 build, project and starting maps swapped, SIMPLE_SEED per replicate; runs ab_seed/<key>_<arm>_s<seed>,
  readouts readout/labels.txt (stage B), labels2.txt (split under the gate env), arm.txt (r131 maps)). Harness: at
  production's seed the replay reproduces production's stage B bit-for-bit (qtoyflex, r131 via write_states; qtoy2/olenull
  replicates within 0.1-5 % of labels); the split step needs SIMPLE_SPLIT_CONTRAST=maps exported (gates/<arm>/split_env.txt),
  without it state_split's choose rule flips to the EM groups and reallocates ~10k particles (an earlier "split not
  reproducible" note was this). Child/fixture projects carry the parent's stacks and prg=selection refuses them: a
  write_states program (tree-style state write, ptcl3D + ptcl2D, no stack check) was added to merge_sm for r131.
  Seeds before refinement: latent beats site_em (qtoy2 NMI 0.571 vs 0.309; r131 r1 closed enrichment 1.18 vs 1.02;
  qtoyflex 0.138 vs 0.078). Finished stage B + split, production / site_em replicate / latent s1 / latent s2:
    qtoy2 NMI 0.656 / 0.657 / 0.649 / 0.650; closed F1 0.55 / 0.55 / 0.61 / 0.61; holo 0.94 all; K* 5, children = states
    qtoyflex NMI 0.168 / 0.168 / 0.172 / 0.172; closed F1 0.17 all; holo 0.62 / 0.62 / 0.63 / 0.63; 3 children all
    r131 r1 closed enrichment 1.247 / 1.245 / 1.023 / 1.021; maps intermediate in both arms (best closed-arm 0.87 vs 0.88,
      none >= 0.95: closed separates at r2 in the mm2 tree, not r1); 2 children all
    qtoynull, olenull: 1 child in every arm and seed.
  VERDICT: latent seeding does not beat site_em seeding on the finished stage B and split: ties on both toys (the one
  gain, qtoy2 closed F1 +0.06, is purity 0.44 vs 0.38 at recall 1.0), loses the closed-row enrichment on r131 r1 at both
  seeds, nulls unchanged. Stage B is nearly seed-invariant and converges to the same basin from seed sets that share only
  NMI 0.03-0.38: the initialisation carries nothing through stage A/B that site_em's does not. No Fortran port.

# Global chart coordinate (f2r_global.py on mest_charts.py's readers and composition; runs math_reworks/latent/global/<tree>/
# report.txt, global.npz; compare2_*.txt and compare_full692_*.txt)
Trees: mm2_r131 (seed 2), mm2_r131_s1, pb15_r131, mm2_qtoy2, pb15_qtoy2, mm2_qtoyflex; 692k: rnasep_klem/stree_port/full692_{m33,base92,base80}_s1.
(a) MINORITY, closed vs rest, cross-half LDA AUROC (d'), global frame / root chart: r131 s2 .689 (.67) / .616 (.40); s1 .671 / .615;
  pb15 .671 / .616 (the composed frame reaches the root-latent supervised ceiling ~.69 the root chart alone misses); qtoy2 .994 (4.5) /
  .982 (3.0); qtoyflex .822 (1.3) / .785 (1.1); 692k m33 .752 (.94) / .772 (1.02), base92 .765 / .772, base80 .819 (1.31) / .772: a wash
  on deep trees. The closed-enriched leaves do NOT stand apart: r131 sibling d' 1.35-1.51 for enriched and plain pairs alike; 692k m33
  closed pair r1111221_s1/_s2 d' 1.12 (AUROC .79, overlap .74), 86 of 112 reported pairs further apart (lower quartile). qtoy2's closed
  leaf is the most isolated sibling (d' 2.7 vs 1.3-1.8).
(b) REPRODUCIBILITY (same particles; leaf-partition NMI | global canonical corr. 1-3 | variate-1 extrema overlap, chance .05 | root chart):
  r131 s2 vs s1 .19 | .89 .88 .87 | .66/.68 | .998; r131 0.3 vs 0.15 .34 | .94 .93 .92 | .80/.77 | 1.000; s1 vs 0.15 .17 | .87 .86 .84 |
  .67/.64 | .998; qtoy2 0.3 vs 0.15 .82 | .99 .99 .98 | .94/.80 | 1.000; 692k m33 vs base92 .45 | .96 .86 .57 | .41/.29 | 1.000;
  m33 vs base80 .24 | .45 .33 .30 | .22/.13 | 1.000; base92 vs base80 .24 | .44 .33 .30 | .23/.14 | 1.000. The reproducible landscape
  is the ROOT CHART (identical fit for the same input, .998 across seeds); the composition adds the children's detail at the cost of
  reproducibility, down to .45 across recipes on 692k. Closed-enriched leaves across trees share Jaccard .27-.45 (r131) and .02-.06 (692k).
(c) MERGE / STOP: on pb15_qtoy2 global d' orders the one truth-equivalent pair (1.17) below the distinct ones (>= 1.8, except a junk leaf
  at .35): AUROC .97; mm2_qtoy2 and the toy have no equivalent pairs. On 692k m33 a latent-distance rule would MERGE the closed pair
  (d' 1.12, lower quartile), as the 32-unit merge-up did and against the hard gate; that merge-up's unions bear no relation to latent
  distance (its united pairs span d' .2-54, the sibling pairs it left apart sit at .2-1.0).
(d) DENSITY MINIMA on the split axis, composed frame vs the ROOT CHART (one affine frame, no seams), valleys (depth @ sd) per half:
  r131 s2 composed .87/.89 @ -.02 vs root NONE (1 mode); r131 s1 none vs none; pb15 .91/.92 @ .00 vs NONE; qtoyflex continuum .80/.81
  @ -.09 vs NONE; qtoy2 (discrete conformers) 1.00 @ -.2 vs root .59/.61 @ -.14/-.15 (tree's cut -.01), reproducible within .01 sd.
  => the composed coordinate is a PIECEWISE-AFFINE ATLAS: particles of different children reach the root frame through different
  maps (R^2 .1-.6), so a valley along a split direction, a cross-node leaf d' (the 20-54 values on 692k) or an extremum can be a SEAM.
  Seam-free statistics: within-node (sibling) pairs, and anything computed in the root chart. The "s2 bimodal / s1 unimodal" contrast
  on r131 was a seam. Real valleys exist only where states are discrete (qtoy2), where the tree already cuts correctly.
VERDICT: the global coordinate answers neither question. It does not find the closed minority on real data (d' <= 1.3, closed
  leaves no more isolated than any sibling pair, the 692k closed pair among the closest), and it cannot decide merges or stops (it
  would merge the closed pair; its between-node distances and valleys are chart seams). What reproduces (root chart, .998-1.000) is
  already an output; the composition is less reproducible than the leaves' nearest-equivalent on 692k. Nothing to port. Possible fix,
  untested: fit the child -> parent maps on overlapping particle sets (LLC-style) so that neighbouring charts agree at their
  boundary; without overlaps the atlas cannot be made continuous.

Owner: the f2r side session. Files: `f2r_calib.py`, `f2r_run.sbatch`, `f2r_{smoke,list,full}.txt` (local
`scripts/hetero_tree_proto_latent/`, cluster `.../math_reworks/latent/`). Outputs: `runs/<key>/f2r_<tag>/f2r.log`, `f2r_summary.json`.

Jobs (batch2): smoke 55002116 (passed), 15k array 55002117 (s3 repro of calib2 t2, s2 on 4 nodes, a3 = stage-A labels+poses on
4 nodes, s3 qtoynull backstop for calib2 t6), full-node array 55002118 (4 nodes, sites basis).

Traps found:
- `prep_crop` accepts any cache whose SHAPE matches: another 15k subsample with `--cache fp_<key>.npy` silently loads the seed-3
  planes. f2r uses its own caches `cache/f2r_<key>_<tag>.npy`; calib2/calib3 caches untouched.
- calib_flexprior2.py rebuilds the references of all arms from the CALIBRATED arm's labels, so its "obj" row after iteration 1
  is one objective-only E-step at the calibrated chain's references, not an objective-only chain. f2r runs independent chains
  C (calibrated), O (objective only), S (stickiness to the input labels, matched at it 1 to C's kept-label share).
- In the sites arm the calibrated (beta, alpha) often has alpha/beta = 1, i.e. the same hard labels as the plain product (1,1);
  calibration then only changes the softness.

Corner job 55000383 (seed 3, 15k, stage-B labels/poses), sites basis: qtoy2 C 0.623/0.620/0.617/0.615 vs obj@C 0.584/0.579/0.576/0.572,
held-out +0.077..+0.062 nats; r131 enrich C 1.392/1.414/1.441/1.468 vs obj@C 1.332/1.420/1.459/1.462, held-out +0.0073 -> +0.0004;
olenull alpha 0 from it 2 (it 1 alpha 0.75 at +0.0002); qtoynull flex basis alpha 0; qtoynull sites (q 49) still fitting > 40 min.
