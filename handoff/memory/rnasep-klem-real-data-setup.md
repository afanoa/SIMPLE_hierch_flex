---
name: rnasep-klem-real-data-setup
description: "How the kernel-latent EM (klem.py --meta) runs on the real RNase P set: 10 stacks denoise_raw_particles_part01..10.mrcs (box 160 @ 1.464 A, 932,496 particles), meta_all.npz (poses, CTF, eo, flags: 692,392 refined, z692 latent), consensus 25_refine3D_states/rec_final_state01.mrc, labels rnasep_softem/states.txt (932k rows); job run_klem_rnasep.sbatch; first gates 2026-10-15 12:40"
metadata:
  type: project
---

Real-data inputs for `klem.py --meta` (`scripts/cls_expansion_proto/toy55/run_klem_rnasep.sbatch`, out
`/mnt/beegfs/elmlund/afan/rnasep_klem/klem<KTAG>`):
- stacks: `/mnt/beegfs/elmlund/rnasep/22_denoise_project/denoise_raw_particles_part01..10.mrcs`, 160x160, mode 2,
  1.464 A, 932,496 particles in total (the stk segment of rnasep_coupled/full692k_2d/rnasep.simple says ctf=yes;
  10 "ctf=flip" strings also occur -> klem decides signed vs |CTF| at 16-12 A by correlation with the consensus).
- `rnasep_coupled/basis_span_all/meta_all.npz`: e1 e2 e3 x y dfx dfy angast eo stkind (1-10) indstk (1-based) for all
  932,496; flags==1 = 692,392 refined (24_refine3D_auto; the other 25.7 % never refined); z692 (692392 x 16) = the
  earlier flex latent of the refined subset; lab692 / parent692 = cls_expansion sub-class / parent labels.
- consensus for calibration: `rnasep_states/25_refine3D_states/rec_final_state01.mrc` (160 @ 1.464).
- external labels for sanity readouts: `rnasep_softem/states.txt` (932,496 rows, S1/S2/S4 soft-EM states; klem maps
  them onto the refined subset; 0 = unlabelled).
- gates 12:40: 53887799 (50k, random init, klem_gate_rand), 53887800 (50k, npz init from z692, klem_gate_npz);
  full 692k run next if the gates behave (expect ~2 h at 16 cores).

FIRST REAL-DATA RESULT (gate npz, 50k, 13:25): CTF model auto; 92 % of particles on the grid, mean max r 0.56.
Nuisance: NMI(cell, tilt bin) 0.006, rotation 0.010, defocus 0.004; |Spearman(z, tilt or defocus)| <= 0.10 -> the
latent is NOT view or defocus. Structure: NMI(cell, cls_expansion sub-class) 0.141 vs NMI(cell, parent class) 0.031;
within-parent NMI 0.218 (median over 111 parents, view fixed) -> the latent separates sub-classes of the same view
without any 2D information in the loop; soft-EM states: NMI 0.042, S1-like state at z1 2.84 vs S2/S4 at 1.95/2.08
(S1 separates along z1, consistent with earlier findings). Full 692k run launched 13:30 (rnasep_klem_full, npz init,
out rnasep_klem/klem_full_npz).

Replication from a RANDOM start (gate rand, 50k, 13:35): nuisance NMI 0.004 / 0.007 / 0.002 (tilt / rotation /
defocus), NMI(cell, sub-class) 0.129 vs parent 0.028, within-parent 0.207 (111 parents), S1 separates by 0.6 cell on
z2 (weaker than the npz start's 0.9 on z1): the structural latent is reproducible across starts. Posterior means pile
up at the grid edges (max r 0.52-0.56, GTM edge effect).

Bin reconstructions (13:45): `klem_bins.py` bins the on-grid particles along the first principal direction of the
klem posterior means into NB equal-count bins -> selection file for rnasep.simple (932,496 rows = meta order);
`run_klem_cells_rec.sbatch` (env KDIR, NB) runs reconstruct3D nstates=NB at full box and compares every bin map with
the 10 docked Degenhardt conformers (corr of bin-minus-mean vs conformer-minus-consensus at 60-20 A, best conformer
per bin, end-to-end difference vs conformer pair differences). First run on the gate latent: job 53888008 (NB 5).

GATE BINS vs DEPOSITED CONFORMERS (13:50, 5 bins of 9.2k particles, maps 11-17 A): the real latent's main axis runs
from the {37, 18, 10} end to the {27, 28, 49} end of the Degenhardt family: bin 1 best = 37 (+0.29; anti 27/28/49
-0.23/-0.26/-0.22), bin 2 = 18, bin 3 = 48, bin 4 = 49 (+0.24; anti 10/18/37 -0.22), bin 5 = 28; end-to-end difference
correlates +0.17..+0.28 with every (27|28|49) - (10|37|18) conformer difference at 60-20 A. Same ordering as the toy's
truth axis = first direct link between a per-particle latent on real data and the paper's conformers (modest
correlations: 9k-particle maps). Figure rnasep_klem/klem_gate_npz/rec/fig_bins_vs_conformers.png. Full-run bins
(7, ~90k particles each) chained after job 53887975.

Full-scale replicate from a random start launched 15:05 (rnasep_klem_fullrand -> rnasep_klem/klem_full_rand, readout + 7-bin reconstructions chained). The npz full run was at it 4/9 at 15:00 (~200 s per iteration at 20 A, ~600 s at 14 A).

Bag-prior gate on real data (15:30): bags = cls_expansion sub-classes lab692 mapped to meta rows (rnasep_klem/bags_lab692.txt), beta 2, 50k, npz init -> rnasep_klem/klem_gate_bags (run_klem_rnasep.sbatch now passes BAGS/BAGBETA); readout chained.

Bag-prior gate verdict (16:05, lab692 bags, beta 2, 50k): NMI(cell, sub-class) 0.318 (circular: the bags are the
sub-classes) but NMI(cell, parent class) 0.130 (vs 0.031 without bags) and nuisance NMI tilt 0.023 / rotation 0.029
(vs 0.006 / 0.010): the bag prior injects the classes' VIEW component into the latent; soft-EM state NMI unchanged
(0.037), S1 separation weaker (0.55 vs 0.9 cells). Not adopted at beta 2; a within-parent-centred bag term would be
the clean version if 2D information is wanted.

FULL 692k RUN (53887975, npz init, 9 iterations, 2.9 h; 16:40): 89 % on the grid, mean max r 0.41; nuisance NMI tilt
0.009 / rotation 0.013 / defocus 0.005; NMI(cell, sub-class) 0.067 vs WITHIN-PARENT SHUFFLE NULL 0.036; within-parent
NMI median 0.065 over 112 parents vs null 0.026 (2.5x); NMI(cell, parent) 0.028 vs null 0.001 (a small view part);
soft-EM state NMI 0.056 vs null 0.000, S1 at z1 2.78 vs 1.9 for S2/S4; 15-NN vs the soft-EM labels 0.47 (chance 0.39).
TRAP: NMI has a large finite-sample bias: the 50k gates' "within-parent NMI 0.21" was mostly bias (450 particles per
parent over 49 cells x 16 sub-classes); klem_real_readout.py now prints shuffle nulls next to every NMI; quote the
excess over the null only.

RUN-TO-RUN REPRODUCIBILITY (klem_match.py, orientation-free: cell volumes matched by deviation correlation, per-particle
agreement through the matching): the two 50k gates (npz vs random init) agree at best-match deviation corr median 0.43,
soft per-particle cross-run cell-deviation corr 0.23 vs 0.005 permuted, 43 % within one cell, junk agreement 97 %; the
full run vs the 50k npz gate: 0.60 / 0.33 / 67 % / 96 %. A run's main axis maps to a coherent but CURVED path through
the other run's grid (grid parametrisation is not shared), so never compare z1 to z1 directly; use the matching.

EMBEDDING FIGURES (klem_umap.py, run locally with the manifoldem python): UMAP of the 51-D responsibility vectors
(Hellinger) = a 2-D sheet with z1 and z2 running across it, two thin strands = the junk flowing into the free
components, tilt/defocus uniform, S1 in one lobe; the flex_pca latent (z692) UMAP of the same particles is a featureless
ball with one pocket = the klem junk. Figures: toy55/figs/fig_rnasep_gate_umap.png, fig_rnasep_flex_vs_klem_umap.png.

CONCURRENT STREAMS launched 16:20-16:50 (user: "continue to work on these issues concurrently"): real 50k gates with the
within-parent-CENTRED bag prior (--bag_parents, 53888927, klem_gate_bagc) and the grid deconvolution (--deconv 1,
53888928, klem_gate_deconv), readouts via run_klem_readout.sbatch (readout + match vs klem_gate_npz); real 50k pose
alternation run_klem_alt_rnasep.sbatch (bins of the gate -> refine3D_states local 20->8 on rnasep.simple with only the
gate particles selected -> poses copied into meta_alt.npz -> klem --subset on the refined poses; 53888937,
rnasep_klem/alt_gate1); random-start full replicate 53888246 still running (readout + bins chained).

FULL-RUN BINS (7 x 88,100 particles, full box, gold-standard FSC 0.143 at 6.2-9.4 A; rnasep_klem/klem_full_npz/rec,
copied to ~/Downloads/rnasep_klem_20261015/full_bins/bin1-7.mrc): bin 1 best conformer 18 (+0.36; 37 +0.32, 10 +0.24;
anti 49 -0.31, 28 -0.23), bin 2 = 37 (+0.30), bins 3-4 = the middle (conformer 1, +0.17/+0.19), bins 5-6 = 49 (+0.25,
+0.29), bin 7 = 28 (+0.20); end-to-end difference +0.23..+0.33 with every (27|28|49) - (10|37|18) difference at 60-20 A
(gate: +0.17..+0.28). The main latent axis of the real data IS the Degenhardt family axis, with sharper maps than the
50k gate; the middle bins are near the consensus, the ends carry the signal.

GATE VERDICTS, real 50k npz init (17:20; plain gate for reference: nuisance 0.006/0.010/0.004, S1 0.9 cells along z1,
junk 8 % in free comps): (1) grid DECONVOLUTION alpha 1: nuisance unchanged, within-parent NMI 0.227 vs null 0.206
(same excess as plain), S1 separation 0.6 cells (weaker), only 3.5 % in the free comps (junk pulled onto the grid, junk
agreement with the plain gate 0.88), mean max r 0.75; toy worse -> NOT adopted. (2) CENTRED bag prior (within-parent
sharing): nuisance tilt 0.016 / rotation 0.017 (plain 0.006 / 0.010, uncentred 0.023 / 0.029), NMI(cell, parent) 0.067
vs null 0.016 (plain 0.028 vs 0.014): the view still leaks at half the uncentred level, because within-parent sub-class
differences are themselves partly view-private (bins80 finding); within-parent NMI 0.374 vs null 0.181 is circular
(bags = sub-classes); S1 separation unchanged; toy worse -> NOT adopted. Bag priors are dead in both forms; 2D
information can only enter as an initialisation. Remaining levers: particles per cell (full scale) and the discrete
split layer for composition (--split_rounds, testing on the toys).

RANDOM-START 692k REPLICATE (53888246, 17:50): nuisance NMI <= 0.008; sub-class NMI 0.041 vs null 0.029, within-parent
0.043 vs 0.025 (weaker than the npz run's 0.067 / 0.065), soft-EM 0.037, max r 0.36 (broader). Bins: the same family
axis (bin 7 = 27/28 at +0.53, bins 3-4 = 10/37 at +0.38 / +0.33; end-to-end +0.26..+0.41) but folded (bin 1 leans 27
too): the axis is found from a random start, its parametrisation is not.

THE CLOSED STATE IS THERE BUT NOT PER PARTICLE (closed_subclasses.py, 18:10): every cls_expansion sub-class average
scored against pose-mean renders of deposited closed 18 vs open 28 (klem convention conv 0 / shift -1 / scale 1, band
40-14 A; averages vs consensus render ncc 0.77 = frame check): 109 of 820 large sub-classes score > 0.5 (max 0.79),
holding 121k particles = 20 % of the data. Along the particle-level main axis the closedness of the sub-class averages
falls monotonically by quintile (0.22 -> 0.05) but the most closed sub-classes sit at t -0.4..-1.0, not at the extreme;
the SUB-CLASS MEAN latent positions separate the closed sub-classes at AUC 0.87 (direction (-0.81, -0.59); cells z1 1-2,
z2 2-4). Particle-level bins cannot purify them: closed-arm occupancy (arm/core density vs conformer 18 = 1.0) is 0.66
for the 14 % end bin, 0.71 for the 5 % end bin, consensus 0.50, open end 0.40; mask 200 vs 220 A changes nothing (the
deposited 18 has 1 % of its mass beyond 100 A). Pooling by sub-class is the lever: sub-class means have sd 0.72 vs
per-particle sd 1.75 (9x their s.e.m.; within a parent 5x). Launched: reconstructions of (a) 6 weighted k-means clusters
of the sub-class means, every member inheriting its sub-class's cluster (unsupervised; 53889414, rec_subpool_k6) and (b)
the closed-scoring sub-classes' particles vs the 109 most open vs the rest (reference, uses the deposited maps for the
selection; 53889415, rec_closedref), both at mask 220; arm occupancy readout chained.

Bins look soft because reconstruct3D output is unsharpened: relative B-factor vs the deposited map ~550 A^2; bin 1
agrees with deposited 18 to 10 A (FSC 0.5); the bin7 - bin1 difference reproduces between halves to 18 A (FSC 0.5) /
10 A (0.143) with difference power 30-40 % of the map power below 20 A. Sharpened copies (B -500, lp 8) in Downloads.

POOLED READOUTS (19:00): (a) sub-class-pooled bins along the main axis (every member inherits its sub-class's bin, 5 x
130k): bin 1 = 18 (+0.36) / 37 (+0.34), bin 5 = 49 (+0.17), end-to-end +0.24..+0.33 = the same as the particle-level
7 bins (no gain from pooling along PC1). (b) unsupervised k-means (6) on the 2-D sub-class means: weak (end-to-end
+0.15, occupancy <= 0.55; clusters cut across the axis). (c) SUPERVISED closed-scoring (109 sub-classes, 130k) vs the
109 most open (37k): closed map = 10 (+0.56) / 37 (+0.52) / 18 (+0.41), open map = 28 (+0.58) / 27 (+0.52), end-to-end
+0.45..+0.50 = 1.5x any latent bin, closed-arm occupancy only 0.60 (the closed family is itself a spread of conformers
at consensus poses). Conclusion: the 2D averages scored against renders carry ~1.5x the conformer signal of the
per-particle latent; the gap is per-particle noise, not information. Next: a sub-class-level E-step inside klem
(--avgs: each sub-class's even/odd average scored against the pose-mean renders of the cell DEVIATIONS, cross-half,
members inherit the responsibilities; the residual-vs-deviation scoring is what kept the closedness scoring view-free).

HARD SUB-CLASS E-STEP ON REAL DATA = VIEW PARTITION (klem_gate_avg, pool mode, 21:15): nuisance NMI tilt 0.066 /
rotation 0.060 (plain gate 0.006 / 0.010), NMI(cell, parent class) 0.39 vs null 0.014, soft-EM state NMI 0.023 (plain
0.042), bins' end correlations 0.17-0.23 (plain 0.29-0.36), match vs the plain gate 0.20 / 0.05. Even with the
residual-vs-deviation scoring, a whole sub-class forced into one cell makes the cells a partition by view: the average's
residual against the consensus render keeps view-specific structure (envelope, anisotropy) that some cell deviation
can mimic. Dead in the hard form. The unmix prior (soft proportions added to the particles' own likelihoods) is the
remaining form (klem_gate_unmix 53891332); judge it by the same nuisance NMIs first.

STATE GRAPHS (draw_state_graph.py, draw_deposited_graph.py; figs fig_rnasep_full_state_graph*.png,
fig_rnasep_state_graph_deposited.png, fig_rnasep_state_graph_joint.png; user: "could you do this on the states
downloaded from the paper", 21:40): (1) our 49 states laid out by MDS of their deviation-map distances (60-20 A) form a
horseshoe/ring: closed-scoring states (closedness +0.32..+0.37) at one end, the most open (+0.05..+0.14) at the other,
a side branch at high z1; nearest-neighbour similarity 0.96; the board's z1 wraps around it twisted. (2) The TEN
DEPOSITED conformers (deviation from their image-weighted mean, 60-12 A) form TWO clusters: open family {27, 28, 48, 49}
(27-28 +0.67, 48-49 +0.69) vs {1, 10, 18, 37, 46, 51} (10-37 +0.60, 18-10 +0.38, 51 loosely attached <= 0.28); between
the clusters strongly anti-correlated (10 vs 49 -0.71, 27 vs 46 -0.65). (3) JOINT embedding on the coarse grid: our
states form the RIM, the deposited conformers sit in the INTERIOR; each of our states correlates at most 0.35 (mean of
best matches) with any deposited conformer while the deposited nearest pairs correlate 0.65-0.69; best-match counts:
27:12, 37:9, 49:9, 51:7, 1:5, 28:3, 10:2, 18:1, 46:1, 48:0. Reading: the lattice model places states at the extremes of
the 2-D landscape (edge effect), the paper's conformers are interior points; pairs like 27/28 (0.67 similar) are the
55-state problem in miniature and far below our 0.35 match level. Supports the graph/density-placed states and
splitting along the path.

UNMIX PRIOR ON REAL DATA (klem_gate_unmix, 22:00): nuisance NMI tilt 0.019 / rotation 0.018 (plain 0.006 / 0.010, hard
pool 0.066), NMI(cell, parent) 0.084 vs null 0.015 (plain 0.028): the view still leaks, at a third of the hard form;
soft-EM NMI 0.038 (plain 0.042); bins' end correlations 0.15-0.23 (plain 0.29-0.36) = the conformer signal is DEGRADED;
match vs the plain gate 0.53 / 0.26. Fourth form of 2D coupling in the E-step to fail by view leakage. Mechanism: a
sub-class average's residual against the consensus render contains view-specific structure (shared by all sub-classes
of the same parent) that some cell deviation can mimic. Last principled form: PARENT-CENTRED unmix (unmixc): unmix the
difference between a sub-class average and its parent's average (same view, the view part cancels) and use only the
within-parent preference among cells as the prior.

PARENT-CENTRED UNMIX ON REAL DATA (klem_gate_unmixc, 00:05): the first 2D-coupled E-step that does NOT leak the view:
nuisance NMI 0.009 / 0.012 / 0.004 (plain 0.006 / 0.010 / 0.004), parent NMI 0.041 vs null 0.015 (plain 0.028), soft-EM
0.041 (plain 0.042); run match vs the plain gate 0.72 / 0.41 / 67 % within one cell = the most consistent variant (the
prior makes members of a sub-class land consistently); identification 9/10 (48 missing) at 0.33-0.56; BUT the PC1 bins'
end correlations 0.14-0.27 vs plain 0.29-0.36: no gain in conformer signal. Verdict: clean and neutral; kept as an
option (--avg_mode unmixc with --bag_parents), not default. 2D information in the E-step is closed as a line of gain.

FOCUS MASK ON THE REAL GATE (klem_gate_focus2, corrected script, 00:35): mask 3.5 % of voxels; nuisance NMI 0.008 /
0.013 / 0.006 (clean), soft-EM NMI 0.051 (plain 0.042), purity 0.505 (0.45), all 49 states reproduce across halves at
median 0.95; identification shifts rather than improves: 48 gains a state (0.40) but 18 loses its only one (0.22), 9/10
either way; the PC1 bins' end correlations drop to 0.15-0.24 (plain 0.29-0.36); match vs the plain gate 0.58 / 0.33.
Verdict: mixed/neutral on real data, as at toy full scale. Not a default; kept for child runs on small subsets.

REAL CLOSED-STATE REFINEMENT, level 1 (run_closed_refine.sbatch, 3 bins of the full latent -> refine3D_states local
20->8 A, 6 its, 1 h; 02:30): refined state 1 (closed-enriched end, 262k particles) = {37 +0.52, 10 +0.48, 18 +0.46; anti
49 -0.49, 48 -0.40} vs its seed bin (0.36 / 0.32 / 0.24): per-state refinement raises the conformer correlations by
~half on real data; state 2 (211k) = 49 / 48 (+0.34 / +0.31); state 3 (143k) = 27 / 28 (+0.26). Closed-arm occupancy of
state 1 stays 0.66: 18 (6.6 % of the paper's images) is a minority within the {37, 10, 18} family, a 3-state split
cannot isolate it -> level 2 (run_closed_level2.sbatch): klem on state 1's particles at the refined poses (no free
comps, pose marg + devw + pnoise, 20->8 A) -> identification -> k-means 3 -> refine3D_states 3 states -> occupancy.

**How to apply:** the readouts on real data are sanity checks, not truth: occupancy spread, mean max r, agreement of
the latent cells with the S1/S2/S4 states and the sub-class labels; final judgement needs reconstructions of the
cells along the latent and their comparison with the Degenhardt conformers. Related: [[klem-kernel-latent-em-2026-10-15]],
[[rnasep-24-refine3d-auto-quarter-unrefined]], [[rnasep-stacks-are-class-prealigned]].

## 2026-10-15 evening: real tree status
- Real expansion (expand1, parent klem_full_npz): family 1 = 410k particles, 25 states, identifies 7/10 deposited
  conformers (37 +0.79, 27 +0.76, 28 +0.71, 10 +0.66, 1 +0.55, 51 +0.48, 18 +0.46); family 2 = 197k, identifies 1/10
  (49 +0.34): the deposited landscape lives in family 1.
- Family-2 pipeline (job 53897698, 1h19): refined states match nothing deposited (best +0.12 to +0.21, closed-arm
  occupancy 0.31-0.54), feedback board 2/10 (37 +0.33, 49 +0.40). Family 1's pipeline (53897255) is the real test;
  the real tree realignment (run_tree_realign_rnasep.sbatch, waiter) follows both. Full 692k best run (klem_full_best,
  devw, 20->8): parent alone identifies 8/10 (mean 0.33); 7 bins at 6-9 A.
- Family-1 pipeline (53897255): refine3D_states 10 states at 7.8-9.8 A; refined maps match the deposited conformers
  only weakly (best +0.10 to +0.35; the BOARDS are the readout, as on the toy where refined hard labels are 0.25-0.37),
  but states 4 (10 +0.35, 37 +0.28, 46 +0.26, 18 +0.24) and 9 (18 +0.20) have closed-arm occupancy 0.73 = the purest
  closed states on real data so far (closed family level 2: 0.68). Feedback board at the refined poses started 19:18;
  a waiter submits the real level 3 (run_tree_level3_rnasep.sbatch: spectral bisection of the feedback board, 2
  grandchildren at the refined poses, focus, PMARG, 12->8, BAL 0.5) when it finishes. The real tree realignment now
  runs its parent and children with BAL 0.5 by default (PBAL/BAL).
- HEADLINE (family-1 pipeline, 2h50, DONE 20:36): the feedback child board at the refined poses (focus, PMARG, devw,
  pnoise, 14->8 A) identifies 10/10 deposited conformers, mean best-match corr 0.46 (1 +0.54, 10 +0.70, 18 +0.50,
  27 +0.79, 28 +0.69, 37 +0.76, 46 +0.41, 48 +0.55, 49 +0.31, 51 +0.52); child at consensus poses 7/10 (0.44), full
  692k best parent 8/10 (0.33). Board: expand1/fam1_pipeline/klem_fb, poses meta_fb.npz. Real level 3 on it
  (BAL 0.5) = job 53901575.
- Real jobs running 20:45: level 3 on fam1's feedback board (53901575, BAL 0.5), tree realignment (53901578: merged
  refined poses of both families -> parent BAL 0.5 devw PMARG 20->8 -> families -> children BAL 0.5), matched
  reconstruction (53901580, run_klem_matched_rec.sbatch: every cell of klem_fb labelled by its best deposited
  conformer, particles inherit, reconstruct3D 10 states at the refined poses, maps scored vs their conformer;
  deposited_match.py --save writes per-cell matches). Outputs under expand1/fam1_pipeline/klem_fb/matched_rec/maps.
- Matched reconstruction (53901580, 4 min): ten full-box maps at 4.8-7.8 A (46: 4.4k particles, 7.8 A; 37: 112k);
  conformer 18's map closed-arm occupancy 0.76 (purest closed state on real data). Scored as (map - toy consensus) the
  own-conformer corr is weak (0.01-0.29) and 46 "wins" several maps = the common real-vs-toy offset biasing every
  deviation; deposited_match subtracts the mean of the cells instead. score_maps.py (maps - mean of maps vs conformer
  - consensus / - mean of conformers, 60-20 and 30-12 A) re-scores both map sets (job 53901591). Always subtract the
  set mean when scoring full-box real maps against the docked conformers.
- CORRECTED SCORES (score_maps.py, set-mean subtracted): the ten matched-state maps are each best matched by their
  OWN deposited conformer at 60-20 A (10/10, mean own corr +0.54: 27 +0.82, 37 +0.79, 10 +0.61, 48 +0.54, 51 +0.54,
  1 +0.53, 28 +0.50, 18 +0.44, 46 +0.37, 49 +0.22) and 9/10 at 30-12 A (18 second to 37; mean +0.41). The ten
  refined states (unlabelled) pick 9 distinct conformers (48 missing, 37 twice) at 0.16-0.80. => The real RNase P
  landscape of the paper's 10 deposited conformers is recovered as ten full-box maps at 4.8-7.8 A from family 1 of
  the klem tree (maps: expand1/fam1_pipeline/klem_fb/matched_rec/maps/recvol_state01..10.mrc, order 1,10,18,27,28,
  37,46,48,49,51).
- CLOSED ZOOM, 116 A mask (53901718, 20 min): a board on the 66k closed-family particles (cells matched to 18/46/10)
  at the refined poses, focus + PMARG + BAL 0.5, 12->8 A: cells at arm occupancy 0.64-0.68 almost uniformly (d'
  0.34, 15-NN at chance): the closed family does NOT split by degree of closure at this noise; full-box maps of its
  purest cells 0.78-0.79 (3.6k / 2.5k particles) = the confidence-ranked 0.81. The classification is no longer what
  caps the arm. Tests: toy true-18 particles reconstructed at true vs realigned vs consensus poses (one-cell board,
  job 53902003 -> pose blur?); real closed-state local refinement against the closed map, mskdiam 240, 12->6 A
  (closed_refine_rnasep.sbatch, 53901992, output matched_rec/closed_refine). Occupancy-layer zoom 53901719 pending.
- CLOSED ARM EXPLAINED (toy job 53902003): PURE true-18 particles (10k, simulated from EMD-74209) reconstructed at
  TRUE poses measure 0.85 closed-arm occupancy by our metric (one-cell board at 8 A, 232 mask), 0.81 at the realigned
  poses, 0.80 at consensus poses. Real closed subsets 0.78-0.81 => the real closed state is at the level of a pure
  simulated closed state at our pose accuracy; the metric's ceiling at 8 A is ~0.85 (band + shrinkage + smoothing).
  Remaining levers: pose refinement against the closed map (job 53901992) for the 0.81 -> 0.85 gap; milder
  sharpening for the flexible arm (B -150) for viewing; the 100 A mask clipping is fixed in the downloads.
- REAL LEVEL 3 on family 1 (53901575, 1h53, BAL 0.5, 12->8 A): sub-family 1 (173k) 7/10 at mean 0.39 (10 +0.66,
  49 +0.46, 46 +0.36), sub-family 2 (237k) 7/10 at mean 0.55 (27 +0.81, 28 +0.70, 1 +0.61, 48 +0.60, 10 +0.56,
  37 +0.56); union 9/10: conformer 18 keeps corr 0.50 but owns no cell at level 3 (absorbed by 10/49 in sub-family
  1). Depth sharpens the open side; the closed conformer is best read at level 2 (klem_fb: 10/10).
- MASK CORRECTION: klem.py --mskdiam only masks the calibration consensus; the references are unmasked. The wide-mask
  toy tree reproduced the 200 run to the digit (0.61 / 0.64). Clipping was only in the sharpened downloads and the
  scoring masks (now 116 A) and reconstruct3D (220 -> 240).
- Closed-state local refinement against the closed map (53901992, mskdiam 240, 12->6 A): arm occupancy 0.76 -> 0.60,
  FSC 7.6 A: pose refinement against a partial-arm reference makes the arm WORSE (body dominates). Occupancy-layer
  closed board (53901719): layer cells with arm "occupancy" 1.2-1.4 but anticorrelated with 18 (blob artefact), tiny
  states, reconstruct failed: drop the occupancy layer. Running: focal P12/P10.1 boards (53902191 closed subset,
  53902192 family 1), defocus-split closed reconstruction (53902244: tip vs base density per defocus third tests the
  box-delocalisation explanation).
- TIP CEILING (toy job 53902293): pure true-18 particles at TRUE poses: arm 0.85, TIP (r > 95 A) 0.64, base 0.88;
  at the REALIGNED poses: tip 0.46, base 0.86; consensus poses: tip 0.46. Real closed cells: tip 0.58-0.61 (coarse)
  and 0.60-0.67 full box => the real closed particles are at/above the pure-particle ceiling at our pose precision.
  The P12-P10.1 contact ("helix cut off") is lost to pose precision at the periphery (3.8 deg at 114 A = 7.6 A), not
  to classification, mask, box or composition. What recovers it: per-state refinement of the closed state at full
  resolution with core-aligned poses and focal refinement of the periphery (the paper: 66k particles, 2.9 A). Box
  delocalisation refuted (defocus split). First closed-state refine3D attempt (12->6 A, local, vs the partial-arm
  map) measured worse on its own masked output; its poses re-measured by plain reconstruction = job 53902316.
- CLOSED-STATE LOCAL REFINEMENT VERDICT (53901992 + plain rec 53902316): refine3D single state, local, 12->6 A,
  mskdiam 240, started from our partial-arm closed map: poses got WORSE (plain reconstruction tip 0.65 -> 0.26, FSC
  5.3 -> 7.6 A). Discard. Per-state refinement of the closed state must follow the paper's order: poses aligned on
  the rigid core (core mask), periphery handled by focal classification/refinement, with the purified docked set.
  Real closed particles at the family poses: arm 0.76, tip 0.65, base 0.78, FSC 5.3 A (the reference numbers).
- TIP CEILING at 6 A (toy 53902343): true poses tip 0.65 / base 0.89; realigned poses tip 0.45 / base 0.86: band
  independent. The metric caps the tip at ~0.65 for a perfect reconstruction of our type (Wiener shrinkage, noise,
  smoothing vs the deposited 2.9 A sharpened map); the real closed particles at the family poses (tip 0.65) sit AT
  that cap => the closed state is recovered as well as this pipeline can render it; "cut off" vs the deposited map is
  a resolution/sharpening gap at the periphery (5.3 A vs 2.9 A), to be closed by per-state refinement in the paper's
  order (core-aligned poses, focal purification/refinement) and local sharpening, not by more classification.
- FOCAL P12/P10.1 BOARD on the closed subset (53902191, 20 min; klem.py --focus_mask = closed-minus-open mask of the
  deposited maps, 83k voxels after dilation): cells spread along docked/undocked (focal score -0.57..+0.73; the
  unfocused closed board was uniform); top-3 docked cells (3.2k particles) corr +0.40 with closed-18 inside the mask
  vs +0.15 for the undocked cells. => focal classification separates the P12-P10.1 contact state. Follow-up: rank all
  66k by responsibility-weighted focal score, reconstruct top 15k / next 10k / bottom 15k (focal_rank_rec_rnasep.sbatch);
  the family-1 focal board (410k, 53902192) pending.
- FOCAL PURIFICATION RESULT (53902471): closed-family particles ranked by responsibility-weighted focal score: top
  15k arm 0.82 / TIP 0.75 / base 0.83 at 6.5 A; next 10k tip 0.70; bottom 15k (undocked) arm 0.46 / tip 0.19. The
  40k "18-matched" set (tip 0.65) was a docked/undocked mixture; the focal board resolves the P12-P10.1 contact.
  Maps: closed_focal_closedfam/rank/maps (state01 top 15k), sharpened copies conf18_docked_top15k_sharpB150.mrc /
  conf18_undocked_bottom15k_sharpB150.mrc delivered to Downloads/rnasep_klem_20261015/recovered_conformers. Next: the
  family-1 focal board (410k, 53902192) for the full docked population (~13 % of apo expected), then per-state
  refinement of the docked set in the paper's order (core-aligned poses, focal periphery).
- SIMPLE (flex2D branch / snap 5e1e6aeee) has NO mask-file option for refine3D (only pcg_mskfile for the PCG solve):
  focal refinement of the periphery is not available; "core alignment" = refine3D with a smaller spherical mskdiam
  (170 A, excludes the arm tip at 95-115 A). run_docked_refine_rnasep.sbatch (53902548): the 15k docked particles,
  core refine3D local 10->8 A, unmasked reconstruction scored vs the family poses (tip 0.75 before).
- USER: the docked top-15k map is "not closed" by eye. Contact VISIBILITY at 8 A (core median = 1, contour 0.5):
  deposited 18 tip p99 2.36 / 81 % of tip voxels above contour; deposited 27 (open) 0.89 / 18 %; docked top 15k
  1.31 / 66 %; all 18-matched 40k 1.15 / 54 %; undocked 0.29 / 0 %. The docked map's contact is present at a low
  contour but its peak is half the deposited map's: a resolution/peakedness gap (6.5 vs 2.9 A), not identification.
  The integrated "occupancy" metric overstated closure; use the visibility metric (contact_visibility.sbatch).
  Route: full docked population from the family-1 focal board -> proper per-state refinement (ML-regularised, to the
  FSC limit, not local-to-6-A) -> postprocessed map. Pose-marginalisation A/B on the toy running (53902669).
- TRAP (2026-10-15 23:20): `pose_policy=local` is a refine3D_STATES parameter only; plain refine3D ignores it and
  runs its default refine=shc GLOBAL search. Both "local" closed-state refinements (240 A / 12->6 and core 170 A /
  10->8) were therefore global re-searches from scratch at lp 10 with 15-40k particles and degraded the poses
  (tip 0.75 -> 0.52, FSC 6.5 -> 7.8 A). A local refine3D needs refine=neigh (or shc_neigh) with athres (deg).
  Per-state refinement is PARKED by the user ("get out all the states first"); when resumed, use refine=neigh.
- FAMILY 2 COMPOSITION CHECK (fam2_composition.sbatch, 23:40): mean cell maps vs the consensus on the coarse grid:
  family 1 keeps 0.97 of the consensus's ordered voxels, family 2 0.92; no large missing region in either (largest
  missing 2503 / 1064 A^3, a domain would be > 10^4 A^3); family 2 has a small extra blob (2065 A^3) at the box centre
  (radius 9 A) and the family 1 - family 2 difference's largest region (1126 A^3) is also central. => family 2 is NOT
  a missing-domain compositional state; it reads as a lower-contrast / more flexible or differently-conformed
  population (the paper's other 45 apo conformers are candidates). The compositional question moves to toy10c.
- USER HYPOTHESIS (23:50): the 2D cls expansion shows the contact because of kernel-weighted averaging or the DENOISED
  particles. Test running (docked_denoised_rnasep.sbatch): the 15k docked particles reconstructed by klem's own
  one-cell cross-half backprojection from the denoised stacks (6 A band) vs the raw reconstruct3D map, contact
  visibility at 8 A + montage (figs_now/docked_denoised_vs_raw_lp8.png). Note: klem's E-step reads the denoised
  stacks; all reconstruct3D maps so far used the project's stacks (raw, to be confirmed by the stk field listing).
  The sub-class closedness scores (closed_subclasses.py) were never saved to disk (no *_closedness.txt found).
- DENOISED HYPOTHESIS EXCLUDED (53902683): both the 2D expansion project and the refinement project reference the
  DENOISED stacks (22_denoise_project/denoise_raw_particles_part*.mrcs), so every 3D map so far already came from
  denoised particles. klem's own one-cell reconstruction (6 A, Wiener shrinkage toward the consensus) UNDERSTATES the
  arm: docked 15k tip p99 0.99 / 32 % above contour vs reconstruct3D 1.31 / 66 % on the same particles -> never judge
  the arm on klem cell maps; use reconstruct3D. Remaining 2D hypothesis: kernel-weighted view-conditioned averaging
  (sub-class purity); the 2D route job (closed_2d_route_rnasep.sbatch, closedness per sub-class average vs
  deposited 18/27 renders, members of the top decile reconstructed) tests it.
- REAL STATE COUNT with identification (53902737; fam1 klem_fb + l3 sub1/sub2 + fam2 klem_fb, 196 surviving cells):
  similarity 0.6 -> 30 states (median 22k particles), 8 deposited conformers identified (1,10,27,28,37,48,49,51),
  8 unmatched; 0.7 -> 39 states, 9 identified (+46), 14 unmatched; 0.8 -> 56, 9 identified, 21 unmatched. 18 never
  dominates a cluster (its 2 cells merge into neighbours). With the toy's 2-3x over-split factor: ~15-20 distinct
  states on 15 % of the data, 9/10 deposited among them, plus ~8-14 unmatched reproducible states (candidates for the
  paper's other 45 apo conformers). Next: reconstruct the unmatched states at full box (deliverable maps).
- 2D sub-class closedness (closed_subclasses.py vs deposited 18/27 renders, 40-14 A): scores span -0.66..+0.80; top:
  sub 37/34/38 (parent 10), 1026/1029 (parent 263), 1537 (parent 449, n 3302), 945 (parent 242), 116 (parent 19),
  1617, 552. The 2D route job (closed_2d_route_rnasep.sbatch, members of the top decile reconstructed) and the
  unmatched-states job (run_unmatched_states_rnasep.sbatch: leaves clustered at 0.6, clusters with best deposited
  corr < 0.3 and >= 5k particles reconstructed per family, maps in rnasep_klem/unmatched_states/fam{1,2}/maps) are
  running. count_states.py --save/--at writes per-cell cluster labels.
- 00:30 resubmissions WITHOUT pose marginalisation (toy A/B: worth 0.02 at the leaves, 7x cost): real tree
  realignment 53902816 (PPMARG=0 CPMARG=0, RTAG=_realign1; the marginalised 53901578 was cancelled after 3 h in its
  parent), real level 3 on family 2 53902817 (PMARG=0; 53902726 cancelled; fam2_l3 dir wiped). run_expand_rnasep /
  run_tree_level3_rnasep now take PMARG; run_tree_realign_rnasep takes CPMARG for the children.
- 2D ROUTE RESULT (23:47 cluster time): sub-class closedness median 0.09, 90th pct 0.44, max 0.80; the top decile =
  179 sub-classes = 182k particles (26 % of the data, far above the ~13 % closed fraction): reconstructed at 4.7 A,
  tip visibility 0.73 / 24 % above contour, weaker than the focal docked 15k (1.31 / 66 %). Sub-classes are mixed
  inside; a decile cut is not a closed selection. klem one-cell maps again understate the arm (0.35): dropped from
  the comparisons. Retry with a strict cut (closedness > 0.6, ~30 sub-classes).
- UNMATCHED STATES (run 53902780): 26 clusters at 0.6 over fam1 l3 sub1/sub2 + fam2 klem_fb; 7 unmatched clusters
  with >= 5k particles, ALL in family 2 (98.7k particles; e.g. 27k best 37 at +0.25, 18.5k best 10 at +0.13):
  the states beyond the deposited ten live in family 2. Maps: rnasep_klem/unmatched_states/fam2/maps.
- UNMATCHED STATE MAPS (53902780): 7 family-2 clusters (5-27k particles) reconstructed at the fam2 refined poses: FSC
  8.4-10.6 A; deviations vs the ten deposited conformers all weak (best +0.08..+0.18, set mean subtracted): genuinely
  outside the deposited ten (other apo conformers or more flexible / partially unfolded particles). Maps in
  rnasep_klem/unmatched_states/fam2/maps/recvol_state01..07.mrc (order = cluster 18, 22, 13, 12, 11, 20, 14).
- 2D ROUTE VERDICT (strict cut, 53902828): 44 sub-classes with closedness > 0.6 = 57k particles, 5.2 A; tip
  visibility 0.73 / 22 % above contour = the decile cut (0.73 / 24 %), far below the 3D focal docked 15k (1.31 /
  66 %). Closedness of a sub-class average (vs deposited 18/27 renders at 40-14 A) does not select docked particles:
  the averages are view-conditioned mixtures and the score tracks a global motion, not the P12-P10.1 contact. DEAD;
  the 3D focal board is the tool for a specific contact. (Denoised klem one-cell maps again understate: 0.39.)
- USER (2026-10-16 00:05): disagrees that the closed contact is a resolution matter; "take the states suspected to be
  closed and refine them". Job 53902xxx closed_states_refine_rnasep.sbatch: the closed family (66k) as 3 states from
  the focal ranking (docked 15k / next 10k / rest), refine3D_STATES pose_policy=local (the mode that worked for the
  family poses; plain refine3D ignores pose_policy) 10->6 A maxits 8 mskdiam 232, then unmasked reconstructions and
  contact visibility per state; the refined docked map goes to Downloads as conf18_docked_refined_state_sharpB150.mrc.
- REAL LEVEL 3 on family 2 (53902817, 18 min WITHOUT pose marginalisation vs hours with): sub-families 102k / 95k,
  49/49 reproducible cells each (0.98), but deposited matches 1/10 (49 +0.31) and 0/10 (means 0.17 / 0.11): family 2
  is entirely outside the deposited ten at every level. Its leaves feed the distinct-state count (rerun queued).
- REAL LANDSCAPE COUNT, all six leaf boards (53903122; fam1 klem_fb + l3 sub1/2, fam2 klem_fb + l3 sub1/2; 294
  surviving cells): similarity 0.6 -> 41 states (8 deposited identified, 18 unmatched), 0.7 -> 61 (9 identified,
  35 unmatched), 0.8 -> 91, 0.9 -> 143. With the toy's 2-3x over-split: ~15-25 distinct states on 15 % of the data,
  9 of the deposited ten plus ~10-15 unmatched (family 2). Paper: 55 apo conformers.
- SUSPECTED-CLOSED STATES REFINED AS STATES (53903047, refine3D_states local, 3 states, 10->6 A, 8 its, 22 min):
  74 % of particles changed state; "docked" state shrank 15.0k -> 11.9k and LOST the contact (tip 0.40 / 0 % above
  contour, arm 0.56, 7.3 A); the two absorbing states (29k, 25k; 5.4-5.7 A) carry it at 0.98-1.07 / 47-54 %, below
  the plain focal docked 15k (1.31 / 66 %). The body dominates the 6 A likelihood and the hard EM reorganises around
  it. NEGATIVE. The refined "docked" map copied to Downloads (conf18_docked_refined_state_sharpB150.mrc) is the WORSE
  one. Last fair variant: the same with a 170 A core mask (CORE_MSKDIAM) so the periphery cannot steer the poses.
- FOCAL BOARD ON ALL OF FAMILY 1 (53902192, 409,609 particles, 2 h): a docked corner of 8 cells with closed-minus-open
  +0.80..+1.11 (closed-family board max was +0.73); top 3 cells 26.5k particles. Ranking + lever-arm test running as
  rnasep_focal_lever (53903807: top 30k / next 20k / bottom 30k reconstructed, then refine3D refine=neigh athres=8 on
  the top 30k as ONE state, full 232 A mask, 10->6 A). Box/sub-class diagnostic 53903808 (arm mass vs radius, pose-blur
  of the deposited map at 2/4/6 deg, sub-class ANOVA of the focal score + montage). Toy pipeline pose errors: angle
  median 3.98 deg (p90 7.35), shift 1.44 px -> 8 A tip displacement at 114 A radius (lever-arm blur candidate).
- The 2026-10-14 focused flex runs (pcg_mskfile region1/2/3, deconvolved mode labels, rnasep_coupled/focus) carry only
  weak docked enrichment (best component 12 % docked-top-30k vs 7 % board average; regions cover 49/13/39 % of the
  P12/P10.1 focal mask): a whole-region focus is not the contact-site focus; the klem focal board separates far better.
- CLOSED DEFICIT DIAGNOSED (2026-10-16 ~01:00, job 53903808): the box is NOT it (90 % of the deposited arm mass inside
  90 A, nothing beyond 114 A; a 232 A soft mask costs tip 0.96, 240 A nothing; a 200 A mask would give 0.17). POSE BLUR
  IS: the deposited closed map blurred by the toy pipeline's pose errors (angle median 3.9 deg, shift 1.4 px) gives
  tip 0.81 / contact p99 1.41 / p90 1.11; our family-1-wide docked top-30k map (score >= 0.84, 5.4 A) measures 0.78 /
  1.46 / 1.10. At 2 deg median the blur would give 0.94 / 1.90. Docked population: 39.7k particles above +0.8, 74.6k
  above +0.6 of the 409k board (paper: ~13 % closed). Top 30k > next 20k (tip 0.78 vs 0.69): mild mixture at the margin.
- REFINE-AS-STATES IS DEAD BOTH WAYS: refine3D_states 3 states local, full mask (53903047) and core mask 170 A (53903562)
  both reshuffle 74-75 % of the particles and the docked state loses the contact (0.40 / 0.52); the absorbing states
  hold it at 0.88-1.07. Never re-classify the docked set with the body likelihood; refine poses only (refine3D
  refine=neigh, nstates=1) = the lever-arm test rnasep_focal_lever 53903807.
- 2D SUB-CLASS PRIOR: 16 % of the per-particle focal-score variance lies between the 1756 expansion sub-classes (null
  0.4 %; particle-vs-subclass-mean corr 0.40), but ranking by own score + shrunk sub-class mean gives a slightly WORSE
  docked map (contact p99 1.36 vs 1.46, tip>0.5 0.62 vs 0.63 at 5.2 vs 5.4 A). The 2D averages see docking, the
  per-particle focal score already contains it; no sub-class route for the closed state.
- LEVER-ARM REFINEMENT NEGATIVE (53903807, 2026-10-16 01:05): refine3D refine=neigh athres=8, one state (docked top
  30k), full 232 A mask, 10->6 A: body FSC 5.4 -> 5.3 A but tip 0.78 -> 0.49, contact p99 1.46 -> 1.10, >0.5 0.63 ->
  0.32. Every pose or state refinement of the docked set has blurred the arm (3-state full mask, 3-state core mask,
  1-state lever arm); the arm is sharpest at the family poses. So the arm is NOT a rigid function of the body pose
  within the docked set: within-docked heterogeneity (paper: 7 closed conformers + continuous). Policy at the docked
  node: do not align to the region, classify on it (focal board within the docked 30k = 53903923, ZTAG _docked30k).
- 2D "CLOSED" SUB-CLASSES ARE ~30 % CLOSED: the 36 sub-classes with shrunk mean focal score >= 0.40 (8652 board
  members) have 29 % of their members in the docked top 30k; the other 71 % (B, own-score median +0.39) reconstruct
  OPEN (contact p99 0.67, >0.5 0.16, like deposited 27 at 0.89/0.18); C (in the top 30k) 1.11/0.48 at 9.4 A. A 2D
  average shows a coherent 30 % docked arm over the smeared 70 % undocked arms, which is why they look closed and why
  the 2D route fails. Montage with ids: toy55/figs/closed_subclasses_montage.png.
- CORRECTION: the refine=neigh "lever-arm" run (53903807) was NOT local: pose change angle median 2.4 deg but p90 77
  deg (max 180), shifts median 5.7 px (p90 16.6): refine=neigh re-searches in-plane/shifts freely (athres only bounds
  the projection direction). Its arm verdict (0.78 -> 0.49) is VOID. refine3D_states needs nstates > 1. The genuinely
  local refiner in this snapshot is refine3D refine=pose_cont (Levenberg-Marquardt, route shift_then_joint|joint):
  rerun as lever2_pose_cont_rnasep.sbatch (53903926). klem's own --pose_every step (_pose_chunk, 3x3x3 rotation grid x
  2 refinements, Gaussian prior, no shifts) is strictly local; --pose_mask/--pose_ref added for the node policy test.
- REAL TREE REALIGNMENT (53902816, no pose marginalisation, 1.5 h): tree regrown at the family-refined poses
  (606k particles moved, |d e2| median 2.0 deg; meta rnasep_klem/realign_realign1/meta_realigned.npz): parent 8/10
  (0.42); family 1 (334k) 9/10 (0.46, 18 at +0.51); family 2 (248k) now identifies 27/28/37/49/51 (5/10, 0.28; before
  the realignment family 2 matched nothing). Family overlap new x old Jaccard 0.73 / 0.56. Children in
  rnasep_klem/expand_realign1; second pass command printed at the end of rnasep_tree_realign.53902816.out.
- WITHIN-DOCKED FOCAL BOARD (53903923, docked top 30k, ZTAG _docked30k, 10 min): the docked set still spans a docking
  gradient (cell focal score +0.52 .. -0.50); the top-3 cells (2.6k particles) have 114 % of the deposited closed map's
  density in the P12/P10.1 region (family-wide docked cells: 102 %, bottom cells 100 %). Ranking + reconstruction of
  its top/next/bottom 10k = focal_rank_rec 53903940 (FD closed_focal_docked30k).
- CLOSED STATE MILESTONE (53903940, 2026-10-16 01:25): ranking the docked 30k by the WITHIN-docked focal board and
  reconstructing its top 10k at the family poses gives arm 0.97 / tip 0.84 / base 0.99, contact p99 1.73, p90 1.31
  (deposited 2.36 / 1.64), FSC 7.1 A; next 10k 1.43/1.07, bottom 10k 1.16/0.91; the family-wide top 30k was 1.46/1.10.
  So the tree policy "classify on the region that reveals itself at the node, keep the poses" is what brings the arm
  out; every pose/state refinement of the docked set blurred it. 10k closed particles at 7 A is the expected yield of
  15 % of the paper's data (66k at 2.9 A from all of it). Map: Downloads/.../conf18_docked_within10k_sharpB150.mrc.
  Next level (docked10k board + rank) running.
- SUB-CLASS 1413 (the user: "sub1413 is what I need in the 3D"): 274 members, focal median +0.64 (36 % above the
  top-30k cut, 26 % in the top 30k, 9 % in the within-10k). Pose-mean renders at its member poses (render_subclass.py):
  our docked maps match its average at ncc 0.97 / arm region 0.85-0.87, the deposited closed map 0.96 / 0.77, open
  0.91 / 0.47, consensus 0.94 / 0.61. The 3D docked maps CONTAIN what the average shows; the difference the user sees
  is display (B -150 at 7 A suppresses the lower-local-resolution arm): delivered lp10 unsharpened + lp8 B-50 versions.
- KLEM POSE-POLICY TEST on the docked 30k (53903931, one-cell maps, identical shrinkage; arm/tip occupancy): family
  poses 0.76/0.75; lever (cell ref, everything; moved 2.3 deg median) 0.77/0.76; core-masked (mask too small, 0.8 % of
  coarse voxels; moved 6.3 deg) 0.71/0.66; REGION-masked (P12 focal mask; moved 4.5 deg) 0.80/0.80. Aligning on the
  arm region sharpens it slightly; aligning on everything is neutral; a too-small core blurs it. Combine with the
  within-node classification (tip 0.78 -> 0.84 unshrunk) = the node policy: classify on the region, then align on it.
- LEVEL 3 INSIDE THE DOCKED 10k (53903966 board + 53904007 rank): the docking gradient continues (top 3k arm 1.05 /
  tip 0.88, next 3k 0.99 / 0.84, bottom 3k 0.87 / 0.83) but the contact visibility no longer improves (top 3k p99 1.55
  / p90 1.23 at 9.0 A vs the 10k's 1.73 / 1.31 at 7.1 A): below ~10k particles the metric is noise-limited. The
  within-10k set is the best purity/particle-count trade-off on this 15 % subset; the tree has resolved the closed
  node as far as the particle count allows.
- refine3D refine=cont (53903939) on the docked 30k: also not a local move (angle median 5.1 deg, p90 19.4; shifts
  median 6.7 px; 30 % moved > 10 deg) although the body FSC stays 5.4 A; tip 0.78 -> 0.55, contact 1.46 -> 1.07. Both
  SIMPLE refine modes re-estimate shifts/in-plane broadly in this snapshot and blur the arm. STOP using refine3D for
  docked-node pose polishing; klem's rotation-only pose step (region-masked) is the only local move that helped
  (region_align_rec job tests it with a proper reconstruction).
- CLOSED STATE CONFIRMED BY THE USER (2026-10-16 ~01:40): conf18_docked_within10k_lp10_unsharpened.mrc and
  conf18_docked_within10k_lp8_B50.mrc = "yup those latest two maps are it". Recipe: family-1-wide focal board (409k)
  -> rank by the responsibility-weighted focal score -> top 30k -> focal board INSIDE the 30k -> rank -> top 10k
  reconstructed at the FAMILY poses (no refinement) -> display without heavy sharpening. See
  [[deliver-maps-unsharpened-for-judging]].
- OVERNIGHT 2026-10-16 (user: "keep pushing this method, work on this overnight", "you would need to derive the correct
  mask"): TRUTH-FREE node recipe launched on family 1. Scripts in toy55: derive_regions.py (spatial components of the
  board's reproducible variability: at thr 0.5 ONE component holds 99 % of the mass, covers 84 % of the P12 reference
  mask), derive_regions2.py (k-means on per-voxel motion profiles over the cells -> 5 regions; region 1 = 1024 coarse
  voxels, 22 % of the mass, covers 43 % of the P12 reference with 59 % of it inside), node_states.py (states = average-
  linkage clusters of surviving cells at sim 0.7, particles ranked by responsibility mass, top min(30k, 0.8 mass)),
  rank_by_target.py (level 2: cells scored against the parent state's deviation), region_board_rnasep.sbatch (board ->
  states -> reconstruct3D -> deposited_match/vis -> level-2 boards with --focus 1 for the 3 largest states),
  run_knobfree_rnasep.sbatch (driver 53904039: regions 1-2, jobs 53904040-41) and run_knobfree2_rnasep.sbatch (driver
  53904042: regions p1-p5, jobs 53904043-47). Outputs under rnasep_klem/knobfree/region_{1,2,p1..p5}/. Deposited maps
  used only in the identification readouts. Toy10c PCA tree: 61 0.72/0.80/0.86, 62 0.35/0.36/0.36 (the added-density
  state is never recovered at any level, random or PCA init).
- TIMINGS of the current method (16 cores unless noted): parent board on 692k 2 h 13 min; family split + children
  1 h 54 min; family-1 pose refinement + feedback board 2 h 50 min (family 2: 1 h 19 min); tree level 3 1 h 53 min;
  tree realignment 1 h 30 min; one focal board on 409k 2 h (on 30k: 10 min, on 10k: 4 min); rank + 3 reconstructions
  2 min (8 cores + 6 workers). Whole tree ~10 h, plus ~2 h per region board; the E-step is Python/Fortran-kernel hybrid.
- REGION-ALIGNED POSES ON THE WITHIN-10k (53904016): klem's rotation-only pose step masked to the P12 region raises
  the one-cell arm occupancy to 0.92 / tip 0.89 (core-masked 0.82 / 0.74; family poses ~0.76-0.80 on the 30k) with
  |d e2| median 1.9 deg. BUT the reconstruct3D maps from the RE-IMPORTED poses (oritab in stack order ->
  import_particles stktab -> reconstruct3D) are wrong for both variants (contact 0.30, arm 0.42, FSC 8.1 A: a
  consistent molecule in the wrong frame/hand) -> the import path, not the poses. import_check_rnasep.sbatch
  (53904048) re-imports the UNMOVED poses (and a shift-sign flip) to find the convention; until it passes, do not
  trust any map reconstructed from a text table through import_particles.
- FAMILY 2 HOLDS LITTLE OF THE CLOSED STATE (53903998 + 53904148): focal board on 197k family-2 particles: cell
  scores max +0.48 (family 1: +1.11); top 15k by score: arm 0.60 / tip 0.30 at 9.0 A (family-1 top 30k: 0.85 /
  0.78). Do not merge family-2 "docked" particles into the closed set; the closed population is in family 1.
- NODE POLICY, ALIGNMENT HALF CONFIRMED (53904147, 2026-10-16 02:50): with the import fixed (selection after import;
  the unmoved control reproduces the family-pose map exactly: 1.73 / 1.31 / 0.60 / arm 0.97 / 7.1 A), klem's
  rotation-only pose step masked to the P12 REGION on the within-10k set gives contact p99 1.89, p90 1.40, arm 1.02
  (|d e2| median 1.9 deg); core-aligned poses 1.37 / 1.08 / 0.89. Recipe: classify on the region at the node (poses
  fixed), then align the selected particles on the region; never on the core, never on everything.
  Maps: klem_fb/region_align_10k/region/maps/recvol_state01.mrc (+ for_download versions in Downloads).
- WIDE-TOY POLICY TEST WITH TRUTH (53903938, toy10w: closed arm intact, consensus-refined poses 3.99 deg median):
  klem rotation-only pose step on the TRUE conformer-18 particles: whole-cell reference ("lever") 3.99 -> 2.89 deg
  (p90 7.41 -> 5.75), arm 0.80 -> 0.83; core-masked 3.48 deg, arm 0.77; REGION-masked 4.59 deg (worse poses) but arm
  0.87 (best map). Same pattern for 51 (lever 2.88) and the open 27 (lever 3.36). POLICY: refine poses against the
  whole cell reference (best body+arm poses, propagate down the tree); the region-masked move is a final map-making
  step for the node's region (like a local refinement), its poses are not to be propagated.
- REGION ALIGNMENT ON THE TOP 30k (53904193): contact p99 1.46 -> 1.52 (p90 1.10 -> 1.15), arm 0.85 -> 0.89, FSC 5.4
  -> 5.9 A (body loses as in the toy); core-aligned 1.18 / 0.76. Smaller gain than on the purer 10k (1.73 -> 1.89):
  the region move pays in proportion to the purity of the set. Best closed maps: within-10k region-aligned (1.89, 7.3 A)
  and the 30k region-aligned (1.52, 5.9 A). Maps under klem_fb/region_align_30k/{region,core2,unmoved}/maps.
- POLICY SETTLED ON REAL DATA (53904234, 2026-10-16 03:25): within-10k docked set, klem rotation-only pose step against
  the WHOLE cell reference (lever): contact p99 1.85 / p90 1.39, arm 0.99, FSC 7.1 -> 6.9 A (body AND arm improve,
  |d e2| median 1.2 deg); REGION-masked 1.89 / 1.40 / 1.02 but FSC 7.3 (arm at the body's expense); core-masked
  1.37 / 0.89. Same ordering as the wide toy with truth. RECIPE: at a node, classify on the region with poses fixed;
  then refine the node's particles' poses against their own cell reference (whole reference, local, rotation-only);
  propagate those poses; the region-masked move is an optional final map polish only.
- OVERNIGHT REGION BOARDS: all seven finished their level-1 states (e.g. region 1: 12 cell clusters -> 8 states of
  24-30k; region p1: 7 clusters -> 5 states of 30k) and ran level 2; their identification/contact readouts crashed
  because region_board_rnasep.sbatch copied maps with `seq -w` (unpadded below 10 -> recvol_state1.mrc never exists):
  fixed (`seq -f %02g`) and a fix-up job (region_readout_fix_rnasep.sbatch, 53904365) redoes the copies, the
  deposited_match identifications and the contact tables for every region/level after the drivers finish.
- OVERNIGHT TRUTH-FREE RESULT, LEVEL 1 (driver 53904042, 5 motion-pattern regions, 04:06): 245 surviving cells (all
  reproducible, self-reproducibility median 1.00), cell clusters -> 5-8 states of 24-30k per region; landscape count
  over the 5 boards: 31 / 36 / 53 states at similarity 0.6 / 0.7 / 0.8 (median 59-73k particles per state). Spatial
  driver 53904039 (regions 1-2): same stages. Identifications vs the deposited maps come from the fix-up job 53904365.
- TEN CONFORMERS AT WHOLE-REFERENCE-REFINED POSES (53904243, 48 min): the family-1 matched states (4k-112k particles
  each) through the klem pose step (|d e2| median 1.27 deg), re-imported (ctf=flip + selection) and reconstructed as
  10 states: identification unchanged (10/10, mean best-match corr 0.53 -> 0.53), FSCs slightly better (28: 6.9 ->
  6.5 A, 37: 5.1 -> 4.9, 46: 9.0 -> 8.1, 10: 6.0 -> 5.9), contact numbers of the broad 18-matched state unchanged
  (1.15 -> 1.17). The pose step pays on PURE sets (docked 10k: 1.73 -> 1.85, 7.1 -> 6.9 A), little on broad mixtures.
  Maps: klem_fb/lever_conformers/maps (recvol_state01-10 = conformers 1 10 18 27 28 37 46 48 49 51),
  for_download/conf*_leveraligned_lp10_unsharpened.mrc (on batch2 only).
- Spatial-region driver 53904039 (regions 1-2, 04:19): 98 surviving cells, 14 / 16 / 20 states at similarity 0.6 /
  0.7 / 0.8 (median ~35k particles per state). Together with the 5 motion regions: 343 reproducible cells over 7
  truth-free boards; identifications from the readout fix-up 53904365.
- TRUTH-FREE REGION BOARDS, IDENTIFICATIONS (fix-up 53904365, 06:10): region_1 (the data-derived moving region, 99 %
  of the variability, no deposited mask) -> 8 level-1 states of 24-30k that identify 8 / 10 deposited conformers
  (mean best-match 0.43; the hand-masked family board gave 10/10 at 0.46); its state 2 (30k) is a docked state with
  contact p99 1.41 / arm 0.79 / 5.4 A (deposited-mask top-30k: 1.46 / 0.85 / 5.4), state 4 1.19 / 0.75; level 2 inside
  state 2 (parent-direction ranking, 10k): 1.23 / 0.75 (the deposited-mask within-10k reached 1.73: the level-2 ranking
  target must be the node's OWN moving region, not the parent cluster direction). Motion-pattern regions: p1 4/10
  (0.32; docked-like states 1.16-1.28 / 0.76-0.81), p2 2/10, p3 2/10 (0.36), p4 2/10. The mask CAN be derived: the
  node's own variability recovers most of the landscape and the closed state at near hand-mask quality.
- 7-BOARD TRUTH-FREE LANDSCAPE COUNT (06:12): 343 reproducible cells over the 7 region boards -> 41 / 50 / 73 states at
  similarity 0.6 / 0.7 / 0.8 (median 64k / 51k / 31k particles per state); the hand-driven tree gave 41 / 61 / 91.
  Region p5 level-2 s1 (10k): contact 1.21 / arm 0.77.
- LEVEL 2 ALONG THE NODE'S OWN MAIN AXIS (53904408): region_1 state 2's within board, pc1 of its cell deviations
  (singular values 1 : 0.83 : 0.52), the minus end 10k: contact 1.47 / p90 1.12 / arm 0.80 (parent-direction top 10k
  1.23 / 0.92 / 0.75; deposited-mask within-10k 1.73 / 1.31 / 0.97). Better than the parent direction, still short of
  the hand mask because the within-state main axis mixes motions (second singular value 0.83 of the first). The
  truth-free way to close the gap = recurse: derive the within-state board's OWN regions (derive_regions2) and run a
  focal board on the arm region there, then rank. region_p1 state 3: plus end 1.17, minus end 1.06 (no gain).
- RECURSION INSIDE THE TRUTH-FREE DOCKED STATE (53904432): the regions derived from the within-state board are
  central (centroids 8-17 A, 44/21/10 % of the P12 reference covered): inside the docked state the variability is the
  body's, not the arm's; a focal board on the derived region splits the 30k into 6 small clusters (1.4-7.7k) with
  contact 1.15-1.44 (best 4.9k: 1.44 / 1.11 / arm 0.80) = no better than the level-1 state (1.41) or the own-axis end
  (1.47), far from the hand mask (1.73). Lesson: the region must be INHERITED down the tree (the parent's region that
  defined the node), re-derived only when a new one reveals itself; test = within-state board with the parent's
  region_1 mask (REG r1s2p).
- INHERITED-MASK WITHIN-STATE BOARD (53904439, region_1 state 2's 30k with the parent's region_1 mask): 12 clusters
  of 0.6-5.2k; best state (4.1k) contact 1.44 / p90 1.11 / arm 0.87 (the best truth-free arm so far; hand-mask
  within-10k 1.73 / 1.31 / 0.97). The whole-periphery region splits the state on several motions at once, so the
  docked-pure subset is small and diluted. Next (running): segment the inherited region by motion pattern INSIDE the
  state (derive_regions2.py --within parent mask) so the arm's motion group separates from the rest of the periphery.
- MOTION-SEGMENTED INHERITED REGION INSIDE THE DOCKED STATE (53904447/8): the main sub-region (534 coarse voxels, 60 %
  inside the P12 reference, covering 25 % of it) gives a 6-cluster within board; best state (6.4k) contact 1.42 /
  1.09 / arm 0.78. Truth-free recursion plateaus at contact 1.44-1.47, arm 0.80-0.87 (hand mask 1.73 / 0.97): the
  data-derived regions find the docked state at the node but do not isolate the docking site sharply enough for the
  final purification; a finer segmentation (reproducibility-weighted variance, more motion groups) is future work.
