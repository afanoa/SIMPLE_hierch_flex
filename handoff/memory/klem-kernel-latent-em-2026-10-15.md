---
name: klem-kernel-latent-em-2026-10-15
description: "klem.py = kernel-latent EM (generative topographic mapping style) prototype the user chose over deformation fields (2026-10-15): 2-D latent grid of free cell volumes smoothed by a grid kernel, soft cross-half E/M steps on a coarse Fourier grid, off-grid junk components; synthetic 2-D family (motion + composition + junk) recovered from a random start; toy10/toy10c runs"
metadata:
  type: project
---

User direction (2026-10-15 ~09:00): "what about a nonlinear model close to what we have like kernel PCA ... I am
not sure I want to pursue deformation fields because they wont do well with compositional heterogeneity" -> agreed;
the kernel goes on the LATENT side as the decoder (GTM): V(z) = kernel-weighted combination of free cell volumes on a
2-D grid; flex_pca is the special case with a linear map. Composition is free (a cell may lack a domain).

Implementation `scripts/cls_expansion_proto/toy55/klem.py` (mirrored on batch2 routes/toy55): n1 x n2 grid + nfree
off-grid components; E-step = band-limited image likelihood vs CROSS-HALF cell volumes (particles of half h scored
against volumes from the other half), uniform mixing; M-step = soft Wiener backprojection of every particle into every
component on a coarse Fourier grid (nearest voxel, Friedel mates, bincount) with weights = responsibilities convolved
with the grid kernel (sigma annealed), shrinkage toward the consensus (eps = 0.02 x median coverage); schedule
20/14/10 A x 3 iterations; per-shell reference scale + noise from the consensus residual; init = two contrast-free
PCs of a flex latent (PCA init, GTM standard) or random; forked multiprocessing (fork context, module globals);
readouts vs truth: posterior-mean (z1, z2) Spearman with the axis, 15-NN, d', closed AUC, cell purity, junk capture.
`klem_synth.py` = synthetic family (blob moving +-6 px with z1, blob present/absent with z2, 10 % junk, random
poses/CTF, box 48): from a RANDOM init, 6 iterations: |Spearman(z1, motion)| 0.80, composition AUC 0.975, 15-NN 0.76
(chance 0.11), junk 84 % in the free component with 0 % conformer particles. `run_klem.sbatch` (16 cores, 96 GB):
jobs 53885926 (toy10, lp 12 init -> lat/klem), 53885927 (random init -> lat/klem_rand), 53885928 (toy10c -> latc/klem).

Results (toy10, 50k subset, PCA init from lp 12, 7x7 + 2 free, 20/14/10 A x 3, 6 min total): one axis orders the
conformers in the truth sequence; Spearman 0.77 at it 6 (14 A) -> 0.73 at it 9 (10 A: with ~1k particles per cell
the 10 A stage is too thin and the posteriors sharpen on noise), 15-NN 0.30 (chance 0.16), purity 0.33, closed AUC
0.74; junk NOT captured by the free components from the PCA start (0.00; the random start captured 96 %), the second
grid axis unused. Random init (50k): fitted Spearman 0.46, purity 0.26 at it 7. Fixes queued 10:50: --init spectral
(kNN Laplacian eigenvectors of the lp 12 latent) and junk-like initial free components (35 A blob, half-contrast).

Fast trio final (10:55): PCA init toy10 0.73 / 15-NN 0.30 / purity 0.33; random init 0.41 / 0.23 / 0.26 (junk 91 %
in free comps); compositional toy 0.70 / 0.24 / 0.27 with the deleted-domain state 61 sitting ON its parent 49 and the
added-blob state 62 ON its parent 46: composition NOT separated, the second grid axis unused (grid effectively 1-D).
Queued: 1-D grid of 15 cells (fast1d) to test finer resolution along the main axis.

Junk-like initial free components work (fastj, 11:05): junk 99 % in the free comps, 1 % of conformer particles;
Spearman 0.78 (it 6) -> 0.75, 15-NN 0.32, purity 0.34, closed AUC 0.80 = best configuration so far (7x7, PCA init
from lp 12, junk init). 1-D 15-cell grid: 0.55-0.59 / 0.27 / 0.33 (finer 1-D scale orders worse than the smoothed
7x7). Spectral init starts slower (fitted 0.31 at it 3 on toy10, 0.67 at it 4 on toy10c). Still effectively 1-D:
queued fastanneal (kernel sigma 2.5 -> 0.7, SOM-style topology formation) on both toys.

fastj final (11:10): Spearman 0.73 (peak 0.78 at it 6), 15-NN 0.32, purity 0.34, closed AUC 0.80, junk 99 % captured.
Latent figure: the vertical axis orders the motion in the truth sequence (27/28 -> 48/49 -> 1 -> 51 -> 18/46 -> 37/10);
the closed conformer 51 separates along the SECOND axis (z1 1.5 vs 3.3-4.0 for the rest; 18 at 2.6): the second
dimension is used for the closed rearrangement, consistent with the model-class result (closed = largest change).
Per-particle posteriors are broad (max r 0.58) and the posterior means pile up at the grid edges.

Composition not separated by the smoothed grid (toy10c fast: 61 sits on 49, 62 on 46): the grid kernel ties
neighbouring cells, so an isolated 8 % compositional state gets averaged into the motion axis. Extension 11:20:
`--nkseed K` = extra FREE (unsmoothed) components seeded from k-means states on the lp 12 latent (half of each
particle's initial mass on its k-means component, half on the grid) = continuum (grid) + discrete candidates (free
comps); per-free-component composition readout added. Runs: toy10c_klem_fastk 53887239, toy10_klem_fastk 53887240.

Spectral (Laplacian) init on toy10c (11:25): 0.68 peak -> 0.62 final, 15-NN 0.22, purity 0.26: no gain over the PCA
init (0.70 / 0.24 / 0.27); the fixed-graph coordinates of the noisy per-particle latent do not seed a better grid.

FULL-SCALE toy10c (200k, PCA init, no junk init; 53885928, ~45 min): Spearman 0.77, 15-NN 0.29 (chance 0.13), purity
0.31, closed AUC 0.79 at it 9: more particles per cell lift every metric over the 50k run (0.70 / 0.24 / 0.27) and
the 10 A stage no longer hurts.

Full-scale random init (200k, 53885927): fitted Spearman 0.80, 15-NN 0.28, purity 0.31, junk 83 % captured, but the
conformer means stay within 0.5 cell of the grid centre (broad posteriors, the axis runs obliquely): the random start
organises at full scale as well as the PCA start, only less legibly; at 50k it reached fitted 0.41. Full-scale toy10c: 61 at (1.93, 2.03) vs parent 49 (2.13, 1.71),
62 at (1.54, 4.63) vs parent 46 (1.57, 4.62): composition still unresolved by the smoothed grid at full scale.

FULL-SCALE toy10 (200k, PCA init from lp 12, no junk init; 53885926, 35 min): Spearman 0.80, 15-NN 0.37 (chance 0.16),
pairwise d' 1.87, purity 0.39, closed AUC 0.76; conformer order along z2 = truth order exactly (27/28 0.7, 48/49 1.8,
1 3.4, 18 4.5, 46 4.3, 37 5.0, 10 5.3), 51 separated on z1 (2.6 vs 3.9-4.7). BEST so far; oracle ceiling 0.84 /
0.59 accuracy. Spectral init on toy10 (50k): fitted 0.48, purity 0.33 (no gain). Full-scale runs with the junk-like
free-component init launched 11:40 (toy10_klem_junk, toy10c_klem_junk -> lat/klem_junk, latc/klem_junk).

Kernel annealing 2.5 -> 0.7 (11:45, 50k): toy10 0.78 peak / 0.75 final, 15-NN 0.32, purity 0.34, closed 0.79 = same
as the junk-init run (no gain); toy10c 0.73 / 0.24 / 0.27, 61 displaced 0.4 cell from its parent (still unresolved).
Annealing does not create the second axis by itself.

k-means-seeded extra free components (nkseed 4, 50k, 11:50): they crystallise conformer-enriched discrete groups
(toy10: free3 = conformer 27 at purity 0.46 vs base 0.15, free4 = 46 at 0.29, free5 = 28 at 0.29; toy10c: free4 = the
deleted-domain state 61 at 0.28 vs base 0.12) but drain the continuum: grid Spearman 0.66 / 0.62, 15-NN 0.25 / 0.19,
purity 0.28 / 0.22 (vs 0.78 / 0.32 / 0.34 without them). Unsmoothed free components compete with the grid for
ordinary particles. Queued: 10x10 grid with the kernel ending at 0.5 cells on toy10c (fast10x10) to see whether a
finer grid lets a compositional cell specialise.

10x10 grid, kernel to 0.5 cells, toy10c 50k (12:10): 0.74 / 15-NN 0.23 / purity 0.28; the deleted-domain state 61 at
(2.5, 3.4) vs parent 49 at (2.8, 3.0) = a 0.45-cell displacement, not a separation. A finer continuum grid does not
resolve composition either. Full-scale junk-init runs at it 7: toy10 0.81 / 0.36 / 0.40 / closed 0.83 (best yet),
toy10c 0.78 / 0.27 / 0.32 / 0.81 with 61 still on 49. Offered to the user (12:10): run the EM on real RNase P (200k,
consensus poses, ~1 h Python) and work composition on the toy with an explicit per-cell occupancy; awaiting the word.

Full-scale toy10c with junk init (53887282, 12:45): Spearman 0.78, 15-NN 0.30 (chance 0.12), purity 0.33, closed AUC
0.84, junk 99 % captured; the deleted-domain state stays 0.5 cell from its parent (unresolved). Real RNase P gates
(50k, random / npz init) submitted 12:40: jobs 53887799 / 53887800.

FULL-SCALE toy10 with junk init (53887281, 12:50): Spearman 0.80, 15-NN 0.39 (chance 0.15), d' 2.13, purity 0.42,
closed AUC 0.88, junk 99 % captured = BEST configuration (7x7 grid, PCA init from lp 12, junk-like free comps,
20/14/10 A x 3, 200k particles). Real-data sanity script: klem_real_readout.py (NMI of latent cells vs tilt /
rotation / defocus bins = nuisance check, vs sub-classes and parents, vs soft-EM states; figure).

Composition work (13:05, user: "still develop the approach on the toy in the meantime"): (1) --nlayer: parallel
unsmoothed grids (block-diagonal kernel); on the synthetic 5x1 grid + 2 layers the layers refined the motion axis
rather than the blob switch (closed AUC 0.77 vs 0.975 for the smoothed 5x2 grid): layers are not the composition
lever by themselves. (2) INSIGHT: the PCA init drops PCs correlated with the per-particle contrast (junk filter), but
a DELETED DOMAIN lowers the contrast, so the filter can drop the compositional axis (toy10c lp 12 PC1 had 61 at its
extreme) -> --init pca12 (raw PCs 1-2, junk left to the free comps). Runs: toy10c_klem_fastpca12, toy10c_klem_fastlayer.

Jobs 13:10: toy10c_klem_fastpca12 53887908 (raw PCs init), toy10c_klem_fastlayer 53887909 (2 layers), toy10_klem_long
53887929 (200k, 12 iterations, kernel to 0.5 = sharpness test), toy10_klem_fasttrue 53887930 (50k with the TRUE poses
= cost of the 4 deg pose error); real RNase P gates 53887799/53887800 with readout jobs after them.

toy10c raw-PC init (pca12, 13:40): 0.71 / 0.25 / 0.27, 61 still on 49: no gain; the contrast filter was not what hid composition.

TRUE-POSE CONTROL (toy10, 50k, 14:00): with the simulation's true poses instead of the refined ones (4 deg median):
Spearman 0.88, 15-NN 0.48, d' 2.53, purity 0.49 (refined poses, same setup: 0.78 / 0.32 / 0.34). The 4 deg pose
error is the single biggest loss: in-loop pose refinement against the latent-decoded references is the next lever.

IN-LOOP POSE REFINEMENT (14:30, --pose_every/--pose_step/--pose_prior): per-particle local search against the argmax
cell volume on a 3x3x3 rotation-vector grid. Findings on the synthetic: an ARGMAX over candidates follows the noise
(from TRUE poses the error jumps to 12.8 deg in one step; exact references + noise-free data do give the identity, so
the scoring is right), a posterior-MEAN update random-walks (4 -> 8 deg from the truth), so a Gaussian prior on the
update (pose_prior) is needed; also: Euler inverse m2euler_simple verified (R = Rz(e3) Ry'(e2) Rz(e1)). Two parallel
layers on toy10c: 50/50 split with no conformer preference (composition unresolved by layers either).

Pose prior results (synthetic, 14:40): prior 2 deg still drifts from the truth (2.3 -> 4.6 deg) and does not improve
perturbed poses (5.5 -> 6.2); prior 1 deg: 5.8 -> 5.5 deg (marginal). Per-particle pose evidence is too weak at this
SNR for an in-loop Python search; the practical route is alternation with SIMPLE's refine3D_states on klem cells
(overnight: local multi-state refinement improved poses 4.0 -> 3.5 deg). Launched toy10_klem_alt (run_klem_alt.sbatch):
7 bins of lat/klem_junk -> refine3D_states local 20->8 -> klem 50k with the refined poses (compare to 0.78/0.32/0.34).

Tied-layer occupancy model (--tie_layers: layer 2 = base + one shared learnt component, base pooled over all particles;
14:55): on the synthetic blob switch it FAILS (closed AUC 0.50, layers split 50/50, component power 4-16 % of base =
noise) whereas the plain smoothed 5x2 grid had found the switch (AUC 0.975): the symmetric layer init does not break
toward composition. 7 composition variants tried, none resolves the toy10c deleted-domain state; parked. Volumes for
the user: ~/Downloads/rnasep_klem_20261015/ (gate_bins/bin1-5.mrc, deposited_conformers/, toy_consensus.mrc, README).

Tied layers on toy10c (15:10): degenerate collapse (all conformer particles into a free component, grid empty): the tied-component M-step is unstable as written; parked with the other composition variants.

BAG PRIOR (15:15, user: 2D information is fine if it helps): --bags FILE (ptcl2D class= per project row or one int per
row) + --bag_beta: each particle's log-likelihood gets + beta x the mean log-likelihood of its 2D sub-class (members
share evidence, keep their own term; a pure bag E-step would cap purity at the bag purity ~0.3-0.4). Test: toy10 50k
with the battery's pose-mean sub-classes (lat/exppm/ptcl2D.txt), beta 2 (toy10_klem_fastbags); compare to 0.78/0.32/0.34.

Bag prior on the synthetic (bags of 20 at 70 % purity, beta 2): 15-NN 0.98, purity 0.80, closed AUC 0.998 (per-particle: 0.76 / 0.65): sharing evidence within impure bags is a large gain when bags are conformer-enriched. toy10 run with the pose-mean sub-classes: job 53888327.

ALTERNATION WORKS (toy10, 15:35): 7 bins of the full-scale klem latent -> refine3D_states (local, 20->8 A, 6 its):
poses 3.98 -> 3.68 deg (p90 7.35 -> 7.00); klem 50k on the refined poses: Spearman 0.80 (it 6) / 0.77, 15-NN 0.39,
purity 0.40, closed AUC 0.87 vs 0.78 / 0.32 / 0.34 / 0.80 before (true-pose ceiling 0.88 / 0.48 / 0.49). One round
recovers ~half the pose-error loss; more rounds next (klem subset EM + final E-step over all particles -> bins ->
refine3D_states from the previous refined project).

--final_all added (EM on the --nsub subset, then one E-step over ALL particles -> latent for everyone; synthetic check ok). Multi-round alternation script run_klem_alt.sbatch (NROUND, KDIR, PROJSRC, ALTDIR): job 53888341 = 2 more rounds from the round-1 refined poses (lat/alt2/round1,2; note the first of these starts from the 50k-only round-1 latent, so its refine3D_states covers 50k particles; round 2 is full).

12-iteration full-scale toy10 run (kernel to 0.5; 15:45): 15-NN 0.41, purity 0.43, closed AUC 0.885, Spearman 0.77 (9 iterations: 0.39 / 0.42 / 0.88 / 0.80): more iterations saturate; the remaining levers are poses (alternation) and particles per cell.

Bag prior on toy10 (pose-mean sub-classes, 1542 bags / 50k, beta 2; 15:50): Spearman 0.65, 15-NN 0.26, purity 0.30 = WORSE than per-particle (0.78 / 0.32 / 0.34): the toy's sub-classes (majority purity 0.27-0.40) are too impure to share evidence; the synthetic gain needed 70 %-pure bags. Real-data verdict pending (gate_bags).

TRAP (15:55): a refine3D_states output project copied elsewhere prints 1 row (relative paths): each alternation round now rebuilds a fresh project from the previous round's refined pose table (POSETAB -> new_project + import_particles + update_project). Resubmitted 2 rounds (lat/alt2).

ALTERNATION PLATEAUS (toy10, rounds 2-3 from rebuilt projects, 16:30): poses 3.68 -> 3.65 -> 3.74 deg (round 2 slightly
WORSE), latent 0.38-0.39 / purity 0.39-0.40 / closed 0.87-0.88 in both rounds = the no-alternation full-scale run. Local
refine3D_states against 7 bin maps (mixtures) refines nothing beyond round 1 (cf. the direction-bin finding: mixture
references do not improve poses, state-pure ones do). Poses parked unless the real 50k alternation shows movement.

New options (16:20): --bag_parents FILE (bag prior centred on the parent class: only the within-view part of the
sub-class evidence is shared; the plain bag prior leaked the view on real data), --deconv ALPHA (one-step grid-kernel
deconvolution of the cell maps after every M-step, V + alpha (V - K V); sharpness lever), --subset FILE (restrict the
EM to listed particles, e.g. a gate's klem_latent.txt). run_klem.sbatch takes AXIS from the environment when there is
no battery log. Toy 50k tests: toy10_klem_fastdeconv 53888931, toy10_klem_fastbagc 53888932 (parents = the `cluster`
field of lat/exppm/ptcl2D.txt, 100 parents x 16 sub-classes; lat/exppm/parents.txt).

COMPOSITION, large-deletion regime (16:45): make_big_variant.py cuts 35 % of the most populous conformer's mass beyond a
plane normal to its principal axis (conf63, 15 % of the particles; the assembly-intermediate case); run_toy10d_sim.sbatch
(wrapper around run_toy10c_sim.sbatch with simd/projd/TAG=d, calibration copied) then run_toy10d_klem.sbatch (truth
axis over 11 states, klem random init 50k, composition readout for conformer 63: junk share, cells holding it, best cell
purity, overlap with the parent). Sim job 53888936; klem submitted by a Mac waiter when the sim ends.

The latent-to-volume map ("preimage") is the M-step itself: a Nadaraya-Watson kernel smoother over the grid whose
averaging operator is the Wiener backprojection, V_m = (sum_i w_im c_i F_i + eps V_cons) / (sum_i w_im c_i^2 + eps),
w_im = sum_m' K(m, m') r_im'. No parametric map; GTM minus the small ridge solve (hence the --deconv option).

Toy 50k verdicts (17:10): grid DECONVOLUTION alpha 1 (fastdeconv): fitted 0.66 / 15-NN 0.30 / purity 0.32 / closed 0.73,
max r 0.75 = WORSE than the baseline 0.78 / 0.32 / 0.34 / 0.80 (sharper references sharpen the posteriors on noise);
not adopted. CENTRED bag prior (fastbagc, within-parent sharing only): 0.67 / 0.27 / 0.29 over 38 used cells = worse
than per-particle, better than the uncentred bag (0.65 / 0.26 / 0.30): on toy10 the sub-classes are too impure to
share evidence in any form. Real-data gates of both variants: readouts 53888980 / 53888981.

SPLIT-LAYER TEST NEGATIVE (--split_rounds 4 on toy10 and toy10c, 50k, 18:20): after the schedule every cell gets a twin
layer with a random per-particle split; the cross-half reproducibility of the layer difference stays at 0.00-0.11 in
every cell (none above 0.5), layer-2 shares 0.47-0.53 for every conformer; on toy10c the extra freedom COLLAPSES the
grid (all conformer means at the centre, 15-NN below chance, junk capture lost). A random split finds no discrete factor
and halves the particles per reference. Composition by layers is dead in all four forms (parallel from the start, tied,
raw-PC init, post-schedule split). Remaining route for composition: hierarchical zoom (within a bin a missing domain
becomes the dominant factor) and the free components for large deletions (toy10d pending).

REAL 50k POSE ALTERNATION (alt_gate1, 18:25): refine3D_states local 20->8 A on the 5 gate bins moved the poses by a
median 3.8 deg (90th pct 9.4) and 0.7 px; klem on the refined poses: soft-EM state NMI 0.052 (plain gate 0.042), 15-NN
vs the soft-EM labels 0.469 (0.456), within-parent NMI 0.218 vs null 0.196 (no change), nuisance unchanged: a small
gain, consistent with the toy plateau. Full-scale alternation (7 bins, 692k) running: job 53889289.

HIERARCHICAL ZOOM TEST (run_klem_zoom.sbatch, toy10, 18:30): 3 bins of the full-scale latent (lat/klem_junk) -> a fresh
klem on each bin's particles alone (--subset, PCA init) -> does within-bin purity beat the global 0.42? This is the
candidate route to the paper's 55-conformer landscape: unsupervised top split, then recursion with per-bin poses.

ZOOM VERDICT (toy10, 3 bins of the full-scale latent -> fresh klem per bin, 20:20): within-bin purity 0.46 / 0.36 /
0.40 vs 0.42 global; every bin still holds all 10 conformers (15-NN 0.46/0.33/0.40 at chance 0.24-0.26). Recursion
on per-particle bins re-finds a weaker copy of the global structure: parents must be purer than particle bins.

SUB-CLASS E-STEP, hard pooling (toy10 50k, --avgs, class-frame convention calibrated at ncc 0.955, renders from all
members' poses): purity 0.20, 15-NN 0.19 (chance 0.15), max r 0.95 = WORSE than per-particle (0.34): the toy's
sub-classes are conformer MIXTURES (purity 0.27-0.40, uniform-view toy), so pooled hard assignment caps purity at the
sub-class purity and the cell maps blur. First version (renders from the subset's members only) was a pose-mismatch
bug, void. Real data (view-tight, state-enriched sub-classes) gate: klem_gate_avg (53890973, pool mode).

UNMIX PRIOR (20:40, --avg_mode unmix, default): each half-average's residual is unmixed by NNLS onto the cell
deviations (whitened by the even/odd noise), the proportions pi (floored at 1e-3) enter as log pi x avg_prior ADDED to
the members' own log-likelihoods: a sub-class may be a mixture, members keep their own evidence. Tests: toy10 50k
(toy10_klem_fastunmix 53891331), real gate (rnasep_klem_gate_unmix 53891332, readout + bins chained).

OTHER CAPACITY TESTS launched 20:10-20:45: 3-D grid (--n3; 5x5x3, toy10_klem_fast3d 53890979), full-scale TRUE poses
(coverage ceiling for poses, 53890618), full-scale 10x10 (53890619), ADAPTIVE kernel width (--adapt 0.5: per-cell
sigma from the cross-half reproducibility of the cell deviation; toy10_klem_fastadapt 53891334). state_coverage.py =
the per-state readout (toy10 global baseline: 10/10 cells owned, one-to-one purity mean 0.48, min 0.36; oracle
per-particle precision 0.38-0.73, mean ~0.6).

BOARD -> STATE GRAPH (21:10; user: "perhaps a board is not the best way to represent this?" after seeing the 7x7
montage fig_rnasep_full_grid_montage.png: rim crowded 17-34k, middle 2-5k, curved paths, run-to-run reparametrisation).
--kernel graph (--graph_after 3, --graph_nn 4): from iteration 3 the M-step kernel is computed from the states
themselves: similarity = cross-half correlation of the state deviations (noise cannot create a neighbour), d^2 = 2(1-S),
Gaussian bandwidth = median distance to the nn_eff nearest states with nn_eff annealed like the lattice sigma, column-
normalised; display coordinates = 2-D spectral embedding of the graph scaled to the board range (readouts, bins
unchanged); cells written as state_###.mrc. K is free of the 2-D geometry (49 and 100 tested: toy10_klem_fastgraph49
53891430, fastgraph100 53891431). Risks: collapse (one map) or fragmentation; guards: annealing + cross-half similarity.
Splitting a state = a new node beside its parent -> the tree / genetic expansion maps onto state_weight_set parents.

GRAPH KERNEL, first runs (22:20): the bandwidth rule (median distance to the nn_eff-th neighbour) gave 30-85 effective
neighbours per state (all states are 0.8-0.9 similar, distances compress) -> over-smoothed: 49 states collapse to 23
used, 100 to 21, purity 0.26-0.28 (lattice 0.34). Fixed by per-state bisection of the bandwidth to nn_eff effective
neighbours (1 / sum p^2); calibrated runs fastgraph49b / 100b (53891468 / 53891469) sit at 0.34 purity at iteration 6 =
lattice. TRUE POSES at full scale, final: 0.91 / 0.587 / 0.609 (refined poses 0.80 / 0.39 / 0.42): poses remain the only
large lever. User (22:30): minor pose refinement inside the loop is acceptable; the model's job is to identify every
state; per-state refinement is theirs; "I want to know stuff is working before pursuing 50 state refinement" -> the
49-state alternation was cancelled; the HANDOFF test runs instead: run_state_refine_toy.sbatch = unsupervised K=10
k-means on the toy posterior means -> refine3D_states local 20->8 -> state_purity.py per conformer before/after
(53892122). Identification readout for real runs: deposited_match.py (in run_klem_readout.sbatch): full run identifies
9/10 deposited conformers at corr >= 0.3 (48 missing, 18 and 46 with one state each), mean best-match 0.35.

FORTRAN KERNELS (23:10, user: "if you need to port things to fortran to speed up testing please do so"):
klem_kernels.f90 (OpenMP, bind(C)) = klem_estep (per-particle likelihoods vs the other half's references, trilinear
slices) and klem_mstep (nearest-voxel backprojection with Friedel mates, parallel over components, in-place
accumulation); built on the batch2 LOGIN node (allowed for builds) with gfortran 15: `gfortran -O3 -fopenmp -shared
-fPIC -o libklem.so klem_kernels.f90` in routes/toy55; klem.py loads it via ctypes (fort_lib / fort_estep /
fort_mstep; KLEM_FORT=path, KLEM_NOFORT=1 forces numpy; KLEM_CHECK=1 prints numpy-vs-Fortran differences on the first
chunk at iteration 0). TRAP: Fortran is case-insensitive -> dummy K/W vs locals k/w clash (renamed kk/wt). Layouts:
numpy C-order arrays passed as-is, Fortran sees reversed index order (refs(n3,n3,n3,K,2) <- (2,K,n3,n3,n3)).
Check jobs: toy10_klem_fortcheck 53892271 (5k, KLEM_CHECK) vs toy10_klem_npcheck (numpy twin) for timing.

PURITY LEVERS for the closed state (user 23:00: "we should make an effort to get purer closed states out of the
model"): --devw (shell weights 2F/(1+F) from the pooled cross-half FSC of the cell deviations, floored; the fine shells
where cell differences do not reproduce stop diluting the evidence) and --focus (E-step references = consensus +
deviation x soft mask of the reproducible variability of the cells, cross-half product in real space; no hand-drawn
region). Toy 50k tests: fastdevw 53892293, fastfocus 53892294, fastdevwfocus 53892295; judge by the purity of
conformers 18 and 51 (state_coverage). Real: closed_refine (run_closed_refine.sbatch, 53892252): 3 states from the full
latent -> refine3D_states local 20->8 -> closed-arm occupancy of the refined closed state (0.66 before).
Parent-centred unmix on the toy (fastunmixc): 0.71 / 0.31 / 0.33 = neutral (no longer hurts); real gate unmixc2 53892124.

COMPOSITION, assembly regime SOLVED (toy10d, 23:40): conformer 63 (conformer 49 minus 35.6 % of its mass beyond a plane,
15 % of the particles; sim 1 h 42, klem random init 50k 85 s with the Fortran kernels): 62 % of its particles in one
cell at purity 0.99, four cells above 0.5 purity, 8 % in the junk comps, none on its parent's cells. Large deletions
separate by themselves; only the 12 % single-domain deletion of toy10c remains unresolved. (Readout label bug: the
"overlap with the parent" line picks the most populous CSV row = 63 itself; ignore it.) The sim chain auto-submitted a
toy10d battery (cancelled, 53892634).

FOCUS MASK AT SCALE: toy full scale 0.82 / 0.367 / 0.402 / closed 0.83 vs baseline 0.80 / 0.39 / 0.42 / 0.88 -> the
+0.05 purity gain at 50k becomes -0.02 at 200k: it helps noisy maps, not good ones. Real gate focus: 15-NN 0.447,
purity 0.476 (readout + bins chained). Fortran timing at 50k: 77 s vs 336 s numpy (10 A iteration 11 s vs 66 s).

GENETIC EXPANSION v1 (user 23:45: "perhaps we can pursue the genetic expansion now"): run_expand_toy.sbatch +
expand_families.py: parent run -> spectral bisection of the state graph into families -> particles follow their argmax
state -> one child klem per family (--subset, band 14 -> 8 A, focus recomputed inside the family, Fortran) -> children
survive if the cross-half reproducibility of their deviation >= 0.5 (klem_states.txt, written by every run now) ->
merged labelling scored by coverage/purity vs the parent (0.42). Job toy10_expand 53892xxx (see log).

TRAP (00:20): in run_klem_rnasep.sbatch the options added after the AVGS block (--adapt, --devw, --focus, --kernel,
--graph_*) had been inserted INSIDE the `${AVGS:+ ...}` conditional, so every real-data run without AVGS silently ran
without them: the first "real gate focus" (klem_gate_focus) was a plain Fortran re-run of the npz gate (argmax agreement
0.996 with it = the Fortran path reproduces numpy on real data). Fixed; rerun = klem_gate_focus2 (53892787). Check with
`awk '/klem.py/{p=1} p{print; if(!/\\$/) exit}' script | fold` that new options sit outside every ${X:+...} block.
Full-scale 10x10 final: 0.81 / 0.365 / 0.410 vs 7x7 0.80 / 0.39 / 0.42 (neutral). Expansion v1 first attempt failed on
free-component argmax indices in the family split (fixed; rerun 53892790).

BATCH OF SIX (toy 50k, 00:40; baseline 0.78 / 0.32 / 0.34 / closed 0.80; conformer 51 purity 0.45): pose
marginalisation 3 deg (7 poses, logsumexp): 0.76 / 0.338 / 0.360 / 0.83, 51: 0.50 (+); tempered finish 0.5: 0.31 /
0.34 / 0.79 (neutral/worse: self-confirmation); per-particle noise: 0.326 / 0.353 / 0.82 (+ tiny); per-particle
contrast folded into C: 0.283 / 0.311, junk capture 0.52 (WORSE: scaling C breaks the junk components); mixing
proportions: 0.312 / 0.334, 31 used cells (states die; neutral/worse); combo with pscale: worse. Small gains that may
add: pose marg (+0.02), devw (+0.02), pnoise (+0.01), focus (+0.05 at 50k only) -> fastbest (all four, 50k),
fullbest (pmarg + devw + pnoise, 200k), fullpmarg (200k). Options: --pose_marg, --temp_final, --pnoise, --pscale,
--mixprop, --devw, --focus (all off by default).

EXPANSION v1 RESULT (toy10, 00:50): families 25 + 24 states (77k + 82k particles), between-family state similarity
-0.31 (the two anti-correlated families, as in the deposited set). Child runs (49 cells each, 14 -> 8 A, focus): merged
98 states over 122k particles; one-state-per-conformer purity mean 0.47 (parent 0.48) but the OPEN family's conformers
rise clearly (27 0.49 -> 0.63, 28 0.49 -> 0.60, 37 0.51 -> 0.62, 46 0.62 -> 0.72) while the closed 51 collapses 0.64
-> 0.13 and 38k real particles land in the children's FREE components ("conformer particles 0.20 / 0.26" in free comps):
an outlier conformer inside a family is absorbed by the junk components, and the focus mask (built from the family's
dominant variability) damps its distinguishing region. Also every child "survives" (cross-half 0.87 median): the
deviation-reproducibility test is too lenient for kernel-smoothed cells. v2 = children with NFREE=0 and no focus
(toy10_expand2, 14 -> 8 A; expand2b at 20 -> 10 A).

FULL-SCALE IDEAS (01:05): 8 A band with deviation-FSC weights (full8A): 0.80 / 0.403 / 0.437 / closed 0.88 vs baseline
0.80 / 0.39 / 0.42 / 0.88; 51 purity 0.67 (0.64), 18 0.40 (0.36): a small consistent gain; the 8 A shells get weight
0.06 (deviation FSC 0.03 there), the gain is the 10-12 A weighting. Junk-free second pass on the parent's on-grid
particles (fullnojunk, nfree 2 kept): WORSE (0.37 / 0.39 / closed 0.70) and the closed 51 collapses (purity 0.10; its
base rate in the subset only 0.03): CORRECTION - in the parent run every conformer incl. 51 sits 99 % on the grid
(the free comps hold the real junk); it is in JUNK-FREE runs (second pass, expansion children) that the free components,
with no junk left to model, turn into outlier sinks and absorb 51 (and 15 % of all conformer particles). Rule: free
components only where junk is present; nfree 0 for child runs and second passes. Test: fullnojunk0 (NFREE=0).

fullnojunk0 (junk-free second pass, NFREE=0, 01:15): 0.81 / 0.357 / 0.401 / closed 0.73; 51 purity 0.40 (recovered
from 0.10 with free comps, below the parent's 0.64), 18 0.46 (parent 0.36). A junk-free second pass gains nothing over
the parent; the free-component outlier-sink mechanism is confirmed (0.10 -> 0.40 by removing them).

HANDOFF TEST RESULT (toy10, 01:20; run_state_refine_toy.sbatch K=10): unsupervised k-means (10) on the posterior means
seeds refine3D_states (local 20->8 A, 8 its, hard labels). SEED one-state-per-conformer purity mean 0.32 (min 0.07,
6/10 conformers own a state) -> REFINED mean 0.43 (min 0.03, 8/10), recall 0.34 -> 0.44; per conformer 1 0.13->0.34,
10 0.26->0.44, 27 0.28->0.48, 28 0.43->0.53, 37 0.35->0.58, 46 0.44->0.59, 49 0.37->0.46, 51 0.49->0.50; LOST: 18
0.07->0.03, 48 0.40->0.31; poses 3.68 -> 3.56 deg. Per-state refinement PURIFIES klem seeds by a third and keeps 8/10:
the identify -> refine pipeline works on a known case; the weakest closed state (18, 6 %) is the casualty. Next: the
49 argmax cells as seeds (purity 0.48 > the k-means 0.32): toy10_state_refine49 (STATEMODE=cells).
EXPANSION v2b (20->10 A children, no free comps, no focus): 70 surviving states (cross-half >= 0.5) over 100k
particles, 9/10, purity mean 0.46 = parent (51 back at 0.63, 46 0.69, 37 0.56, 27 0.55; 18 0.25): removing the free
comps restores the closed state but the recursion itself adds nothing at the parent band. v2 (14->8 A) log checked
separately. fastbest (pmarg+devw+pnoise+focus, 50k): 0.81 / 0.327 / 0.379 / 0.83 = focus alone; gains do not add.

FULL-SCALE LEVERS (01:50): fullpmarg (pose marginalisation 3 deg alone): 0.81 / 0.405 / 0.434 / closed 0.89 (baseline
0.80 / 0.39 / 0.42 / 0.88), 51 0.67, 18 0.39; fullbest (pmarg + devw + pnoise): 0.82 / 0.393 / 0.437 / 0.90, 18 0.45,
51 0.66. Each lever +0.015-0.02 purity, they do not add beyond ~0.44 (true poses 0.61): the per-particle model is at
its limit with refined poses; the rest is pose accuracy = refinement's job.
49-STATE REFINEMENT FAILS IN SIMPLE (toy10_state_refine49, 01:45): refine3D_states nstates=49 -> every PROB_TAB_NEIGH
worker exits 1 after ~30 s with empty stdout/stderr (the probabilistic assignment table scales with nstates x nspace;
48 GB per worker); nstates=10 ran fine (the handoff). Practical limit ~10-20 states per refinement -> refine
HIERARCHICALLY: per family, <= 10-20 states each (toy10_state_refine_fam1/2 from the expand2b children).

EXPANSION v3 (14 -> 8 A children, no free comps, no focus, survival judged at <= 12 A; 02:05): 28 + 38 of 98 children
survive; merged purity mean 0.46 (parent 0.48), 8/10 own a state; 51 0.64, 46 0.63, 37 0.59, 27 0.53; 18 0.28, 1 0.30.
A fresh board per family does not beat the parent at either band; v1's open-family gains (27 0.63, 28 0.60) came with
the FOCUS mask (children are the noisy-map regime where focus helps) -> v4 = focus on, free comps off.
UMAP, toy full scale (figs/fig_toy10_umap_baseline_vs_best.png): the best configuration's responsibility manifold is
visibly articulated into lobes (51 its own lobe, 27/28 an arm, 10/37/46 a lobe, 48/49 a lobe) where the baseline is one
continuous sheet; in-plane 15-NN accuracy 0.40 -> 0.41 only. Real 692k run with the best configuration: 53893365.

LEVERS ROUND 2 (02:40, user: "implement all these ideas and submit jobs"): --pose_marg_n 27 (full rotation-vector grid),
--shift_marg px (five shifts), --pose_marg_mstep (joint pose-state posterior from the E-step's per-pose log-likelihoods
backprojects each particle at its neighbourhood poses), --temp_start T0 (deterministic annealing of the E-step
temperature over the first two thirds), --pnoise_m (per-particle noise factor as the M-step Wiener weight). Toy 50k:
fastpm27 53893661, fastpmshift 53893663, fastpmM 53893664, fasttanneal 53893665, fastpnoiseM 53893666, fastcombo2
(pmarg + pmM + pnoise + pnoise_m + devw) 53893667. Lever 5 = family-level pose feedback (run_family_feedback.sbatch:
each family's 10-state refinement poses -> child klem on the family subset, NFREE=0, best config) launched by a waiter
when toy10_state_refine_fam1/2 finish. Real data: full 692k best configuration 53893365 (readout + 7 bins + identification
chained), closed family level 2 53893778 (run_closed_level2.sbatch).

EXPANSION v4 = FIRST CLEAR GAIN (toy10, 02:50; run_expand_toy.sbatch ETAG=4 NFREE=0 FOCUS=1, children 14 -> 8 A):
98 children, all survive (cross-half 0.95 median at <= 12 A), 10/10 conformers own a state, one-state-per-conformer
purity mean 0.54 (parent 0.48): 18 0.36 -> 0.62, 51 0.64 -> 0.66, 27 0.49 -> 0.63, 28 0.49 -> 0.58, 37 0.51 -> 0.61,
46 0.62 -> 0.72, 48 0.45, 49 0.48; 1 0.38 and 10 0.31 slightly down. RECIPE: families by spectral bisection of the
parent's state graph; child board per family; focus ON (the children's maps are the noisy regime where focus helps);
free components OFF (they become outlier sinks without junk); survival judged at 12 A and coarser. v5 = v4 + pose
marg / devw / pnoise in the children (53893xxx); real data: run_expand_rnasep.sbatch (families of klem_full_npz, same
recipe, identification per family; rnasep_expand1).

LEVERS ROUND 2 RESULTS (toy 50k, 03:10; baseline 0.78 / 0.32 / 0.34 / closed 0.80, 51 0.45, 18 0.20): pose-marginalised
M-step (7 poses, E + M): 0.348 / 0.369 / 0.84, 51 0.53 (+0.01 over E-step-only marginalisation); deterministic annealing
T0 = 3: neutral (0.321 / 0.347 / 0.79); noise-weighted M-step: 0.323 / 0.351 / 0.82 (tiny); COMBO2 (pmarg E + M, pnoise E +
M, devw): 0.357 / 0.393 / closed 0.855, 51 0.61, 18 0.29, 10/10 own a cell, one-to-one mean 0.43 = the best single-board
configuration at 50k (+0.05 purity, closed conformer 51 +0.16). These pose/noise levers DO add (unlike the earlier
set). Full-scale confirmation: toy10_klem_fullcombo2. 27-pose and shift variants still running (slow: 27x / 35x E-step).

FAMILY-WISE HANDOFF (toy10, 03:30; expand2b children -> k-means 10 per family -> refine3D_states 10 states per family,
8 its): family 1 (open): 27 0.57 -> 0.72, 28 0.49 -> 0.61, 48 0.38, 49 0.43; family 2 (closed): 18 0.49 -> 0.68, 51
0.61 -> 0.65, 37 0.56 -> 0.62, 46 0.51 -> 0.62, 1 0.19 -> 0.41, 10 0.21 -> 0.37. Best-per-conformer across families:
1 0.41, 10 0.37, 18 0.68, 27 0.72, 28 0.61, 37 0.62, 46 0.62, 48 0.38, 49 0.43, 51 0.65 = mean 0.55 (parent board 0.48,
flat 10-state handoff 0.43, expansion v4 alone 0.54). The pipeline board -> families -> child board -> ~10 states per
family -> refinement is the best identification so far; the closed 18 reaches 0.68. Poses after the family refinements
3.75 / 3.80 deg (no gain in poses; the gain is the assignment). Lever 5 (family pose feedback) launched from these.
UMAPs (figs/fig_toy10_umap_50k_baseline_vs_combo2.png, fig_toy10_umap_expand4_families.png): the best single board
at 50k is more articulated than the baseline (27/28 arm, 10/37 lobe, 51 lobe); inside family 2's child board the closed
conformers 18 and 51 are separate lobes (in-plane 15-NN 0.48 vs chance 0.25).

FAMILY POSE FEEDBACK = BEST RESULT (toy10, 04:40; run_family_feedback.sbatch: each family's 10-state refinement poses
-> child board on the family, no free comps, pmarg + devw + pnoise): merged one-state-per-conformer purity mean 0.62
(min 0.43), 10/10: 18 0.65, 51 0.73, 27 0.76, 28 0.70, 37 0.67, 46 0.78, 1 0.48, 10 0.50, 48 0.43, 49 0.50 (parent board
0.48, expansion v4 0.54, family-wise handoff 0.55, true-pose single board 0.61). The refined poses' global error is
unchanged (3.75-3.80 deg): the gain comes from poses made consistent with the state references. Expansion v5 (pose
marg inside children WITHOUT feedback): 0.51, 18 0.41 -> pose marginalisation hurts the children unless the poses are
refined. Full-scale combo2 (pose-marg M-step added): 0.436 = fullbest; the M-step marginalisation adds nothing at 200k.
27 poses: 0.367 vs 0.360 (7 poses); shifts: nothing. RECIPE (toy): parent board -> families -> child boards (focus,
no free comps) -> k-means 10 per family -> refine3D_states 10 states -> child boards again at the refined poses.
REAL DATA: level 2 of the closed family (klem on the 262k refined-state-1 particles, k-means 3, refine3D_states):
sub-states = {18 (+0.18)}, {46 (+0.22)}, {10 (+0.29)} with closed-arm occupancy 0.68 / 0.68 / 0.61 (level 1: 0.66):
the closed arm does not purify at level 2 on real data (yet); real expansion family 1 identification (child board with
focus): 37 +0.79, 27 +0.76, 28 +0.71, 10 +0.66, 1 +0.55, 51 +0.48, 18 +0.46 (parent 0.35 mean; 37 0.64) = the child
boards match the deposited conformers far better; family 2 pending.

**How to apply:** this is the model class to develop; next levers if the toy runs confirm it: more cells / 3-D grid,
a bag-level (sub-class) E-step for the noise, a final local multi-state pose refinement, then Fortran. Related:
[[toy10-latent-night-2026-10-15]], [[toy10-model-class-verdict-2026-10-15]], [[prototype-in-python-before-fortran]].

## 2026-10-15 evening: one-pass tree realignment (user: "with this tree you could realign things in one pass")
- The tree = the pose refinement: each family's 10-state refine3D_states realigns its particles once; the merged
  refined tables (lat/fb_fam1,2/refined_state1.txt, state=1 everywhere) are ONE realigned pose table for the toy.
- run_tree_realign_toy.sbatch (job 53897503): merged table -> parent regrown at realigned poses (nfree 2) -> split ->
  children at realigned poses with the feedback recipe (NFREE=0 PMARG=3 DEVW=0.05 PNOISE=1, default 20->10 schedule,
  no focus) -> merge vs the 0.62 old-family feedback boards; LOOP=1 adds the second pass (k-means 10 -> refine3D_states
  started from the matched old family's refined project via PROJSRC -> feedback boards -> merge).
- run_tree_level3_toy.sbatch (job 53897504): each feedback child board bisected again -> 4 grandchildren at the
  family's refined poses (focus 1, PMARG 3, 12->8 A) -> merged purity vs level 2 (targets: 48/49 at 0.43/0.50, 1/10 at
  0.48/0.50).
- Real data: run_tree_realign_rnasep.sbatch (merged meta_realigned.npz from expand1/fam{1,2}_pipeline/refined_ptcl3D.txt
  -> parent INIT=npz NFREE=2 DEVW 0.05 PMARG 3 PNOISE, 20->8 schedule -> run_expand_rnasep at META=realigned ->
  identification); a local waiter submits it once both family pipelines print DONE. run_klem_rnasep.sbatch and
  run_expand_rnasep.sbatch now take META=; run_expand_toy.sbatch passes TABLE/PMARG/DEVW/PNOISE; run_state_refine_toy
  takes PROJSRC.
- The real best run (klem_full_best) = DEVW + 20/14/10/8 schedule, no focus (log markers), so the realigned real parent
  uses that recipe.
- Literature levers implemented 2026-10-15 ~18:40 (user: "try these on the toy in the meantime", in order): (1) annealing
  split criterion: klem.py now saves klem_halves.npz (even/odd cell deviations at <= 12 A, soft populations);
  anneal_split.py splits recursively while the cross-half reproducible variance along the family's principal mode
  exceeds tau x the noise along it (Rose's critical temperature), number of families from the data; run_expand_toy
  SPLITMODE=anneal; job run_levers_anneal.sbatch (parent re-run klem_junk_h, anneal vs spectral-2 control, children
  at realigned poses). (2) density modes: density_modes.py = extreme deconvolution mixture (per-particle posterior
  covariance S_i) + mean-shift modes + noise-aware assignment, and Richardson-Lucy + scale-space modes; vs k-means /
  argmax cells / undeconvolved mixture; job run_levers_density.sbatch on klem_junk, klem_fb_fam1/2. (3) OT-balanced
  E-step: klem.py --balance lambda (--balance_until k): Sinkhorn in log space, grid columns -> equal shares of the grid
  mass, free comps unconstrained, rb ~ exp(L + lambda v); env BAL/BALU; job run_levers_balance.sbatch (lambda 1, 0.5).
  Smoke test job run_levers_smoke.sbatch (synthetic 3 families / 6 modes) gates the submissions (waiter chain).
- LEVEL 3 RESULT (job 53897504, 44 min): each feedback child board bisected again, 4 grandchildren at the refined poses
  (focus 1, PMARG 3, 12->8 A): merged purity 0.66 (level 2: 0.62), 10/10 owned; 1 0.48->0.60, 10 0.50->0.63,
  18 0.65->0.75, 48 0.43->0.46, 49 0.50->0.54, 28 0.70->0.65. Caveat: 193 cells vs 98, the max-over-cells purity has
  a selection bias with more cells. Grandchild 1.2 (within-family similarity 0.59, the 48/49 side) has 15-NN below
  chance (0.378 vs 0.402): the near-identical pair does not separate at level 3 either. Depth helps 1/10/18.
- LEVER VERDICTS (toy, 2026-10-15 19:00): (3) balance: lambda 1 / 0.5 -> one-cell purity 0.51 / 0.51 (plain 0.48),
  populations 2559-3740 (plain 375-10152), closed AUC 0.83 / 0.86, d' 1.85 / 1.97; cheap, adopt 0.5 as the hedge for
  unequal real populations. (2) density modes: on the parent, deconvolved mixture = 12 modes at 0.42 purity (k-means
  10: 0.32, 6/10 owned; argmax cells 0.48 over 49); on balanced parents 15-16 modes at 0.45-0.49 with 9-10/10 owned;
  on the PMARG children the posteriors are too broad (sd 1.3 cells) and the deconvolution collapses to 5-6 modes; RL +
  scale-space finds 2-4 modes only (useless). Mode reading >> k-means, ~ cells with 3-4x fewer states. (1) anneal:
  root ratio 12.6, every node ratio > 3 -> a motion continuum never stops splitting (tau=1 is no stop rule); 8 families
  of 10-28k particles, pairs grouped (48/49, 27/28) but 7x7 children have 0/49 reproducible cells (~200 particles per
  cell) -> merge empty. Follow-ups running: anneal2 = same split with boards sized to the family (PERCELL 1500,
  job 53899002); balanced realigned tree (BAL 0.5 parent + children at realigned poses, job 53899003).
- REGROWN TREE (job 53897503): parent at realigned poses one-cell purity 0.57 (consensus parent 0.48), closed AUC 0.92,
  15-NN 0.47; families = the old ones (Jaccard 0.86/0.85); regrown children merged 0.62 = old tree 0.62 (51 0.73->0.60,
  46 0.78->0.82): iterating the tree at fixed poses saturates; second pass (refine again) pending.
- Spectral 2-way control on the re-run parent (klem_junk_h) with feedback-recipe children at realigned poses: 0.62
  (1 0.47, 10 0.54, 18 0.67, 27 0.78, 28 0.68, 37 0.70, 46 0.83, 48 0.45, 49 0.48, 51 0.62). 0.62 is the plateau of
  "2-way tree at realigned poses" whatever the parent (consensus / realigned / re-run); only depth (0.66) moved it.
  anneal_g (4x4 boards, ~1500 per cell): cells still 0-1 of 16 reproducible at 0.5 (median 0.37): small homogeneous
  families have deviations near the noise; survival gate, not particles per cell, is what fails there.
- anneal_g loose merges (all cells / 0.3): 0.58 / 0.57 over 107 / 90 states (48/49 0.41/0.43): the 8 annealed
  families on small unfocused boards are WORSE than the 2-way tree (0.62). anneal3 (job 53899260) = same families with
  the level-3 recipe (focus 1, 12->8 A, PMARG 3, ~1000 per cell) is the last word on lever 1.
- BALANCED REALIGNED TREE (job 53899003, BAL 0.5 parent + children at realigned poses): parent one-cell purity 0.61
  (unbalanced realigned 0.57, consensus 0.48; 51 0.74, 18 0.70, 27 0.71, 37 0.72, 28 0.66, 46 0.68), closed AUC 0.90,
  populations 1341-5274; children merged 0.64 (2-way plateau 0.62): 27 0.82, 28 0.75, 51 0.79, 46 0.79, 37 0.72,
  18 0.68, 10 0.54, 1 0.37 (lost), 48/49 0.44/0.46. Balance compounds with poses; the parent alone now reads as well
  as the true-pose board (0.61). Next: level 3 on this tree with BAL 0.5 (stack balance + depth).
- SECOND REFINEMENT PASS (realigned tree, job 53897503): k-means 10 on the regrown children -> refine3D_states again
  (poses unchanged in angle, 3.75-3.81 deg) -> feedback boards: merged 0.64 (first pass 0.62): 1 0.52, 10 0.54,
  18 0.68, 27 0.79, 28 0.72, 37 0.71, 46 0.80, 48 0.45, 49 0.50, 51 0.68. refine3D_states' own hard labels stay poor
  (0.37 / 0.25 per family): its value is the poses. Toy ledger at realigned poses: 2-way tree 0.62, + second pass 0.64,
  + balance 0.64, + level 3 0.66; 48/49 never above 0.50 (floor pair).
- LEVER 1 FINAL (anneal_f, job 53899260): the 8 annealed families with the level-3 recipe (focus, 12->8, PMARG, 3x3-5x5
  boards) all survive (median reproducibility 0.72-0.82, 157 states) but merge at 0.51 (2-way tree 0.62; 48/49
  0.47/0.51). Fine splits decided on the parent's 10 A board cut across conformers (each family spans 3-7 conformers)
  and the children cannot undo them; the 2-way split commits only to the coarsest division. VERDICT: the annealing
  criterion = a direction finder equivalent to the spectral bisection, no stop rule on a continuum, and multi-way
  parent splits hurt. Depth belongs to the children's own boards (one bisection per board, level 3 = 0.66). Closed.
- LEVER 2 FINAL (density2, job 53899621): on the balanced realigned parent, argmax cells 0.61 (49 states, recall
  0.19), k-means 10 0.40 (7/10), plain mixture modes 0.56 over 19 (9/10), DECONVOLVED modes 0.55 over 18 with 10/10
  owned and recall 0.34; unbalanced realigned parent: plain modes 0.61 over 18, deconvolved 0.53 over 14 (10/10,
  recall 0.38) vs cells 0.57. VERDICT: read states as mixture modes at the parent (= cells' purity with 3x fewer
  states, 10/10 coverage, double the recall); k-means is the wrong readout. On PMARG children (posterior sd 1.3
  cells) the deconvolution over-merges (4-7 modes): use plain modes there. density_modes.py is the readout tool.
- BEST TOY RESULT (job 53899620, level 3 + balance on the balanced realigned tree): merged 0.71 over 196 states,
  10/10 owned: 1 0.71, 10 0.72, 18 0.88, 27 0.79, 28 0.58, 37 0.80, 46 0.82, 48 0.47, 49 0.50, 51 0.88 (level 3
  without balance 0.66 over 193 states, like for like). Full stack: consensus parent 0.48 -> 2-way tree 0.62 ->
  realigned poses + balance 0.64 -> level 3 + balance 0.71. RECIPE TO PORT: balanced E-step (0.5), tree of one
  bisection per board with focus at the fine band, pose refinement per family between levels (poses, not labels),
  density-mode reading at the parent. 48/49 (0.67-similar) stay at the per-particle floor at every level.
- HEADROOM (job 53901705, 2026-10-15 22:20): hard family splits keep 0.88-0.89 of particles in their conformer's
  family at level 2 (conformer 1 only 0.70) and 0.79-0.81 at level 3 (10: 0.55, 18: 0.52 inside family 1's split):
  0.89 x 0.80 = 0.71 = the level-3 purity -> the family split IS the limit. Oracle (true refs) at the REALIGNED poses:
  0.613 / 0.618 / 0.619 at 10 / 8 / 6 A (consensus poses 0.58-0.59); nothing below 8 A (signal fraction < 1 % per
  shell beyond ~12 A); 48 recall 0.30, 49 0.53. Levers running: soft families (klem.py --weights = parent family
  responsibility in the M-step; expand_families softsplit / merge --weights; run_klem2.sbatch, run_tree_level3c_toy
  = soft level 3 on the hard children (53901829), run_tree_soft_toy = full soft tree (53901830)); leaf pose
  feedback (run_tree_leaf_feedback_toy, 53901831); final kernel width 0.5 at level 3 (53901832); toy tree with the
  232 A mask (53901722). Band lever (6 A) DEAD.
- LEVER: final kernel width 0.5 instead of 0.7 at level 3 (job 53901832): merged 0.71 = the 0.7 run (0.71), per
  conformer within 0.03. NEUTRAL; keep 0.7.
- USER DIRECTION (2026-10-15 23:30): get ALL states out before per-state refinement (refinement parked); keep developing
  on the toy; suspects a missing PRIMARY COMPOSITIONAL state. Running: toy10c (12 states incl. 61 deleted domain,
  62 added density) through the full current recipe (run_toy10c_tree.sbatch, job 53902677: parent random init ->
  families -> children -> per-family refinement -> feedback children BAL 0.5 -> level 3); the real family 2 (197k,
  nothing deposited) tested as a compositional state (fam2_composition.sbatch: mass / missing / extra density vs the
  consensus). Tree questions: realignment happens ONCE (family-wise refine3D_states between level 1 and 2, poses
  only); the leaf feedback job (53901831) tests re-refining at the leaves; "traversal back up" = regrow at refined
  poses (parent 0.48 -> 0.61 with balance, children unchanged 0.64) and global re-assignment at the leaves
  (run_tree_upglobal_toy.sbatch 53902671, --final_all + merge --bymaxr) and the soft trees (53901829/30: soft level 2
  = 0.64, neutral so far). Scripts for toy10c: run_expand_toyc / run_state_refine_toyc / run_tree_level3b_toyc
  (TAG/SIMDIR/PJ plumbed copies; originals untouched because jobs were running them).
- SOFT LEVEL 3 on the hard level-2 children (53901829, 1h17): merged 0.71 = hard level 3 (0.71); per conformer
  1 0.71->0.75, 46 0.82->0.85, 28 0.58->0.61, but 18 0.88->0.84, 51 0.88->0.79, 10 0.72->0.70. NEUTRAL: weighting the
  unions does not repair the stranded 20 % by itself. Pending: the full soft tree (53901830) and the global
  re-assignment at the leaves (53902671).
- TOOL count_states.py (2026-10-15 23:50): distinct reproducible states across boards = surviving cells (cross-half
  >= 0.5) clustered by cross-half deviation similarity (average linkage, thresholds 0.6-0.9) from klem_halves.npz;
  with --gt reports conformers represented, with cell_matches.txt the deposited conformers identified per cluster.
  run_count_states.sbatch: toy level-3 leaves (calibration vs 10 conformers) and real family 1 (klem_fb + l3 sub1/2).
  Purpose: the unsupervised "how many states" readout for the landscape goal (paper: 55 apo conformers).
- STATE COUNT CALIBRATION (count_states, 53902704): toy level-3 leaves (196 surviving cells, 4 boards): similarity
  0.6 -> 22 states (8/10 conformers), 0.7 -> 32 (9/10), 0.8 -> 48 (10/10), 0.9 -> 161; self-reproducibility median
  0.90, between-cell similarity median 0.13. The counter over-splits a discrete conformer ~2-3x at 0.6-0.7 (46 and 49
  split most). REAL family 1 (klem_fb + l3 sub1/2, 147 cells): 0.6 -> 20 states (median 29k particles), 0.7 -> 27,
  0.8 -> 39, 0.9 -> 66; self-reproducibility 0.99. With the toy's factor, family 1 holds ~10-15 distinct states at
  this resolution and subset size (paper: 55 apo conformers). Rerun with deposited identification per cluster and
  family 2's leaves = run_count_states_real.sbatch.
- POSE MARGINALISATION A/B (toy, 53902669, 29 min for children + level 3): WITHOUT --pose_marg: children 0.64 (with:
  0.64), level 3 0.69 (with: 0.71; 18 0.87, 27 0.85, 51 0.84). Worth 0.02 at the leaves at 7x the E-step cost.
  VERDICT: off by default (real runs 7x faster; the realignment parent had run 3 h with it); optional at the leaves.
  Real realignment resubmitted without it (PPMARG=0 CPMARG=0).
- TOY10C (compositional, 12 states) first stages (job 53902677, random init): parent 9/12 owned, 0.26; children at
  consensus poses (focus, no free comps) 11/12 owned, mean 0.42: deleted domain 61 (12 % mass) 0.44 (3.6x base),
  added density 62 0.32 (4.3x base), conformers 0.31-0.55. The tree DOES give the compositional states their own
  states at the second level; purity moderate; feedback + level 3 pending. PCA-init parent (53902719): 11/12 owned,
  0.37 (61 0.50, 62 0.27) - better than random (0.26).
- LEAF POSE FEEDBACK (53901831, 1h45): 10-state refinements seeded from each level-3 grandchild (their hard labels
  0.14-0.30: too few particles per state at ~4k), poses merged (angle 3.66 deg vs 3.81), grandchildren rebuilt:
  merged 0.71 = before (10 0.72->0.76, 1 0.71, 28 0.58->0.57, 49 0.50). NEUTRAL: pose feedback pays at the family
  level (0.48 -> 0.62) and not again at the leaves. Toy recipe ledger now closed: poses (family level), balance,
  depth; neutral: soft level 3, leaf feedback, sharper kernels, second family pass (+0.02), PMARG (+0.02 at 7x).
- FULL SOFT TREE (53901830, 2h): soft level 2 = 0.64 (hard 0.64), soft level 3 = 0.71 (hard 0.71): NEUTRAL at
  both levels. Weighting the unions by the parent's family responsibilities does not recover the stranded particles;
  the only remaining test of the family-boundary loss is the global re-assignment at the leaves (53902671).
- UP-TRAVERSAL global re-assignment (53902671): each rebuilt grandchild scored ALL 160k particles (--final_all, ~12
  min each); merge --bymaxr: 0.73 (hard 0.71) over 165 states, BUT 18/27/48/49 at "1.00" = tiny pure states (size-
  blind purity); 1 0.40, 46 0.50. Re-scored with --minsize 1000/2000 (expand_families merge --minsize, job
  toy10_merge_minsize) before any verdict. TRAP: the one-state-per-conformer purity rewards tiny states.
- UP-TRAVERSAL VERDICT (size-floored, job toy10_merge_minsize): hard level 3 = 0.71 (all) / 0.69 (states >= 1000
  particles, 67 states, 101k) / 0.48 (>= 2000: only 12 states left, cells hold ~800 median); global re-assignment by
  max responsibility = 0.73 (all) / 0.45 (>= 1000) / 0.35 (>= 2000): its gain was tiny pure states; cross-board
  responsibilities are not comparable and particles collapse into over-confident cells. NEGATIVE. All three "up"
  variants (regrow, soft families, global re-assignment) fail to recover the stranded particles: the family-boundary
  loss is a limit of the per-particle evidence at this SNR, not of the assignment rule. The hard tree stands; quote
  purities with a 1000-particle floor from now on (0.69 at level 3).
- TOY10C PCA-init children (53902719, consensus poses, focus): 11/12 owned, mean 0.51 (random init 0.42; toy10's
  equivalent stage 0.54): deleted domain 61 = 0.72 (the purest state), added density 62 = 0.35 (4.7x base),
  conformers 0.38-0.57. Composition is carried; the 6 %-mass added density is the weak one. Feedback + level 3 pending.
