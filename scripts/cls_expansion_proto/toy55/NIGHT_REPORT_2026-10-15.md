# Latent-variable night on the RNase P toy (2026-10-15, 03:00-08:00)

Question: which latent model class can represent the RNase P conformational landscape, and can it be inferred per
particle at RNase P signal with realistic poses. Everything below is scored against truth on a toy built from the ten
deposited conformers (real poses and CTF, calibrated SNR, 20 % junk) and on a second toy with two compositional
variants (a conformer with its mobile domain deleted, 12 %; a conformer with an added surface density, 8 %).

## Verdict in five lines

1. The landscape is not low-rank in volume space (4 to 5 dimensions per band) but it is low-rank in deformation space
   (2 to 3 modes): a deformation latent plus per-component occupancies is the right model class; a linear volume
   basis is not.
2. The information is there per particle: with the true maps as references a single particle picks its conformer with
   accuracy 0.53 to 0.59 (chance 0.15) and the posterior coordinate tracks the truth at Spearman 0.84.
3. The linear latent is the bottleneck, and its band matters: lp 20 is near chance (0.17), lp 12 reaches 0.25 with one
   component ordering particles at 0.63.
4. EM against the particles climbs from any seed that is not a single junk axis: supervised windows reach Spearman
   0.76 and purity 0.31, truth-free k-means seeds reach purity 0.28 in three rounds and are still rising; one canonical
   map warped along a learnt flow mode replaces eight reconstructed references.
5. What stops it short of the ceiling is the one-dimensional state space and the 2-iteration per-window pose step
   (which hurts); a full local multi-state refinement instead improves poses to 3.5 deg.

## 1. Model class (truth maps only)

Unexplained fraction of each conformer's difference from the consensus, mean over conformers, per band
(117-30 / 30-20 / 20-12 / 12-8 A):

| model | 117-30 | 30-20 | 20-12 | 12-8 |
|---|---|---|---|---|
| linear volume basis, K=1 | 0.70 | 0.75 | 0.76 | 0.80 |
| linear volume basis, K=5 | 0.65 | 0.64 | 0.64 | 0.67 |
| global shift only | 1.00 | 1.00 | 1.01 | 1.00 |
| smooth deformation field (elastic 3 px) | 0.37 | 0.38 | 0.49 | 0.60 |
| deformation + occupancy | 0.34 | 0.36 | 0.49 | 0.59 |

- The deformation space of the ten conformers has effective dimension 2.6 (first mode 57 %); the volume-difference
  space needs 3.5 to 4.7 dimensions per band. Displacements reach 5 to 12 A. The two closed conformers are the
  worst-explained (largest rearrangement).
- Regularisation decides: elastic sigma 6 px leaves 0.61 to 0.65 unexplained, 8 to 14 px 0.7 to 0.9. The field must be
  weakly regularised or piecewise; a stiff global field collapses to the same "breathing" as the linear basis.
- Composition: the deleted-domain state is unexplained by deformation (0.64) and explained by deformation plus a
  multiplicative occupancy (0.08 at low resolution, occupancy map finds the hole). The added density is not helped by
  an occupancy clipped at one: addition needs an additive component. The deformation coordinate of each variant equals
  its parent conformer's, so conformational and compositional axes separate cleanly.

## 2. Per-particle inference ceiling (oracle: true conformer maps as references)

| band | accuracy over 10 conformers (chance 0.15) | Spearman vs conformer axis | closed recall |
|---|---|---|---|
| 20 A | 0.53 | 0.81 | 0.70 |
| 14 A | 0.57 | 0.83 | 0.74 |
| 10 A | 0.58 | 0.84 | 0.76 |
| 8 A | 0.59 | 0.84 | 0.77 |

Compositional toy (12 states, chance 0.12): 0.49 to 0.53, closed 0.71 to 0.77, deleted domain 0.56 to 0.58, added
density 0.41 to 0.43. Most of the information sits at 20 A and coarser. Pose error after refinement: median 4.0 deg,
90th percentile 7.4 deg, shifts 1.4 px.

## 3. Linear latent per particle (flex_pca, 16 components)

At lp 20: 15-NN accuracy 0.167 (chance 0.149); no single component correlates with the conformer axis above 0.20; the
best supervised linear combination reaches Spearman 0.49; closed-state AUC 0.67; the first component is a junk and
contrast axis (junk AUC 0.81); delivered states: one of 158k plus three tiny ones.

At lp 12 the same model is far more informative: one unsupervised component orders particles at Spearman 0.63, the
supervised combination reaches 0.71, 15-NN accuracy 0.25, pairwise separation 1.54 (lp 20: 0.88), closed AUC 0.75.
The lp 20 basis is simply too coarse; the discriminative information sits at 20 to 12 Å, where the oracle gains as
well. A truth-free EM seeded by k-means on the lp 12 latent is running.

## 4. Coarse-to-fine windows with a cross-half E-step (EM over an ordered latent)

- Seeded from PC1: stays at chance through three rounds (window purity 0.15 to 0.16).
- Seeded from the supervised linear combination (window Spearman 0.49, purity 0.21): one E-step at 20 A lifts the
  posterior-mean coordinate to Spearman 0.72; rounds at 14 and 10 A reach 0.74 and 0.75, purity 0.30, Hungarian
  accuracy 0.28, max posterior 0.81. EM climbs from a weak seed but converges slowly below the oracle because the
  eight references are mixtures along one coordinate while the conformers span two to three deformation dimensions.
- Per-window pose refinement (two shc iterations against a 25k-particle window map) made poses worse on random
  windows (median 4.0 to 6.5 deg) and slightly worse on coherent windows (4.0 to 4.7, then 4.3). Without the pose
  step the same seed reaches Spearman 0.77 at 14 A and 0.76 at 10 A with purity 0.31 and four of eight window maps
  matching their dominant conformer (one of eight with the pose step). Both plateau near 0.76: the one-dimensional
  window limit, not the pose step.

- Compositional toy, same supervised seed without the pose step: Spearman 0.68, 0.73, 0.72 over the rounds,
  purity 0.25. The deleted-domain state ends up dominating its own window (share 0.31) next to its parent conformer's
  window, and junk collects in one end window (65 %): an 8 % mass loss is a strong low-resolution signal that even a
  one-dimensional coordinate picks up.

## 5. Deformation latent on the window maps

The eight window maps form a one-mode deformation trajectory (87 % of the flow variance, monotone with the
coordinate). Nine references made by warping one canonical map along that mode discriminate particles as well as the
eight separately reconstructed maps (Spearman 0.746 vs 0.752, purity 0.28 vs 0.30) and reach the closed state better
(recall 0.42 vs 0.22). This replicates on the no-pose-step arm (grid Spearman 0.767 vs windows 0.758, closed recall 0.37 vs 0.28) and on
the compositional toy (0.741 vs 0.719, closed 0.36 vs 0.21). A second flow mode from one-dimensional windows is weak: on the noisier arm it is noise (pairwise conformer
separation 1.42 in 2-D vs 1.40 in 1-D), on the cleaner arm it adds a little (1.57 vs 1.33) and splits the conformers
into two groups orthogonal to the main axis. A second latent dimension has to come from a two-dimensional state
space in the EM itself, not from one-dimensional windows.

Figures: `figs/fig_model_class_s3.png` and `figs/fig_model_class_c3.png` (model class per band), `figs/fig_windows_lat_sup3_r3.png`
(window maps of the no-pose-step arm beside their dominant truth conformers; the end windows carry a coherent
difference, the middle ones are noise-dominated mixtures).

## 6. Sub-class pooling and a truth-free seed

- The pose-mean 2D expansion gives 1600 sub-classes (median 54 members) with majority-conformer purity 0.27; junk is
  isolated almost perfectly (84 % of junk in junk-majority sub-classes).
- Pooling the linear latent per sub-class doubles its discriminability (5-NN accuracy 0.34 vs 0.17 per particle,
  pairwise separation 1.74 vs 0.88, Spearman up to 0.45 on the second to fourth components) but stays far below the
  per-particle oracle.
- Seeding eight states by k-means on the pooled latent (no truth) gives purity 0.19 with half of the closed particles
  in one state; three rounds of hard EM lift purity to 0.21, 0.22 and 0.23 (Hungarian accuracy 0.23, closed recall up
  to 0.57). A truth-free seed climbs slowly but monotonically, far below the supervised axis, because the k-means
  states partly follow view and junk structure.

### Truth-free seeds from k-means on the particle latent

Eight k-means states on the lp 12 latent (no truth anywhere) followed by three rounds of hard EM against the
particles: purity 0.22, 0.24, 0.26, 0.28, Hungarian accuracy 0.28, one state collecting 53 % of the closed particles at
purity 0.33. The same on the lp 20 latent: purity 0.18 to 0.27. Both are still rising at round three and already match
the supervised-seed window arms. The cross-half E-step does the work once the seed is not a single junk axis. Four
more rounds at 10 and 8 Å plateau at purity 0.29 (Hungarian accuracy 0.28) with the closed state crystallising in one
state (56 % of the closed particles at purity 0.41): eight hard states on a ten-conformer continuum plus junk saturate
there. On the compositional toy the lp 12 latent reaches 15-NN accuracy 0.21 (chance 0.12) with the
deleted-domain state at one extreme of its first component, the pooled seed climbs from purity 0.16 to 0.21 over three
rounds, and the k-means-seeded EM reaches purity 0.24 in three rounds with two states that are pure junk (95 to 98 %)
and one state owned by the deleted-domain conformer at 37 % purity: without any truth the EM isolates the junk and the
compositional state by itself, while the conformational states stay mixtures.

## 7. Discrete multi-state refinement (refine3D_states, four states, local poses, 20 to 8 Å)

From four random coarse windows: weighted purity 0.16 to 0.19, Hungarian accuracy 0.09 to 0.11, half of the closed
particles collected in one state at purity 0.16, poses slightly degraded (4.0 to 4.4 deg). Discrete EM from a random
seed moves little and mostly on junk and view structure. Seeded from the supervised coarse windows it keeps purity
flat (0.25), raises closed recall from 0.34 to 0.54 and improves the poses from 4.3 to 3.7 deg: a full local multi-state
refinement fixes poses once the states are coherent, but the discrete E-step does not sharpen the partition. From the
no-pose-step arm's windows: purity 0.27 to 0.29, closed recall 0.47 to 0.60, poses 4.0 to 3.5 deg (the best of the night).
On the compositional toy the same step improves poses (4.0 to 3.6 deg) and closed recall (0.45 to 0.58) while purity
slips from 0.21 to 0.19.

## 8. Status at 08:50

All streams complete. The per-group local flex runs produced no embeddings and were dropped.

## 10. Kernel-latent EM (the nonlinear model chosen over deformation fields), 15 October daytime

Model: a 7 by 7 latent grid of free cell volumes made smooth by a radial-basis kernel on the grid, two off-grid
components for junk, soft cross-half E-step from the band-limited image likelihood, soft kernel-smoothed M-step on a
coarse Fourier grid, bands 20, 14, 10 Å with the kernel annealed, initialised from two contrast-free components of the
lp 12 linear latent. flex_pca is the special case with a linear map from the latent to the volume; composition is free
because each cell is a free volume. Python only (`klem.py`), about 30 seconds per iteration on 50k particles.

- Synthetic 2-D family (6 px motion, present-or-absent blob, 10 % junk), random start: first axis follows the motion at
  Spearman 0.80, second axis separates the composition at AUC 0.975, 15-NN accuracy 0.76 (chance 0.11), junk 84 % in
  the free component.
- toy10, 50k particles, six minutes: Spearman 0.78 at the 14 Å stage (0.73 after the 10 Å stage, where 1k particles per
  cell are too few), 15-NN 0.32 (chance 0.15), cell purity 0.34, closed AUC 0.80, 99 % of the junk in the free
  components once they start junk-like. The vertical axis orders the ten conformers in the truth sequence and the
  closed conformer 51 separates along the second axis. This matches the best overnight EM in a hundredth of the time.
- Compositional toy, 50k: Spearman 0.70, 15-NN 0.24 (chance 0.13), purity 0.27; the deleted-domain and added-density
  states still sit on their parent conformers, so composition is not yet resolved by this configuration.
- Full scale (200k particles, 35 minutes): Spearman 0.80, 15-NN 0.37 (chance 0.16), purity 0.39, closed AUC 0.76; the
  conformer order along the main axis is exactly the truth order and the closed conformer 51 separates on the second
  axis. The compositional toy at full scale: 0.77 / 0.29 / 0.31, composition still unresolved. Oracle ceiling for
  comparison: coordinate 0.84, per-particle accuracy 0.59.
- With the simulation's true poses instead of the refined ones (median error 4 deg) the same 50k run reaches Spearman
  0.88, 15-NN 0.48 and purity 0.49 against 0.78, 0.32 and 0.34: pose error is the largest single loss, so pose
  refinement against the latent-decoded references is the next lever.
- Alternating with SIMPLE's local multi-state refinement (seven latent bins as states) recovers about half of that
  loss in one round: poses 4.0 to 3.7 deg, and the latent on the refined poses reaches Spearman 0.80, 15-NN 0.39,
  purity 0.40 and closed AUC 0.87 (same 50k setup before: 0.78, 0.32, 0.34, 0.80). More rounds are running.
- Random start: organises too (fitted 0.80 at full scale, purity 0.31) but less legibly. A 15-cell one-dimensional grid
  orders worse (0.55) than the smoothed 7 by 7; the spectral (Laplacian) initialisation brings no gain. Kernel-annealing
  brings nothing either. Extra free components seeded by k-means do crystallise conformer-enriched discrete groups
  (conformer 27 at purity 0.46 against a base rate of 0.15; the deleted-domain state at 0.28 against 0.12) but drain
  the continuum, so the grid metrics drop. A finer 10 by 10 grid is being tested for the compositional state.

## 11. First run on the real RNase P set (50k-particle gate, 15 October 13:25)

The kernel-latent EM on 50k real particles (consensus poses, latent initialised from the earlier flex embedding) gives
a latent that is structural, not nuisance: mutual information between latent cells and tilt, rotation or defocus bins
is 0.004 to 0.010, while the cells agree with the cls_expansion sub-classes within a parent class (within-parent NMI
0.22, median over 111 classes) and separate the S1-like soft-EM state from the S2/S4 pair along one axis. No 2D
information enters the loop. A second gate from a random start reproduces the same structure (within-parent NMI
0.21, nuisance NMI below 0.01). Five bins along the gate's main latent axis, reconstructed at full box (9k particles
each, 11 to 17 Å), resemble the deposited conformers in order: bin 1 conformer 37 (then 18, 48, 49), bin 5 conformer
28, with the bin 1 map anti-correlated with 27, 28 and 49 and bin 4 anti-correlated with 10, 18 and 37; the end-to-end
difference correlates positively (0.17 to 0.28 at 60 to 20 Å) with every (27, 28, 49) minus (10, 37, 18) conformer
difference. This is the same ordering as the toy's truth axis. The full 692k-particle run and its bins are in progress.

Sharing evidence within the cls_expansion sub-classes (a bag prior in the E-step) was tried on both toy and real data: on
the toy it hurts (the sub-classes are too impure), on the real data it raises the agreement with the sub-classes but also
pulls in their view component (tilt and rotation information triple), so it is not adopted in this form.

## 12. Afternoon: embeddings, reproducibility, the full run, alternation plateau, concurrent streams (16:50)

- Embeddings (`klem_umap.py`, figs/fig_rnasep_gate_umap.png, fig_rnasep_flex_vs_klem_umap.png): a UMAP of the 51-D
  responsibility vectors of the 50k gate is a two-dimensional sheet with z1 and z2 running across it; the two thin
  strands are the junk flowing into the free components; tilt and defocus are uniform over the sheet; the S1 soft-EM
  state sits in one lobe. The flex_pca latent of the same particles is a featureless ball whose only pocket is the junk.
- Reproducibility (`klem_match.py`, orientation-free matching of cell volumes): the two 50k gates agree at a median
  matched-cell deviation correlation of 0.43, with 43 % of particles within one cell of the matched cell and a soft
  per-particle cross-run correlation of 0.23 against 0.005 for permuted particles; the full 692k run against the 50k
  gate: 0.60, 67 %, 0.33. Per-particle positions are noisy, the manifold is reproducible and sharpens with particles.
- Full 692k run (2.9 h): 89 % on the grid; nuisance NMI 0.005 to 0.013; within-parent sub-class NMI 0.065 against a
  within-parent shuffle null of 0.026; soft-EM state NMI 0.056 against 0.000; S1 separates along z1. Correction: the
  50k gates' within-parent NMI of 0.21 was mostly finite-sample bias; every NMI is now printed with its null.
- Alternation plateau on the toy: rounds 2 and 3 leave the poses at 3.65 to 3.74 degrees and the latent at
  0.39 / 0.40 / 0.87, equal to the run without alternation. Local multi-state refinement against bin maps gives
  nothing beyond the first round. Poses are parked unless the real 50k alternation (running) shows movement.
- Composition clarified: nothing compositional has been tested on RNase P (all ten conformers are full-length). The
  toy's deleted mobile domain (12 %) stays on its parent cell because the parent's cell map is an 88:12 mixture that
  fits it better than any neighbour and the grid has no discrete direction. The assembly-intermediate regime (35 % of
  the mass missing, 15 % of particles; toy10d) is being simulated to test the claim that large deletions are easy.
- Running: centred bag prior and grid deconvolution on the toy and the real 50k gate, the real 50k pose alternation,
  the toy10d chain, the full-run 7 bins at full box, the random-start 692k replicate.

## 13. Evening: the ledger of levers, the state graphs, the deposited landscape, the handoff (22:30)

Objective reset by the user: the model's job is to identify every state; pose refinement is a separate, later step
(minor in-loop refinement acceptable). Readouts for that objective: `state_coverage.py` (toy, per conformer: cells
owned, purity, recall) and `deposited_match.py` (real: for each of the ten deposited conformers the best-matching state).

- Coverage now: toy10 full scale 10/10 conformers own a cell, one-to-one purity mean 0.48 (min 0.36); oracle per-
  particle precision 0.38-0.73 (figs/fig_state_toy10_coverage.png). Real: 9/10 deposited conformers identified at
  correlation >= 0.3 (48 missing; 18 and 46 with one state), mean best match 0.35.
- Levers tested and neutral or worse: 3-D grid, 10x10 grid, adaptive kernel width, grid deconvolution, learned state
  graph (over-smoothed first version; calibrated version equal to the lattice at iteration 6), zoom on particle bins.
- 2D coupling in the E-step: four forms (mean-log-likelihood bag prior, its parent-centred form, hard sub-class
  pooling, NNLS unmixing prior) all leak the view into the latent on real data and lower the bins' conformer
  correlations; the parent-centred unmixing (difference to the parent average, view cancels) is the last form, running.
- True poses at full scale: 0.91 / 0.59 / 0.61 against 0.80 / 0.39 / 0.42 with refined poses: the only large lever.
- State graphs (figs/fig_rnasep_full_state_graph*.png, fig_rnasep_state_graph_deposited.png, _joint.png): our 49
  states form a horseshoe with the closed family at one end; the ten deposited conformers form two anti-correlated
  families with tight pairs (27/28 at 0.67, 48/49 at 0.69); in a joint layout our states sit on the rim and the
  deposited conformers in the interior, 0.35 correlation at best. The lattice places states at extremes; the paper's
  conformers are interior points.
- Handoff test (running): unsupervised 10-cluster partition of the toy latent -> SIMPLE multi-state refinement ->
  per-conformer purity before and after. This is the "does the pipeline work" check before any 50-state refinement.
- Fortran port plan written: doc/implementation_notes/klem_fortran_port_plan_2026_10_15.md (master's state weight
  set with parent-state genealogy + the flex-softem E-step; the lattice kernel is the one new piece).

## 14. Night: Fortran kernels, levers, the handoff, and the expansion that works (03:00)

- Fortran OpenMP kernels for the E-step and M-step (`klem_kernels.f90`, ctypes from `klem.py`): identical results,
  50k run 336 s -> 77 s; a full 692k run in ~2 h instead of 9.
- Handoff (unsupervised 10-cluster partition of the toy latent -> SIMPLE multi-state refinement): purity 0.32 -> 0.43,
  8/10 conformers kept, poses 3.68 -> 3.56 deg: identify-then-refine works on a known case. refine3D_states fails at
  49 states (assignment table scales with states x orientations), so refinement is hierarchical: ~10 states per family.
- Real closed family: 3-state refinement from the full-run bins raises the conformer correlations by half
  (state 1 = 37/10/18 at 0.52/0.48/0.46) but the closed arm stays at 0.66 occupancy: 18 is a minority in its family;
  level 2 (klem inside the refined family at refined poses -> 3-state refinement) is running.
- Levers: pose marginalisation (+0.015 purity), deviation-FSC shell weights (+0.02), per-particle noise (+0.01),
  focus (+0.05 at 50k, neutral at full scale), 8 A band with weights (+0.015); tempered finish, mixing proportions,
  per-particle contrast, junk-free second pass, 10x10, 3-D grid, learned state graph: neutral or worse. The gains do
  not add beyond ~0.44 purity with refined poses (true poses 0.61).
- Free components turn into outlier sinks in junk-free runs (they absorbed 60 % of the closed conformer 51 in a second
  pass and in expansion children): none below the top level.
- EXPANSION v4 (families from the state graph, child boards at 14-8 A with focus, no free components, survival at
  <= 12 A): purity 0.48 -> 0.54, 10/10 states, closed 18 0.36 -> 0.62, 51 0.66, open family +0.1 each. The first
  clear model-side gain; v5 (pose levers inside the children) and the real-data version are running.
- All five remaining levers (27-pose neighbourhood, shifts, pose-marginalised M-step, annealing, noise-weighted M-step,
  family pose feedback) implemented and queued.

## 9. Where everything is

batch2 `/mnt/beegfs/elmlund/afan/routes/toy55/`: toy `sim/` and `simc/`, refined projects `proj/`, `projc/`, truth maps
`maps160/docked/` and `maps160/variants/`, model class `model_class_s3/`, `model_class_c3/`, arms `lat/` (PC1 seed),
`lat_sup2/`, `lat_sup3/` (supervised seed with and without the pose step), `lat/kseed/`, `lat/kseed20/` (k-means
seeds), `lat/pooled/` (sub-class-pooled seed), `latc*` the same on the compositional toy; logs `toy10*_*.out`;
`collect_results.sh` prints the digest. Scripts mirrored in `scripts/cls_expansion_proto/toy55/`.

## 15. One-pass tree realignment (submitted 2026-10-15 ~18:00)

The tree is itself the realignment: every on-grid particle has been refined once against its family's ten state
references. Two toy jobs test what the tree gives when that is used in one pass:

- `run_tree_realign_toy.sbatch` (53897503): the two families' refined pose tables merged into one; the parent board
  regrown on it, re-split into families, child boards at the realigned poses (the feedback recipe), merged purity against
  the 0.62 of the feedback boards on the consensus-pose families. Then a second pass (k-means 10 per regrown family,
  refine3D_states started from the previous refinement project, feedback boards, merge).
- `run_tree_level3_toy.sbatch` (53897504): a third tree level, every feedback child board bisected again into two
  sub-families with a grandchild board each (focus on, pose marginalisation, 12 to 8 A), merged purity against level 2.
  Targets the near-identical pairs 48/49 (0.43/0.50) and 1/10 (0.48/0.50).
- Real data: `run_tree_realign_rnasep.sbatch` is queued by a waiter once both family pipelines finish: merged realigned
  poses over the 692k, parent regrown (deviation weights, pose marginalisation, 20 to 8 A), families and child boards at
  the realigned poses, identification against the deposited conformers at each level.

## 16. Literature levers on the toy (submitted 2026-10-15 ~18:40)

Three levers from the literature search, each aimed at a measured failure, implemented and submitted on the toy:

- **Annealing split criterion** (Rose): `anneal_split.py` reads the even/odd cell deviations that `klem.py` now saves
  (`klem_halves.npz`) and splits a family while the cross-half reproducible variance along its principal mode exceeds
  the noise along it; the number of families comes from the data. `run_levers_anneal.sbatch` re-runs the parent, splits
  it both ways (anneal vs spectral two-way) and runs the children at the realigned poses with the same recipe.
- **Density-mode state reading** (RECOVAR): `density_modes.py` fits an extreme-deconvolution mixture to the posterior
  means with each particle's posterior covariance, finds the modes by mean-shift, assigns particles with the noise-aware
  posterior; a Richardson-Lucy and scale-space variant alongside. `run_levers_density.sbatch` on the parent and the two
  feedback child boards, against k-means, argmax cells and the undeconvolved mixture.
- **OT-balanced E-step** (Sinkhorn): `klem.py --balance` constrains the grid cells to equal shares of the grid mass,
  free components unconstrained. `run_levers_balance.sbatch` runs the parent at lambda 1 and 0.5 against the plain parent.

**Level 3 result (18:40).** Four grandchildren, merged: purity 0.66 against 0.62 at level 2, all ten conformers owning a
state. Conformers 1, 10 and 18 gain most (0.48 to 0.60, 0.50 to 0.63, 0.65 to 0.75); 48/49 stay at 0.46/0.54 and the
grandchild board holding them reads at chance. Caveat: 193 cells against 98, so the best-cell purity carries a larger
selection bias at level 3.

**Lever results (19:00).** Balance: one-cell purity 0.51 at both strengths against 0.48 plain, cell populations
2559-3740 against 375-10152, closed separation 0.83-0.86 against 0.74. Density modes: 12 modes on the parent at 0.42
purity against k-means ten at 0.32 (six of ten owned) and the 49 argmax cells at 0.48; on the balanced parents 15-16
modes at 0.45-0.49 with nine or ten conformers owned; on the pose-marginalised children the posteriors are too broad and
the deconvolved density collapses to five or six modes. Annealing split: every node of the board has reproducible
variance well above its noise (root ratio 12.6, all nodes above 3), so on a motion continuum the criterion gives the
direction but no stopping rule; the eight small families it produced starve 7x7 children (no reproducible cells).
Follow-ups: the same eight families with boards sized to their mass, and the balanced E-step through the realigned tree.

**Regrown tree (18:40).** The parent at the realigned poses reads 0.57 one-cell purity against 0.48 at consensus poses,
closed separation 0.92; its families are the old ones; the regrown children merge at 0.62, the same as before. Iterating
the tree at fixed poses saturates; the second refinement pass is pending.

**Balanced realigned tree (19:08).** With the balanced E-step at half strength, the parent at the realigned poses reads
0.61 one-cell purity (0.57 unbalanced, 0.48 at consensus poses) and its children merge at 0.64 against the 0.62 plateau;
27, 28, 51 and 46 reach 0.75 to 0.82, conformer 1 drops to 0.37, 48/49 stay at 0.44/0.46. The two levers that moved
the toy are poses and balance at the parent, depth at the leaves; next is their combination.

**Second refinement pass (19:14).** Refining the regrown families again and rebuilding their boards gives 0.64 against
0.62, with the poses unchanged in angular error; the refinement's hard state labels themselves stay at 0.25 to 0.37
per family, so its contribution is the poses. Toy ledger at realigned poses: two-way tree 0.62, second pass 0.64,
balanced E-step 0.64, third level 0.66; 48/49 never pass 0.50.

**Annealing split, final (19:17).** With the focused fine-band recipe every cell of the eight annealed families
reproduces (0.72 to 0.82), yet the merge reads 0.51 against 0.62 for the two-way tree. Fine splits decided on the
parent's 10 Å board cut across conformers and the children cannot undo them; the two-way split commits only to the
coarsest division and leaves the rest to the children. The criterion is a direction finder equivalent to the spectral
bisection and offers no stopping rule on a continuum. Closed; depth comes from one bisection per child board.

**Density modes, final (19:45).** On the balanced realigned parent the deconvolved mixture reads 18 modes at 0.55
purity with all ten conformers owning a mode and recall 0.34, against the 49 argmax cells at 0.61 purity and recall
0.19, and k-means ten at 0.40 with seven owned. On the unbalanced realigned parent the plain mixture's 18 modes reach
0.61, the cells' level. The mode reading is the state readout at the parent level; on the pose-marginalised children
the posteriors are too broad for the deconvolution, which over-merges to four to seven modes, so the plain modes
serve there.

## 17. Best toy result: level 3 with the balanced E-step (19:58)

Four grandchildren of the balanced realigned tree, each with the balanced E-step, focus mask and the 12 to 8 Å band,
merge at 0.71 with all ten conformers owning a state: 1 0.71, 10 0.72, 18 0.88, 27 0.79, 28 0.58, 37 0.80, 46 0.82,
48 0.47, 49 0.50, 51 0.88. The same depth without balance gave 0.66 over the same number of states. The full stack on
the toy: consensus parent 0.48, two-way tree 0.62, realigned poses with balance 0.64, third level with balance 0.71.
The near-identical pair 48/49 stays at the per-particle floor at every level. This is the recipe to carry to the
real data and to the Fortran port: balanced E-step, one bisection per board with the focus mask at the fine band, pose
refinement per family between levels, density-mode reading at the parent.

## 18. Real data: all ten deposited conformers identified (20:36)

Family 1's pipeline on the real RNase P set (410k particles): ten-state refinement from the child board, then the
child board again at the refined poses. That feedback board identifies all ten deposited conformers, mean best-match
correlation 0.46: 1 +0.54, 10 +0.70, 18 +0.50, 27 +0.79, 28 +0.69, 37 +0.76, 46 +0.41, 48 +0.55, 49 +0.31, 51 +0.52.
The child board at consensus poses identified seven (0.44) and the full 692k parent eight (0.33). The real third
level on this board, with the balanced E-step, is job 53901575; the real tree realignment follows.

## 19. The ten deposited conformers as real-data maps (20:50)

Every cell of family 1's feedback board labelled by its best deposited conformer, particles inheriting the label, one
full-box reconstruction per conformer at the refined poses: ten maps at 4.8 to 7.8 Å (conformer 46 from 4.4k particles
at 7.8 Å, 37 from 112k at 4.8 Å). Scored with the set mean subtracted, every map is best matched by its own deposited
conformer at 60 to 20 Å (mean +0.54; 27 +0.82, 37 +0.79, 10 +0.61) and nine of ten at 30 to 12 Å. Conformer 18's map
has closed-arm occupancy 0.76, the purest closed state on the real data. Scoring trap recorded: subtract the set mean,
or the common real-versus-toy offset makes one conformer win everything.

## 20. The closed arm, the mask, and the next levers (22:30)

The deposited closed conformer's arm reaches 114 Å from the box centre; klem's default reference mask (radius 100 Å)
and the 100 Å mask on the sharpened maps clipped its tip at every stage. The closed boards now run with a 232 Å
mask diameter, the sharpened maps were remade at 116 Å, and the toy's best tree is being rerun with the wide mask.
Ranking the closed particles by confidence gives arm occupancy 0.81 / 0.77 / 0.73 for the top 10k / next 10k / next
20k, so purity is only part of the shortfall. In the deposited maps the ordered mass is within 5 to 12 % of conformer
18 for every conformer: the arm moves rather than disappears, with some disorder at the 48/49 end.

Headroom: the hard family splits keep 0.89 of the particles in their conformer's family at level 2 and 0.80 at
level 3; their product is the level-3 purity, so the split is the limit. The oracle at the realigned poses reads
0.613 / 0.618 / 0.619 at 10 / 8 / 6 Å: the band lever is dead. Running: soft families (weights from the parent's
responsibilities in the M-step, every particle scored in every leaf), pose feedback from the leaves, a sharper leaf
kernel, and the wide-mask tree.

**The closed arm, resolved (22:30).** On the toy, pure conformer-18 particles reconstructed at their true poses measure
0.85 arm occupancy by our metric, 0.81 at the realigned poses, 0.80 at consensus poses; the real closed subsets measure
0.78 to 0.81. The real closed state therefore matches a perfectly pure closed state at our pose accuracy; the metric's
own ceiling at 8 Å is about 0.85. What made the helix look short was the 100 Å mask at the reference and viewing
stages (fixed) and pose blur at the arm tip, which the closed-state local refinement now running addresses.

## 21. What "closed" is, the paper's pipeline, and the box (23:00)

Degenhardt et al. (Nature, 9 Oct 2026): closed = the tetraloop-tetraloop-receptor interaction between P10.1 and P12
in the S-domain, with P12's interdigitated T-loop motif properly folded; open conformers lack it. Closed is about 13 %
of the apo ensemble (43 % with rnpA). Their pipeline: 21.5 M picks over 78,777 micrographs, iterative cleaning by 2D
classification and heterogeneous refinement with decoy classes down to ~50 particles per micrograph, poses aligned on
the rigid core, 3D classification with full solvent masks and focal masks on the periphery iterated with pruning,
3DVA/3DFlex for subsets that resist discrete classes, per-conformer non-uniform refinement; conformer 18 (EMD-74209)
from 66k particles at 2.9 Å in a 480-pixel box at 0.732 Å, 351 Å across.

Our particles are the same collection binned twice in a 234 Å box, and the arm reaches 114 Å from the centre; at
1.5 µm defocus the 8 Å signal is delocalised by ~37 Å, so the tip's signal leaves our box. A local refinement of the
closed particles against the closed map made the arm worse (0.76 to 0.60), and the occupancy-layer board produced
blob artefacts; the P12/P10.1 focal boards and a defocus-split reconstruction (tip density per defocus third) are
running. The micrographs (31k at 1.464 Å) and picks are on batch2, so a re-extraction at a 240-pixel box is possible.

**Box explanation refuted (23:10).** The 40k closed particles split into defocus thirds (0.98, 1.39, 1.70 µm) give
tip density 0.60, 0.65, 0.67 of the deposited map against base density 0.75, 0.78, 0.79: the tip does not fall with
defocus, so its signal is not leaving the box and no re-extraction is warranted. The tip is uniformly weaker than the
base; the remaining candidates are a mixture at the tip, which the P12/P10.1 focal boards test, and intrinsic
flexibility, bounded by the pure-particle tip ceiling on the toy.

**The tip ceiling (23:20).** Pure simulated closed particles reconstructed at true poses reach tip density 0.64 of
the deposited map (base 0.88); at our realigned poses the same particles reach 0.46 (base 0.86). The real closed
particles reach 0.60 to 0.67 at the tip, above the pure-particle ceiling at our pose precision. The P12–P10.1 contact
is therefore lost to pose precision at the periphery, where a 3.8° error displaces the tip by 7.6 Å, and not to the
classification, the mask, the box or composition. Recovering it is the per-state refinement of the closed state at
full resolution with core-aligned poses, as the paper did with 66k particles at 2.9 Å.

**Closed-state refinement verdict (23:30).** A plain reconstruction at the poses from the closed-state local
refinement (12 to 6 Å, 240 Å mask, started from our partial-arm map) gives tip 0.26 against 0.65 at the family poses
and 7.6 Å against 5.3 Å: that refinement degraded the poses and is discarded. The per-state refinement of the closed
state has to follow the paper's order, poses aligned on the rigid core and the periphery handled by focal
classification and refinement on the purified docked set. The pure-particle tip ceiling at a 6 Å band is being
measured on the toy to match the real maps' resolution.

## 22. The docked closed state by focal classification (23:05)

A klem board whose E-step references differ only inside the P12–P10.1 region (the closed-minus-open mask of the
deposited maps) spreads the 66k closed-family particles along a docked/undocked axis that the plain board could not
see. Ranked by their responsibility-weighted focal score and reconstructed at the family poses: top 15k arm 0.82, tip
0.75, base 0.83 at 6.5 Å; next 10k tip 0.70; bottom 15k arm 0.46, tip 0.19. The 40k particles matched to conformer 18
(tip 0.65) were a docked/undocked mixture. The same board on all 410k of family 1 is running for the full docked
population; the per-state refinement of that set, core-aligned then focal, is the step that follows.

**Family 2 is not a missing-domain state (23:40).** Its mean map keeps 0.92 of the consensus's ordered voxels
against 0.97 for family 1, with no missing region larger than about 1000 Å³ (a domain would be an order of magnitude
larger) and only a small extra blob at the box centre. It reads as a lower-contrast or differently conformed
population, plausibly among the 45 apo conformers that are not deposited among our ten. The compositional question is
being asked of the toy instead (toy10c through the full recipe).

**Soft level 3 and the denoised hypothesis (23:35).** Soft families at level 3 merge at 0.71, the same as the hard
tree, with per-conformer shifts both ways: weighting the unions does not by itself repair the stranded fifth of the
particles; the full soft tree and the global re-assignment at the leaves are still running. Both the 2D expansion
project and the refinement project already read the denoised stacks, so denoising is not what separates the 2D
averages from the 3D maps; klem's own cell maps understate the arm because of their shrinkage toward the consensus
(tip 0.99 against 1.31 for reconstruct3D on the same particles), so the arm is judged on reconstruct3D maps only. The
2D route, closedness scored per sub-class average and the top decile's members reconstructed, is running.

## 23. Autonomous run toward the landscape (23:55)

Standing goal from the user: as much of the conformational landscape and all state families as the 15 % subset
allows; per-state refinement parked. In flight: toy10c through the full recipe with random and with PCA initialisation
(the compositional states 61 and 62 are the readout), the pose-marginalisation A/B, the global leaf re-assignment and
the full soft tree (the stranded-particle question), the leaf pose feedback, the family-1 focal board, the real
realignment, real level 3 on family 2, the 2D sub-class route to the contact, and a distinct-state counter
(count_states.py: surviving cells clustered by cross-half similarity) calibrated on the toy and applied to the real
leaves against the paper's 55 apo conformers.

**Counting states (00:05).** Surviving cells clustered by cross-half similarity: on the toy's level-3 leaves a 0.6
threshold gives 22 states for 10 conformers (8 represented), 0.7 gives 32 (9 represented), 0.8 gives 48 (all 10), so
the counter over-splits a discrete conformer two- to three-fold at the thresholds where it loses none. Real family 1's
leaves give 20, 27, 39 and 66 states at 0.6, 0.7, 0.8 and 0.9, which with the toy's factor is 10 to 15 distinct
states at this resolution and subset size, against the paper's 55 apo conformers from seven times the data at 2.9 Å.

**The real landscape, counted (00:15).** Across family 1's feedback board and its two grandchildren and family 2's
feedback board, 196 reproducible cells cluster into 30 states at similarity 0.6 (39 at 0.7, 56 at 0.8). Eight of the
ten deposited conformers dominate a cluster at 0.6 and nine at 0.7, with 8 to 14 reproducible states matching none of
the ten. Corrected by the toy's over-split factor, the subset resolves roughly 15 to 20 distinct states against the
paper's 55 apo conformers. The unmatched states are being reconstructed at full box as deliverable maps.

**Pose marginalisation (00:25).** Without it the children merge at 0.64, the same as with it, and the third level
at 0.69 against 0.71, for one E-step instead of seven; the whole toy tree ran in 29 minutes instead of over an hour.
It is now off by default; the real tree realignment, three hours into a marginalised parent, was cancelled and
resubmitted without it.

**Where the unmatched states live, and the 2D route (23:50).** Clustering the leaves at 0.6 gives 26 states; the
seven unmatched ones with more than 5k particles all lie in family 2 (99k particles together), so the states beyond
the deposited ten are family 2's, and they are being reconstructed. The 2D route's top decile of closedness is too
broad a cut, 179 sub-classes and 182k particles against a closed fraction near 13 %, and its map carries the contact
more weakly (0.73) than the focal docked set (1.31); a strict cut at closedness 0.6 is running.

**The unmatched states as maps (23:48).** Seven family-2 clusters of 5k to 27k particles reconstruct at 8.4 to
10.6 Å, and their deviations correlate at most +0.18 with any of the ten deposited conformers, so they are genuinely
outside the deposited set: candidates for the other 45 apo conformers, or the more flexible and partially unfolded
part of the population that the paper's pruning discarded.

**Composition on the toy, first stages (23:50).** On toy10c the children at consensus poses give 11 of 12 states a
home at mean purity 0.42: the deleted-domain state at 0.44 and the added-density state at 0.32, inside the range of
the ten conformers (0.31 to 0.55). The tree does carry compositional states as states at the second level; the pose
feedback and the third level, which took the ten-conformer toy from 0.54 to 0.71, are running on it, with a PCA-
initialised parent in parallel (0.37 against 0.26 for the random start).

**The 2D route, closed (23:51).** A strict cut at closedness 0.6 keeps 44 sub-classes and 57k particles, and their
5.2 Å map carries the contact no better than the decile cut (0.73 at the tip, 22 % above contour) and far below the
focal docked set (1.31, 66 %). Scoring averages against deposited renders tracks a global motion, not the contact;
the sub-classes are mixtures inside. For a specific contact the 3D focal board is the tool.

**Leaf pose feedback (23:58).** Refining poses against the level-3 leaves and rebuilding them merges at 0.71, the
same as before; the leaf-level refinements have too few particles per state (about 4k) to produce useful references,
and the poses move by a tenth of a degree. Pose feedback pays once, at the family level. The toy recipe ledger is now
closed: family-level poses, balance and depth are the levers; soft families, leaf feedback, sharper kernels, a
second family pass and pose marginalisation are neutral or not worth their cost.

## 24. Refining the suspected-closed states as states (00:05)

The user disagrees that the contact is a resolution matter and asks for the suspected-closed states to be refined.
Done with the tool that produced the family poses: refine3D_states with the local pose policy on the closed family,
three states from the focal ranking (docked 15k, intermediate 10k, the rest), 10 to 6 Å, wide mask; each refined state
is reconstructed unmasked and scored for the contact, and the refined docked map goes to the user's folder.

**Soft tree, final (00:02).** Soft families at level 2 merge at 0.64 and at level 3 at 0.71, identical to the hard
tree at both levels. Weighting the unions by the parent's responsibilities does not recover the stranded particles;
the global re-assignment at the leaves is the last test of that loss.

**Family 2's third level (00:05).** Without pose marginalisation the two grandchild boards of family 2 (102k and 95k
particles) ran in 18 minutes where the marginalised version would have taken hours; every cell reproduces (0.98), and
none matches a deposited conformer (means 0.17 and 0.11). Family 2 lies entirely outside the deposited ten at every
level; its leaves now enter the distinct-state count.

**The landscape, counted with family 2 (00:15).** All six leaf boards, 294 reproducible cells: 41 states at
similarity 0.6 (8 deposited conformers identified, 18 unmatched), 61 at 0.7 (9 identified, 35 unmatched), 91 at 0.8.
Corrected by the toy's over-split factor, 15 to 25 distinct states on the 15 % subset, nine of the deposited ten plus
10 to 15 unmatched states that all belong to family 2.

**Up-traversal, raw (00:08).** Each rebuilt grandchild scored every particle in the tree once, and sending each
particle to its best leaf anywhere merges at 0.73 against 0.71, but four conformers read exactly 1.00 because the
size-blind purity rewards tiny pure states. The hard and the globally re-assigned leaves are being re-scored with a
minimum state size before any verdict; the same floor is being applied to the earlier numbers for comparison.

**Up-traversal, verdict (00:20).** With states below 1000 particles excluded the hard third level reads 0.69 over
67 states and the global re-assignment 0.45 over 54: its 0.73 came from tiny pure states, and across boards the
responsibilities are not comparable, so particles collapse into over-confident cells. Together with the regrown
parent and the soft families, every way of going back up the tree fails to recover the stranded particles; that loss
is a limit of the per-particle evidence at this noise, not of the assignment rule. The hard tree stands, and purities
are quoted with a 1000-particle floor from here on.

**Composition on the toy, PCA start (00:12).** From the PCA-initialised parent the second level gives 11 of 12
states a home at mean purity 0.51, with the deleted-domain state the purest of all at 0.72 and the added-density
state, 6 % of the mass, at 0.35. The pose feedback and third level are running on both compositional trees.

**Suspected-closed states refined as states (00:25).** refine3D_states with the local policy over three states
from the focal ranking reshuffled 74 % of the particles; the docked state shrank to 11.9k and lost the contact (tip
0.40, 7.3 Å), while the states that absorbed its particles carry it below the plain focal set (0.98 to 1.07 against
1.31). The body dominates the likelihood at 6 Å and the labels reorganise around it. The core-masked variant, the
paper's order, is the last fair test.

## 26. Closed state, next tests (2026-10-16 ~01:00)
- Toy caveat: simulate_particles masks the input map at mskdiam (200 A in run_toy55_sim.sbatch), so the toy's conformer-18
  particles never contained the arm beyond 100 A radius; the "pure simulated 18 tip 0.64" ceiling was that mask. The real
  docked tip 0.75 vs deposited 1.00 is a genuine deficit.
- Family-1-wide focal board (409k): docked corner of 8 cells at +0.80..+1.11; top 3 cells 26.5k particles.
- Focal ranking shows the docked set is homogeneous across ranks (top 15k tip 0.75, next 10k 0.70, bottom 15k 0.19): the
  deficit is not a mixture signature; pose blur (toy pipeline: 4 deg median, 8 A at the tip) is the candidate.
- Running: focal rank + lever-arm local refinement (53903807), box/sub-class diagnostic (53903808), core-masked 3-state
  refinement (53903562). Focused-flex 2026-10-14 region runs: weak docked enrichment only (12 % vs 7 %).

## 27. Closed state: lever arm and the 2D sub-classes (2026-10-16 ~01:10)
- Lever-arm local refinement of the docked top 30k (refine3D refine=neigh, full mask, 10->6 A): body 5.4 -> 5.3 A, arm
  tip 0.78 -> 0.49, contact p99 1.46 -> 1.10. Third refinement in a row that blurs the arm (3-state full, 3-state core,
  1-state lever). The arm is not rigidly attached to the body pose within the docked set: heterogeneity inside the
  docked population. Next: focal board within the docked 30k (53903923) to split it by arm position.
- Sub-classes the 3D focal score calls most docked look closed to the user; their members are 29 % in the docked set
  and the remaining 71 % reconstruct open (0.67 / 0.16). A 2D average renders a coherent minority arm.
- Box/mask ruled out (232 A mask costs 4 % of the tip); pose blur of 4 deg reproduces the docked map's numbers, but
  pose refinement does not recover them, so the blur is conformational, not alignment.
- CORRECTION (01:20): the refine=neigh lever run was not local (p90 77 deg, shifts 5.7 px): verdict void. Rerun with
  refine=pose_cont (53903926); klem --pose_mask/--pose_ref prototype of the node policy on the docked 30k follows.

## 28. Closed state recovered by a second focal level (2026-10-16 ~01:30)
- Focal board INSIDE the docked top 30k, ranked, top 10k reconstructed at the family poses: arm occupancy 0.97, tip 0.84,
  contact p99 1.73 / p90 1.31 (deposited 2.36 / 1.64; family-wide top 30k 1.46 / 1.10), 7.1 A. Delivered as
  conf18_docked_within10k_sharpB150.mrc.
- Policy confirmed on the real data: at the node where the arm reveals itself, classify on it and keep the poses;
  the refine-as-states and (invalid, non-local) neigh refinements all blurred it. Local-refiner tests (refine=cont,
  klem pose step with core / region / full masks) and the wide-toy policy test are still running.
- klem pose policy on the docked 30k: family poses arm 0.76 | lever 0.77 | core (too small a mask) 0.71 | REGION 0.80.
- sub-class 1413 (user) vs 3D: our docked maps render its average at arm-region ncc 0.85-0.87 (deposited closed 0.77, open 0.47).
- USER CONFIRMED (01:40): the within-10k docked map at 10 A unsharpened / 8 A B-50 "are it". Closed state recovered on the 15 % subset.

## 29. Overnight: the truth-free node recipe (launched 2026-10-16 01:44)
- Regions of the family-1 node derived from its own feedback board: spatial components (1 big) and motion-pattern
  k-means (5 regions); one focal board per region over all 409k family-1 particles; states = clusters of reproducible
  cells, particles ranked by responsibility mass, top sets reconstructed; level 2 inside the largest states with the
  data-derived focus and the parent state's deviation as the ranking target. Deposited maps only for identification.
- Drivers 53904039 (regions 1-2) and 53904042 (regions p1-p5); results under rnasep_klem/knobfree/.
- 02:40 Fortran port started: ~/SIMPLE_klem worktree, prg=klem_estep (cross-half band-limited E-step twin of the Python), dump/compare harness on the toy; first build failed on name clashes, rebuilding.
- Re-import diagnosis: pairing and poses are exact; the remaining suspect is the stacks (family project vs denoised) - import_check2 running.
- 02:55 PORT: the Fortran klem E-step (prg=klem_estep, ~/SIMPLE_klem2) matches the Python E-step exactly (ll corr 1.0000, argmax agreement 1.000) on the toy; 20k x 51 cells in 40 s on 16 threads.
- 02:50 REGION ALIGNMENT WORKS: within-10k at region-aligned poses contact p99 1.89 / p90 1.40 / arm 1.02 (family poses 1.73 / 1.31 / 0.97; core-aligned 1.37); import path validated by an unmoved control.
- 04:00 Fortran board: load-once-per-schedule rebuilt; toy test rerun; FIRST REAL-DATA RUN = closed_focal_fortran_rnasep.sbatch (family 1, 409k, focal mask, 9 rounds) to time against the 2 h Python board and to compare the docked corner.
- 03:55 region alignment of the docked top 30k running (region_align_rec TAG=30k); within-10k region-aligned maps delivered (lp10 unsharpened, lp8 B-50).
- 04:05 wide-toy policy (truth): lever (whole reference) best poses 3.99->2.89 deg; region best arm map (0.87) but worse poses (4.59); core worst arm. Policy: lever for poses, region only for the node map.
- 04:15 region alignment on the top 30k: contact 1.46 -> 1.52, arm 0.85 -> 0.89, FSC 5.4 -> 5.9 A (gain scales with purity; the 10k set gained 1.73 -> 1.89).
- 04:25 POLICY SETTLED: whole-reference local pose step on the within-10k: contact 1.85 / arm 0.99 / FSC 6.9 A (both improve); region-masked 1.89 / 1.02 / 7.3 A; core 1.37. Matches the wide toy with truth. Classify on the region, refine poses on the whole cell reference, region-mask only for a final map polish.
- 05:20 Fortran pose step validated on the wide toy (3.81 -> 2.90 deg, Python 2.89); board v2 (devw, pnoise, latent init) built; real-data v2 focal board rerunning after a private-allocatable fix.
- 05:40 overnight truth-free boards done: 5 motion regions -> 245 reproducible cells -> 31/36/53 states at sim 0.6/0.7/0.8; identifications pending the readout fix-up.
- 05:50 ten conformers at whole-reference-refined poses: identification unchanged (10/10, 0.53), FSCs slightly better (28 6.9->6.5, 37 5.1->4.9, 46 9.0->8.1); the pose step pays on pure sets, little on broad matched states.
- 06:10 TRUTH-FREE: the data-derived moving region (region_1) identifies 8/10 deposited conformers (0.43) with no deposited mask, and holds a docked state at contact 1.41 / arm 0.79 (hand-mask 1.46 / 0.85). Level-2 ranking by the parent direction is too generic (1.23 vs 1.73); next: rank by the node's own moving region.

## 30. MORNING SUMMARY (2026-10-16 06:25)

**Closed state.** Confirmed by you from the within-10k map (10 A unsharpened / 8 A B-50). Since then the alignment half
of the node policy was settled on the real data and on the wide toy with truth: at a node, classify on the region with
the poses fixed, then refine the node's particles' poses against their OWN whole cell reference (klem's local rotation
step): contact 1.73 -> 1.85, FSC 7.1 -> 6.9 A on the docked 10k (toy: pose error 3.99 -> 2.89 deg); region-masked
alignment gives a slightly fuller arm (1.89) at the body's expense (7.3 A); core-masked alignment hurts (1.37). The
whole-reference-aligned and region-aligned maps are in Downloads/.../recovered_conformers (README has the list). On the
broad matched states (20-110k) the pose step only trims FSCs (28: 6.9 -> 6.5, 37: 5.1 -> 4.9, 46: 9.0 -> 8.1).

**Deriving the mask (your point).** The family-1 node's own reproducible-variance map, with no deposited map involved,
gives one big moving region (99 % of the variability) and five motion-pattern regions. The big region's board recovers
8 / 10 deposited conformers (mean 0.43; the hand-masked board: 10 / 10 at 0.46) and holds a docked state at contact
1.41 / arm 0.79 (hand mask: 1.46 / 0.85) -> Downloads/.../closed_state_truthfree_region1_state2_lp10_unsharpened.mrc.
Seven truth-free boards: 343 reproducible cells, 41 / 50 / 73 states at similarity 0.6 / 0.7 / 0.8 (hand-driven tree:
41 / 61 / 91). Level 2 ranked by the parent direction is too generic (1.23 vs 1.73); ranking along the node's own main
deviation axis is being tested (job 53904408).

**Fortran port.** prg=klem_estep (exact match to the Python E-step), prg=klem_board (board driver: bands, kernel,
Wiener shrinkage, focus, balance, devw, pnoise, latent init, pose step, outputs readable by the Python tools via
klem_fortran2npz.py) in ~/SIMPLE_klem2 (snapshot commit 5e1e6aeee). Toy: matches the Python board; pose step
3.81 -> 2.90 deg (Python 2.89). Real data: a 409k focal board in 41 min (Python 2 h), of which 17 min is loading; v1
(random init, no devw) under-performed the Python board; v2 (devw, pnoise, latent init): 31 min, reproducibility 0.99, the same docked cells as the Python board and an identical docked reconstruction (corr +0.42 / -0.50, 98 % density vs +0.41 / -0.51 / 102 %): a drop-in at 4x the speed.

**Traps found tonight (in memory):** simulate_particles masks at mskdiam (the toy's closed arm had been clipped);
import_particles sets state=1 on every row (selection needed after import); the denoised stacks are ctf=flip; SIMPLE's
forward FFT divides by N and centres the phase origin; refine3D neigh/cont are not local in this snapshot.

**Not done / next:** parallel particle reads in the Fortran loader; the within-node level-2 target (own axis);
the toy10c added-density state (62) is never recovered (random or PCA init); the Fortran tree driver (families,
children) is still bash + Python around the Fortran board.
- 06:35 truth-free level 2 along the node's own main axis: region_1 docked state, minus end 10k = contact 1.47 / arm 0.80 (parent direction 1.23; hand mask 1.73). Next: recurse the region derivation inside the state.
- 06:55 recursion with a region re-derived INSIDE the docked state fails (body variability dominates; best cluster 1.44): the region must be inherited from the parent node; inherited-mask test running (r1s2p).
- 07:20 FINAL: truth-free recursion (inherited region, motion-segmented) plateaus at contact 1.42-1.47 / arm 0.78-0.87 vs the hand mask 1.73 / 0.97; the node-level derivation works (8/10 conformers, docked state 1.41), the final purification still needs a contact-site-sharp mask. All jobs finished.
