---
name: rnasep-closed-state-not-pure-in-klem-leaves
description: arm readout (closed-arm mask occupancy vs the deposited closed map, conf 18 = 1.00, open conformers 0.11-0.34): no f692c leaf or consolidated state is pure closed (max 0.71; the docking "18 identified" was loose); mixtures of 60-70 % closed; purification pass under test
metadata:
  type: project
---

2026-10-17: arm_readout.py (closed-arm mask = where the deposited closed conformer 18 has density the consensus lacks,
10 A, standardised maps; occupancy relative to the deposited closed map). Calibration: conf18 1.00, open conformers 0.11-0.34
(conf10/51 0.45, partially closed), consensus 0.29. f692c: best leaf maps 0.63-0.71 (family 1 motion regions, level 1 and 2),
consolidated state 4 (56k particles, docking said 18 at corr 0.46) 0.62; 0 of 117 leaf maps and 0 of 12 consolidated states
above 0.71 -> states that are ~60-70 % closed, i.e. mixtures, no pure closed state; 26/117 intermediates, no isolated arm
position. The user saw this by eye ("I only see one or two closed", "many are blurry") and the docking correlation hid it.
Suspected cause: Sinkhorn balance + lattice kernel spread a 13 % minority over neighbouring cells (toy leaf purity was only
0.21-0.29). Test running: purify.sbatch = consolidated states as starting labels -> refine3D_states re-assignment (prob_neigh,
local, 10->8 A, maxits 5) on the 150k subset -> arm readout of the output maps (job 53917716).

Purification result (14:52): refine3D_states re-assignment of the 150k subset among the 12 consolidated states moved 79k of 126k particles but the best closed-arm occupancy stayed 0.63 -> a whole-map likelihood does not separate the arm. Next test: arm_board.sbatch = klem board on family 1 with focus = the arm's travel mask (make_arm_mask.py), 3x3 lattice, balance 0 vs 0.5, closed_arm_cells.py per cell (job 53918079).

Arm-focused boards (15:05-15:09): 3x3 lattice on the arm travel mask, balance 0 or 0.5 (dev family 1, 68k): max cell 0.51 on the cell scale (consensus 0.53); seeded 2-cell arm EM on the full family 1 (335k, 20k closed-ish seeds vs rest, no balance): converges to 185k/150k, max cell 0.65. Verdict: the closed population exists (best states 2-3x above any open conformer in the arm) but at this SNR the arm (1.7 % of the box) gives a per-image likelihood difference at the noise floor, so every per-particle assignment splits ~randomly and every pooled map is a 50-65 % mixture. Purification needs pooled evidence before the decision (2D-expansion averages) or more closed particles.

RESOLVED 17:30: the closed map of 2026-10-16 came from the docked-node route (f692_docked.sh: cells matched to conformer 18 -> top ~10k particles -> zoom board -> small states); today's rerun on f692c family 1 reproduces it: docked states 02 and 05 score 0.94 and 0.95 closed-arm occupancy with 4.9k and 1.5k particles (clusters of 9 and 4 cells, mass 6.2k / 2.6k). So the closed population is ~2 % of family 1 (not 13 %), fills <= 1 lattice cell, and the truth-free leaves fold it into open groups because leaf clustering and the merge tree use whole-map similarity and 30k/union selections. Fix (truth-free): zoom on every reproducible cluster regardless of mass, region-local merge similarity, state maps from the confident core (top particles ~ cluster mass).

Sharp final kernel test (motion_1 board, sigma 0.7 -> 0.3 -> 0.15): max cell closed occupancy 0.66 -> 0.70 -> 0.71 (cell scale): marginal; the kernel is not the main dilution, the re-partition of a subset (zoom) is. Zoom on every leaf state implemented (zoom_one.sh / zoom_leaves.sbatch: 7x7 board + 1 free comp on each level-1 state's particles, region focus, supported sub-states reconstructed; outputs fam?/zoomregion<reg>_s<NN>/level1; merge/gallery/readout scripts take LEAFGLOB='zoomregion*'; zoom_analysis.sbatch), launched on f692c 17:38 (jobs 53918402/3, 47 + 50 leaves, 4 zooms at a time per 48-core node).

ZOOM RESULT 18:24 (truth-free): zoom of f692c fam1 regionmotion_1 state 5 -> 4 sub-states (6.1k/5.3k/4.8k/4.1k), sub-state 3 = 0.92 closed-arm occupancy with 4.8k particles (docked route 0.94); neighbouring leaves' zooms give 0.75-0.76 sub-states. The closed state is recoverable without the deposited map: leaf -> zoom board -> supported sub-states.

Seeded diagnostic (closed_seed_refine.sbatch, 18:38): refine3D nstates=2 on the 150k subset seeded with vol1 = consensus, vol2 = deposited closed map (full-SNR references), prob, 10->8 A, 6 iterations: the closed state fills with 14.5k particles (9.7 %) but its map scores only 0.66 closed-arm occupancy (state 1: 0.15, pure open). Even template-seeded per-image classification gives a mixture; the pooled zoom (0.92) beats it.

**Why:** "closed identified" by docked-map correlation is not evidence of a pure state; the arm occupancy is.
**How to apply:** report closed/open/intermediate occupancies for every delivered state; a purification E-step after
consolidation is probably a required last step of the method.

**Never-refined quarter tested 2026-10-18 (unrefined_closed.sbatch, unrefined_board.sbatch):** the two-reference
refine3D (consensus vs deposited closed, refine=prob) drains 90 % into the closed-seeded state in BOTH the 240k
never-refined rows and a 240k refined control (0.66 occupancy) - reference sharpness, not closed content; useless as a
test. The chain's own machinery (single-reference global poses, then the arm-focused free-cell board) gives ONE
supported state on the 240k (166k ptcls, closed occupancy 0.23 = open consensus): the missing closed particles are NOT
in the never-refined quarter; the 692k set's closed shortfall (5-16k vs ~90k expected at the paper's 13 %) is in the
picking/cleaning upstream. Closed in f692e: zoom state_01 4.9k at 0.94 (f692e_closed_substates in the user's folder).
