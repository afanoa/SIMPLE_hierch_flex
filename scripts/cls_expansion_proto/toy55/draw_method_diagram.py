#!/usr/bin/env python3
"""Method diagram of klem as a PNG (matplotlib, wrapped text): the board (left), the tree (middle), safeguards (right)."""
import textwrap, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt; from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
fig = plt.figure(figsize=(24, 14), dpi=110); ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 24); ax.set_ylim(0, 14); ax.axis('off')
def box(x, y, w, h, title, body='', fc='#eef3fb', ec='#2b4c7e', width=78, fs=10.5):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0.02,rounding_size=0.15', fc=fc, ec=ec, lw=1.5))
    ax.text(x + w / 2, y + h - 0.28, title, ha='center', va='top', fontsize=12.5, fontweight='bold')
    if body: ax.text(x + 0.18, y + h - 0.62, '\n'.join(textwrap.fill(p, width) for p in body.split('\n')), ha='left', va='top', fontsize=fs, color='#222', linespacing=1.35)
def arrow(x, y1, y2, label=''):
    ax.add_patch(FancyArrowPatch((x, y1), (x, y2), arrowstyle='-|>', mutation_scale=18, color='#2b4c7e', lw=1.6))
    if label: ax.text(x + 0.12, (y1 + y2) / 2, label, fontsize=9.5, color='#2b4c7e', va='center')
fig.text(0.5, 0.975, 'klem: hierarchical kernel-latent EM.  One node = one board.  Evidence at every step = agreement of the two independent halves.', ha='center', fontsize=16, fontweight='bold')
# left: the board
box(0.4, 12.0, 7.2, 1.15, 'INPUT', 'particle images at the consensus poses, per-image CTF (|CTF| for phase-flipped stacks), per-particle contrast scale; flex_pca latent only as the initial placement on the lattice', fc='#f7f7f7', ec='#888')
box(0.4, 5.3, 7.2, 6.4, 'ONE BOARD = one node', 
'7 x 7 lattice of free Fourier-cell volumes V_k plus 2 free junk components. No basis, no deformation: every cell is an unconstrained reconstruction.\n'
'E-step (per image, pose + CTF aware): responsibility r_ik of image i for cell k from the likelihood of the image against the OTHER half\'s cell projection; smoothed by the lattice kernel (neighbouring cells share particles), Sinkhorn balance, per-particle noise factor. Softmax posteriors: no hard labels, particles never leave.\n'
'M-step: per-cell Wiener reconstruction with the responsibilities as weights, shrinkage (2 %) toward the mean cell.\n'
'Bands 20 -> 14 -> 10 (-> 8) A, 3 iterations each.\n'
'Output per cell: two half-volumes -> deviations D0, D1 = cell minus mean.\n'
'EVIDENCE: a cell survives only if corr(D0, D1) >= 0.5 on the band.', width=70)
box(0.4, 4.0, 7.2, 1.05, 'distributed board', 'master + N persistent workers (10-20 x 8 threads), half-sphere Fourier cubes exchanged through files in single precision; ~5 min for 692k particles at 20-10 A', fc='#f7f7f7', ec='#888')
box(0.4, 1.0, 7.2, 2.7, 'readouts (never inside the fit)', 'identification of cells/states against deposited maps (correlation of deviations at 60-20 A); the whole-tree embedding (UMAP of the concatenated responsibilities of every board a particle visited); per-island truth on toys', fc='#f7f7f7', ec='#888')
# middle: the tree
X, W = 8.4, 7.6; ys = [12.1, 10.75, 9.1, 7.75, 6.3, 4.6, 3.25, 1.75]
box(X, ys[0], W, 1.0, 'PARENT BOARD', 'all particles: coarse ordering and junk removal')
arrow(X + W / 2, ys[0], ys[1] + 1.0)
box(X, ys[1], W, 1.0, 'FAMILIES', 'spectral bisection of the surviving-cell deviation graph (2 families)')
arrow(X + W / 2, ys[1], ys[2] + 1.3)
box(X, ys[2], W, 1.3, 'FAMILY NODE', 'child board -> k-means seeds (10, never reported) -> refine3D_states (local poses per seed, 20 -> 8 A) -> FEEDBACK BOARD at the refined poses')
arrow(X + W / 2, ys[2], ys[3] + 1.0)
box(X, ys[3], W, 1.0, 'REGIONS (data-derived)', 'cross-half product of the deviations -> spatial components + motion-pattern groups, written as soft masks')
arrow(X + W / 2, ys[3], ys[4] + 1.1)
box(X, ys[4], W, 1.1, 'REGION BOARD (focus mask = the region)', 'states = clusters of surviving cells (average linkage at cross-half similarity 0.7)')
arrow(X + W / 2, ys[4], ys[5] + 1.35)
box(X, ys[5], W, 1.35, 'SUPPORT TEST (no particle counts anywhere)', 'a state exists iff its mean deviation reproduces between the halves (>= 0.7) AND is separable from every other state (< 0.7); an unsupported cluster merges into its nearest neighbour only if the merged cluster passes, otherwise no state')
arrow(X + W / 2, ys[5], ys[6] + 1.0)
box(X, ys[6], W, 1.0, 'LEVEL 2 (largest supported states, inherited region)', 'board on the state\'s particles, ranked along the node chart (cell PPCA); same support test')
arrow(X + W / 2, ys[6], ys[7] + 1.15)
box(X, ys[7], W, 1.15, 'LANDSCAPE', 'leaves deduplicated by cross-half similarity: one count per threshold 0.6 / 0.7 / 0.8 / 0.9; states reconstructed at the refined poses')
# right: safeguards
box(16.6, 9.0, 7.0, 4.15, 'WHAT KEEPS IT FROM OVER-SPLITTING',
'1. Every split is made from the two independent halves and kept only where they agree (cells 0.5, states 0.7, regions from the cross-half product).\n'
'2. The lattice kernel keeps a continuum as neighbouring cells; it is not cut into pseudo-states.\n'
'3. Separability: two clusters more alike than 0.7 across halves are one state.\n'
'4. Junk components absorb outliers instead of seeding states.\n'
'5. The landscape is reported at four thresholds and deduplicated across leaves.\n'
'6. Level 2 runs only below supported states and stops when its own support fails or the band is exhausted.\n'
'7. The k-means seeds are seeds for refinement, never reported states.', fc='#fff4e5', ec='#b36b00', width=66)
box(16.6, 6.5, 7.0, 2.3, 'NOT YET TESTED FOR SUPPORT', 'the family split itself (always 2 families); the number of regions per family (capped by maxreg); the band at which a node is called exhausted (8 A on this data)', fc='#fdecea', ec='#a12', width=66)
box(16.6, 3.6, 7.0, 2.7, 'POSES AND REFINEMENT', 'consensus poses at the parent; per-seed local refinement at the family node; inherited below. Toy (poses already within 1.4 deg): every local refinement against mixture references made them worse (2.8 deg) -> refine late, against pure leaf states only. Real data: A/B running (f692c with, f692d without).', fc='#eaf6ea', ec='#2a7', width=66)
box(16.6, 1.0, 7.0, 2.4, 'WHERE flex_pca FITS (per-image model)', 'now: initial placement + the node chart from the cells (PPCA on cell volumes, BIC rank). Proposed: mixture of PPCAs on the lattice: each node a flex_pca basis fitted to the per-image residuals; responsibilities from the per-image evidence marginalised over that basis; per-image coordinates within the node.', fc='#eaf6ea', ec='#2a7', width=66)
fig.savefig('/mnt/beegfs/elmlund/afan/routes/toy55/figs/klem_method_diagram.png', dpi=110, facecolor='white'); print('written')
