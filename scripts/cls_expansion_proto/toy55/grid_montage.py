#!/usr/bin/env python3
"""The 7x7 grid as what it is: a chessboard of maps. Each panel = one cell's map projected along z (left figure) and
its deviation from the mean cell (right figure), placed at the cell's grid position; particle counts per cell written in.
usage: grid_montage.py CELLS_DIR LATENT OUT.png [--n1 7 --n2 7]"""
import sys, glob, re, argparse, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
ap = argparse.ArgumentParser(); ap.add_argument('cells'); ap.add_argument('latent'); ap.add_argument('out'); ap.add_argument('--n1', type=int, default=7); ap.add_argument('--n2', type=int, default=7); a = ap.parse_args()
def read_mrc(p):
    h = np.fromfile(p, dtype=np.int32, count=256); nx, ny, nz = h[0], h[1], h[2]; return np.fromfile(p, dtype=np.float32, offset=1024 + h[23]).reshape(nz, ny, nx).astype(float)
V = {}
for f in glob.glob(a.cells + '/cell_??_??.mrc'):
    i, j = map(int, re.findall(r'cell_(\d+)_(\d+)', f)[0]); V[(i, j)] = read_mrc(f)
mean = np.mean(list(V.values()), 0); L = np.loadtxt(a.latent, comments='#'); am = L[:, 3].astype(int) - 1; gm = L[:, 5]; cnt = np.bincount(am[gm > 0.5], minlength=a.n1 * a.n2)
P = {k: v.sum(0) for k, v in V.items()}; D = {k: (v - mean).sum(0) for k, v in V.items()}; vmax = max(np.abs(d).max() for d in D.values()); pmax = max(p.max() for p in P.values())
fig, axes = plt.subplots(a.n2, 2 * a.n1 + 1, figsize=(2.0 * (2 * a.n1 + 1), 2.0 * a.n2), gridspec_kw=dict(width_ratios=[1] * a.n1 + [0.25] + [1] * a.n1))
for i in range(a.n1):
    for j in range(a.n2):
        r = a.n2 - 1 - j; ax = axes[r, i]; ax.imshow(P[(i, j)], cmap='gray', vmin=0, vmax=pmax); ax.set_xticks([]); ax.set_yticks([]); ax.text(2, 6, '%dk' % (cnt[i * a.n2 + j] // 1000), color='yellow', fontsize=8, va='top')
        ax2 = axes[r, a.n1 + 1 + i]; ax2.imshow(D[(i, j)], cmap='RdBu_r', vmin=-vmax, vmax=vmax); ax2.set_xticks([]); ax2.set_yticks([])
        if i == 0: ax.set_ylabel('z2 = %d' % j, fontsize=9)
        if r == a.n2 - 1: ax.set_xlabel('z1 = %d' % i, fontsize=9); ax2.set_xlabel('z1 = %d' % i, fontsize=9)
for r in range(a.n2): axes[r, a.n1].axis('off')
axes[0, a.n1 // 2].set_title('cell maps (projection along z); yellow = particles per cell', fontsize=11); axes[0, a.n1 + 1 + a.n1 // 2].set_title('cell minus the mean of all cells (red = extra density, blue = missing)', fontsize=11)
plt.tight_layout(); plt.savefig(a.out, dpi=90); print('figure ->', a.out)
