#!/usr/bin/env python3
"""Render a dumped pair difference (merge_cross_tree.py DUMPPAIR): central slices of the tree-A and tree-B differences
(raw and after the nuisance projection, 15 A) and the per-shell cross-tree FSC. usage: pairdump_plot.py NPZ OUT_PNG"""
import sys, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
z = np.load(sys.argv[1]); c = z['dA'].shape[0]; f = np.fft.fftfreq(c); kx, ky, kz = np.meshgrid(f, f, f, indexing='ij'); G = np.exp(-0.5 * (kx**2 + ky**2 + kz**2) * (15. / (1.464 * 2) / 1.18) ** 2)
lp = lambda v: np.real(np.fft.ifftn(np.fft.fftn(v) * G))
fig = plt.figure(figsize=(12, 9)); rows = [('raw A', z['rawA']), ('raw B', z['rawB']), ('projected A', z['dA']), ('projected B', z['dB'])]
for r, (nm, v) in enumerate(rows):
    w = lp(v); m = np.abs(w).max(); h = c // 2
    for j, sl in enumerate((w[h], w[:, h], w[:, :, h])):
        ax = fig.add_subplot(5, 3, 3 * r + j + 1); ax.imshow(sl, cmap='RdBu_r', vmin=-m, vmax=m); ax.set_xticks([]); ax.set_yticks([])
        if j == 0: ax.set_ylabel(nm, fontsize=8)
ax = fig.add_subplot(5, 1, 5); k = np.arange(len(z['fsc'])); ax.plot(k[2:41], z['fscraw'][2:41], label='raw'); ax.plot(k[2:41], z['fsc'][2:41], label='projected'); ax.axhline(0, color='k', lw=0.5)
ax.set_xlabel('shell (box 160 at 1.464 A: 234/k A)'); ax.set_ylabel('cross-tree FSC'); ax.legend(fontsize=8); ax.set_title('pair statistic %.2f' % float(z['stat']), fontsize=9)
plt.tight_layout(); plt.savefig(sys.argv[2], dpi=80); print('png', sys.argv[2])
