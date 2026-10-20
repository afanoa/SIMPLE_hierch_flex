#!/usr/bin/env python3
"""Figure of the apo/holo x open/closed toy boards: for each board, the lattice cells coloured by their holo fraction and
by their closed fraction (cell size = particle count), plus the per-cell dominant conformer; one row per board.
usage: toyah_fig.py GT CLOSED_IDS OUT.png BOARD:LABEL [BOARD:LABEL ...]"""
import sys, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
gt_fn, cl_fn, out = sys.argv[1:4]; boards = [a.split(':', 1) for a in sys.argv[4:]]
gt = {}
for l in open(gt_fn):
    if l.strip() and not l.startswith('#'): w = l.split(); gt[int(w[0])] = int(w[1])
cl = set(int(v) for v in open(cl_fn).read().strip().split(','))
fig, axes = plt.subplots(len(boards), 3, figsize=(13, 3.6 * len(boards)), squeeze=False)
for row, (bd, label) in enumerate(boards):
    L = np.loadtxt(bd + '/klem_latent.txt', comments='#'); pind = L[:, 0].astype(int); am = L[:, 3].astype(int) - 1; gm = L[:, 5]
    conf = np.array([gt.get(p, 0) for p in pind]); on = (gm > 0.5) & (conf > 0); holo = conf >= 101; closed = np.array([c in cl for c in conf])
    # grid geometry from the resp file (i1, i2 per cell)
    with open(bd + '/klem_resp.bin', 'rb') as f:
        n, ncl, ngrid, n1, n2 = np.fromfile(f, np.int32, 5); _ = np.fromfile(f, np.int32, n); i1 = np.fromfile(f, np.int32, ngrid); i2 = np.fromfile(f, np.int32, ngrid)
    cnt = np.array([(on & (am == c)).sum() for c in range(ngrid)]); fh = np.array([holo[on & (am == c)].mean() if cnt[c] else np.nan for c in range(ngrid)]); fc = np.array([closed[on & (am == c)].mean() if cnt[c] else np.nan for c in range(ngrid)])
    dom = []
    for c in range(ngrid):
        m = on & (am == c); dom.append(np.bincount(conf[m]).argmax() if m.any() else 0)
    for col, (vals, title, cmap, vmin, vmax) in enumerate(((fh, 'holo fraction per cell', 'coolwarm', 0, 1), (fc, 'closed fraction per cell', 'viridis', 0, max(0.5, np.nanmax(fc) if np.isfinite(np.nanmax(fc)) else 0.5)), (np.array(dom, float), 'dominant conformer per cell', 'tab20', 0, 110))):
        ax = axes[row, col]; H = np.full((n1, n2), np.nan)
        for c in range(ngrid): H[i1[c], i2[c]] = vals[c]
        im = ax.imshow(H, cmap=cmap, vmin=vmin, vmax=vmax, origin='lower')
        for c in range(ngrid):
            if cnt[c] > 0: ax.text(i2[c], i1[c], ('%d' % dom[c]) if col == 2 else ('%.2f\n%dk' % (vals[c], cnt[c] // 1000)), ha='center', va='center', fontsize=6, color='k' if col != 2 else 'w')
        ax.set_title('%s: %s' % (label, title), fontsize=9); ax.set_xticks([]); ax.set_yticks([])
        if col < 2: plt.colorbar(im, ax=ax, fraction=0.046)
    purity_h = np.nanmean(fh[fh > 0.5]) if (fh > 0.5).any() else float('nan'); cover_h = holo[on & np.isin(am, np.where(fh > 0.5)[0])].sum() / max(holo[on].sum(), 1)
    axes[row, 0].set_xlabel('holo-majority cells: %d, hold %.2f of holo at purity %.2f' % ((fh > 0.5).sum(), cover_h, purity_h), fontsize=8)
fig.suptitle('apo/holo x open/closed toy: 10 apo + 8 holo deposited conformers (closed: 18, 51, 102, 105, 106), 160k conformer particles, real poses/CTF/SNR', fontsize=10)
plt.tight_layout(); plt.savefig(out, dpi=130); print('->', out)
