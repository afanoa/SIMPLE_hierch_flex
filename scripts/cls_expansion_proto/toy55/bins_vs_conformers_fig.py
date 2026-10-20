#!/usr/bin/env python3
"""Figure: klem latent bins (full-box reconstructions) against the docked deposited conformers. Left: correlation of
(bin - mean of bins) with (conformer - consensus) in the 60-20 A band. Right: projections of the band-limited
difference maps for the end bins and the matching conformer differences. usage: bins_vs_conformers_fig.py BINS_DIR CONF_DIR CONSENSUS OUT.png"""
import sys, glob, os, re, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
bins_dir, conf_dir, cons_f, out = sys.argv[1:5]; SMPD = 1.464
def read_mrc(p):
    h = np.fromfile(p, dtype=np.int32, count=256); nx, ny, nz = h[0], h[1], h[2]; return np.fromfile(p, dtype=np.float32, offset=1024 + h[23]).reshape(nz, ny, nx).astype(float)
bf = sorted(glob.glob(bins_dir + '/bin*.mrc'), key=lambda f: int(re.findall(r'bin(\d+)', f)[0])); B = [read_mrc(f) for f in bf]; n = B[0].shape[0]
cf = sorted(glob.glob(conf_dir + '/conf*_m160.mrc')); cid = [int(re.findall(r'conf(\d+)', os.path.basename(f))[0]) for f in cf]; C = [read_mrc(f) for f in cf]; cons = read_mrc(cons_f)
zz, yy, xx = np.indices((n, n, n)); r = np.hypot(np.hypot(xx - n / 2, yy - n / 2), zz - n / 2); msk = 1. / (1. + np.exp((r - 100 / SMPD) / 2.))
F = lambda v: np.fft.fftshift(np.fft.fftn(np.fft.ifftshift(v * msk))); band = (r >= n * SMPD / 60.) & (r <= n * SMPD / 20.)
Bm = np.mean(B, 0); DB = [F(b - Bm) for b in B]; DC = [F(c - cons) for c in C]
def bc(x, y): return float(np.real(np.sum(np.conj(x[band]) * y[band])) / np.sqrt(np.sum(np.abs(x[band]) ** 2) * np.sum(np.abs(y[band]) ** 2)))
R = np.array([[bc(d, t) for t in DC] for d in DB]); print('corr matrix (bins x conformers):'); print(np.array2string(R, precision=2))
lp = lambda X: np.real(np.fft.fftshift(np.fft.ifftn(np.fft.ifftshift(X * (r <= n * SMPD / 20.)))))
fig = plt.figure(figsize=(17, 8)); gs = fig.add_gridspec(2, 6, width_ratios=[2.2, 1, 1, 1, 1, 1])
ax = fig.add_subplot(gs[:, 0]); im = ax.imshow(R, cmap='RdBu_r', vmin=-0.4, vmax=0.4, aspect='auto'); ax.set_xticks(range(len(cid))); ax.set_xticklabels(cid); ax.set_yticks(range(len(B))); ax.set_yticklabels(['bin %d' % (i + 1) for i in range(len(B))])
ax.set_xlabel('deposited conformer'); ax.set_title('corr(bin - mean, conformer - consensus), 60-20 A'); plt.colorbar(im, ax=ax, fraction=0.04)
for i in range(len(B)):
    for j in range(len(cid)): ax.text(j, i, '%+.2f' % R[i, j], ha='center', va='center', fontsize=7)
picks = [0, len(B) // 2, len(B) - 1]; best = [cid[int(np.argmax(R[p]))] for p in picks]
for k, p in enumerate(picks):
    a = fig.add_subplot(gs[0, 1 + k]); d = lp(DB[p]); a.imshow(d.sum(0), cmap='RdBu_r', vmin=-np.abs(d.sum(0)).max(), vmax=np.abs(d.sum(0)).max()); a.set_title('bin %d - mean (proj z)' % (p + 1)); a.set_xticks([]); a.set_yticks([])
    a = fig.add_subplot(gs[1, 1 + k]); d = lp(DC[cid.index(best[k])]); a.imshow(d.sum(0), cmap='RdBu_r', vmin=-np.abs(d.sum(0)).max(), vmax=np.abs(d.sum(0)).max()); a.set_title('conformer %d - consensus' % best[k]); a.set_xticks([]); a.set_yticks([])
a = fig.add_subplot(gs[0, 4]); a.imshow(lp(F(B[0])).sum(0), cmap='gray'); a.set_title('bin 1 map (proj z, 20 A)'); a.set_xticks([]); a.set_yticks([])
a = fig.add_subplot(gs[1, 4]); a.imshow(lp(F(B[-1])).sum(0), cmap='gray'); a.set_title('bin %d map' % len(B)); a.set_xticks([]); a.set_yticks([])
a = fig.add_subplot(gs[0, 5]); a.imshow(lp(F(C[cid.index(best[0])])).sum(0), cmap='gray'); a.set_title('conformer %d map' % best[0]); a.set_xticks([]); a.set_yticks([])
a = fig.add_subplot(gs[1, 5]); a.imshow(lp(F(C[cid.index(best[-1])])).sum(0), cmap='gray'); a.set_title('conformer %d map' % best[-1]); a.set_xticks([]); a.set_yticks([])
plt.tight_layout(); plt.savefig(out, dpi=110); print('figure ->', out)
