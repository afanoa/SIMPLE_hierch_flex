#!/usr/bin/env python3
"""Is a group selected on the particle images' own statistics? Per image: mean, background (outer ring) mean and std,
centre-disk mean and std, power in shells 1-3 / 4-10 / 11-25 (relative to total); univariate AUC between two groups.
usage: image_stats_discrim.py STACK LABELS G1 G2 [GT]   (LABELS: one label per stack row)"""
import sys, numpy as np
stk, labf, g1, g2 = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
lab = np.loadtxt(labf).astype(int); conf = np.loadtxt(sys.argv[5], comments='#')[:, 1].astype(int) if len(sys.argv) > 5 else np.ones(len(lab), int)
h = np.frombuffer(open(stk, 'rb').read(1024), np.int32); nx, ny, nz, nsym = [int(v) for v in (h[0], h[1], h[2], h[23])]
I = np.memmap(stk, np.float32, 'r', offset=1024 + nsym, shape=(nz, ny, nx)); n = nx
rows = np.where(((lab == g1) | (lab == g2)) & (conf[:len(lab)] > 0))[0]; y = (lab[rows] == g2).astype(float)
x = np.arange(n) - n // 2; rr = np.hypot(x[:, None], x[None, :]); ring = (rr > 0.42 * n) & (rr < 0.5 * n); cen = rr < 0.25 * n
fy, fx = np.meshgrid(np.fft.fftfreq(n) * n, np.fft.fftfreq(n) * n, indexing='ij'); sh = np.round(np.hypot(fx, fy)).astype(int)
bands = {'pow 1-3': (sh >= 1) & (sh <= 3), 'pow 4-10': (sh >= 4) & (sh <= 10), 'pow 11-25': (sh >= 11) & (sh <= 25)}
F = {k: np.zeros(len(rows)) for k in ['mean', 'ring mean', 'ring std', 'centre mean', 'centre std', 'centre-ring', *bands]}
for b0 in range(0, len(rows), 2000):
    blk = np.array(I[rows[b0:b0 + 2000]], np.float64); s = slice(b0, b0 + len(blk))
    F['mean'][s] = blk.mean((1, 2)); F['ring mean'][s] = blk[:, ring].mean(1); F['ring std'][s] = blk[:, ring].std(1)
    F['centre mean'][s] = blk[:, cen].mean(1); F['centre std'][s] = blk[:, cen].std(1); F['centre-ring'][s] = F['centre mean'][s] - F['ring mean'][s]
    P = np.abs(np.fft.fft2(blk)) ** 2; tot = P.reshape(len(blk), -1).sum(1)
    for k, m in bands.items(): F[k][s] = P[:, m].sum(1) / tot
def auc(v, y):
    o = np.argsort(v); rk = np.empty(len(v)); rk[o] = np.arange(1, len(v) + 1); n1 = y.sum(); n0 = len(y) - n1
    return (rk[y == 1].sum() - n1 * (n1 + 1) / 2) / max(n0 * n1, 1)
print('%d (group %d) vs %d (group %d) images' % ((y == 0).sum(), g1, (y == 1).sum(), g2))
for k, v in sorted(F.items(), key=lambda kv: -abs(auc(kv[1], y) - 0.5)): print('  %-12s AUC %.3f   means %.5g vs %.5g' % (k, auc(v, y), v[y == 0].mean(), v[y == 1].mean()))
