#!/usr/bin/env python3
"""Montage of maps as three orthogonal projections at a common low-pass: OUR reconstructions (klem latent bins) in one
block, the deposited conformers in another, same frame and contrast scaling. usage: montage_maps.py OUT.png --lp 10 map1.mrc:label map2.mrc:label ..."""
import sys, argparse, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
ap = argparse.ArgumentParser(); ap.add_argument('out'); ap.add_argument('--lp', type=float, default=10.); ap.add_argument('--smpd', type=float, default=1.464); ap.add_argument('maps', nargs='+'); a = ap.parse_args()
def read_mrc(p):
    h = np.fromfile(p, dtype=np.int32, count=256); nx, ny, nz = h[0], h[1], h[2]; return np.fromfile(p, dtype=np.float32, offset=1024 + h[23]).reshape(nz, ny, nx).astype(float)
items = [(m.split(':')[0], m.split(':')[1] if ':' in m else m) for m in a.maps]; V = [read_mrc(f) for f, _ in items]; n = V[0].shape[0]
zz, yy, xx = np.indices((n, n, n)); r = np.hypot(np.hypot(xx - n / 2, yy - n / 2), zz - n / 2); msk = 1. / (1. + np.exp((r - 100 / a.smpd) / 2.)); lpm = r <= n * a.smpd / a.lp
def lp(v): return np.real(np.fft.fftshift(np.fft.ifftn(np.fft.ifftshift(np.fft.fftshift(np.fft.fftn(np.fft.ifftshift(v * msk))) * lpm))))
P = [[lp(v).sum(ax) for ax in (0, 1, 2)] for v in V]
fig, axes = plt.subplots(3, len(V), figsize=(2.6 * len(V), 7.8), squeeze=False)
for j, (f, lab) in enumerate(items):
    for i in range(3):
        img = P[j][i]; vmax = np.percentile(img, 99.7); axes[i, j].imshow(img, cmap='gray', vmin=0, vmax=vmax); axes[i, j].set_xticks([]); axes[i, j].set_yticks([])
        if i == 0: axes[i, j].set_title(lab, fontsize=10)
    axes[0, 0].set_ylabel('proj z'); axes[1, 0].set_ylabel('proj y'); axes[2, 0].set_ylabel('proj x')
plt.suptitle('projections at %g A low-pass, same frame and contrast scaling' % a.lp); plt.tight_layout(); plt.savefig(a.out, dpi=110); print('figure ->', a.out)
