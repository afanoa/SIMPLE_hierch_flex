#!/usr/bin/env python3
"""Shape of a coordinate stage's leaf path (leaves in coordinate order, independent half-maps): a mixture axis between discrete structures
is a straight line in map space, a motion (density moving) bends. Cross-half covariance of the leaves' centred maps over the 40-8 A band
(unbiased: halves' noise independent) -> eigenvalue shares; straightness = cross-half end-to-end distance / sum of consecutive distances.
usage: path_shape.py RECONSTRUCT3D_DIR [--nb 10 --mskdiam 220 --hp 40 --lp 8]"""
import sys, argparse, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('rd'); ap.add_argument('--nb', type=int, default=10); ap.add_argument('--mskdiam', type=float, default=220.)
ap.add_argument('--hp', type=float, default=40.); ap.add_argument('--lp', type=float, default=8.); a = ap.parse_args()
def rd(p):
    h = np.fromfile(p, np.int32, 256); n = int(h[0]); return np.fromfile(p, np.float32, offset=1024 + int(h[23]), count=n ** 3).reshape(n, n, n).astype(np.float64)
f0 = '%s/recvol_state01_even.mrc' % a.rd; n = int(np.fromfile(f0, np.int32, 1)[0]); smpd = float(np.fromfile(f0, np.float32, 13)[10]) / n
g = np.arange(n) - n // 2; R = np.sqrt(g[:, None, None] ** 2 + g[None, :, None] ** 2 + g[None, None, :] ** 2); r0 = a.mskdiam / 2 / smpd
M = np.where(R <= r0, 1., np.where(R <= r0 + 6, 0.5 + 0.5 * np.cos(np.pi * (R - r0) / 6), 0.))
f = np.fft.fftfreq(n); kx, ky, kz = np.meshgrid(f, f, np.fft.rfftfreq(n), indexing='ij'); res = np.where(np.sqrt(kx ** 2 + ky ** 2 + kz ** 2) > 0, smpd / np.maximum(np.sqrt(kx ** 2 + ky ** 2 + kz ** 2), 1e-9), 1e9)
band = (res <= a.hp) & (res >= a.lp); w = np.where((kz == 0) | (np.abs(kz) == 0.5), 1., 2.)[band]
E = np.array([np.fft.rfftn(M * rd('%s/recvol_state%02d_even.mrc' % (a.rd, s)))[band] for s in range(1, a.nb + 1)]); O = np.array([np.fft.rfftn(M * rd('%s/recvol_state%02d_odd.mrc' % (a.rd, s)))[band] for s in range(1, a.nb + 1)])
E -= E.mean(0); O -= O.mean(0); C = np.real((E * w) @ O.conj().T); C = (C + C.T) / 2; ev = np.sort(np.linalg.eigvalsh(C))[::-1]; pos = ev[ev > 0]
d = lambda i, j: max(C[i, i] + C[j, j] - 2 * C[i, j], 0.) ** 0.5; st = d(0, a.nb - 1) / max(sum(d(i, i + 1) for i in range(a.nb - 1)), 1e-30)
print('%s: eigen shares %s | straightness %.2f (1 = straight line)' % (a.rd.split('rnasep_klem/')[-1], ' '.join('%.2f' % v for v in pos[:4] / pos.sum()), st))
