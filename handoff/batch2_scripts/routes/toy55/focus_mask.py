#!/usr/bin/env python3
"""Focus mask of one tree node from its members' half-maps: where the members differ between independent halves. Every member's
even / odd map is regressed per Fourier shell on the node's weighted mean (cross-half coefficient; removes scale, resolution and
blur differences), low-passed to LPMASK; v(x) = sum_j w_j dE_j(x) dO_j(x) / sum_j w_j (the node's cross-half variance = its
reproducible principal components summed) is zero in expectation where the members agree and symmetric about zero there, so the
voxels above the (1 - ALPHA) quantile of -v are significant; the mask keeps the significant voxels holding MASS of the significant
variance, dilated 3 px, soft 2 px, inside the sphere.
usage: focus_mask.py OUT.mrc --list members.txt (lines: name even.mrc odd.mrc particles) [--mskdiam 220 --lpmask 20 --alpha 0.001 --mass 0.8]"""
import argparse, warnings, numpy as np
warnings.filterwarnings('ignore')
ap = argparse.ArgumentParser(); ap.add_argument('out'); ap.add_argument('--list', required=True); ap.add_argument('--mskdiam', type=float, default=220.)
ap.add_argument('--lpmask', type=float, default=20.); ap.add_argument('--alpha', type=float, default=0.001); ap.add_argument('--mass', type=float, default=0.8)
a = ap.parse_args()
def rd(p):
    h = np.fromfile(p, np.int32, 256); n = int(h[0]); return np.fromfile(p, np.float32, offset=1024 + int(h[23]), count=n ** 3).reshape(n, n, n).astype(np.float64)
def wr(p, x):
    h = np.zeros(256, np.int32); h[0:3] = x.shape; h[3] = 2; h[7:10] = x.shape; hf = h.view(np.float32); hf[10:13] = np.array(x.shape) * smpd; hf[13:16] = 90.
    h[16:19] = [1, 2, 3]; hf[19], hf[20], hf[21] = x.min(), x.max(), x.mean(); h[52] = int.from_bytes(b'MAP ', 'little'); h[53] = 0x4444
    with open(p, 'wb') as f: f.write(h.tobytes()); f.write(x.astype(np.float32).tobytes())
L = [l.split() for l in open(a.list) if l.strip() and not l.startswith('#')]; cnt = np.array([float(l[3]) for l in L]); w = cnt / cnt.sum()
n = int(np.fromfile(L[0][1], np.int32, 1)[0]); smpd = float(np.fromfile(L[0][1], np.float32, 13)[10]) / n
g = np.arange(n) - n // 2; R = np.sqrt(g[:, None, None] ** 2 + g[None, :, None] ** 2 + g[None, None, :] ** 2); r0 = a.mskdiam / 2 / smpd
M = np.where(R <= r0, 1., np.where(R <= r0 + 6, 0.5 + 0.5 * np.cos(np.pi * (R - r0) / 6), 0.)); inm = M > 0.5
f = np.fft.fftfreq(n); kx, ky, kz = np.meshgrid(f, f, np.fft.rfftfreq(n), indexing='ij'); kr = np.sqrt(kx ** 2 + ky ** 2 + kz ** 2)
band = kr / smpd <= 1.5 / a.lpmask; shb = np.round(kr * n).astype(int)[band]; nb = shb.max() + 1; wb = np.where((kz == 0) | (np.abs(kz) == 0.5), 1., 2.)[band]
lpb = np.exp(-(kr / smpd * a.lpmask) ** 2 * np.log(2))[band]
EF = np.array([np.fft.rfftn(M * rd(l[1]))[band] for l in L]); OF = np.array([np.fft.rfftn(M * rd(l[2]))[band] for l in L])
def ssum(z): return np.bincount(shb, weights=wb * z, minlength=nb)
def back(Fb): x = np.zeros(kr.shape, complex); x[band] = Fb * lpb; return np.fft.irfftn(x, s=(n, n, n))[inm]
EuF, OuF = w @ EF, w @ OF; den = ssum(np.real(EuF * OuF.conj())); okd = den > 1e-6 * den.max(); den = np.where(okd, den, 1.)
def coef(z): return np.where(okd, z / den, 1.)[shb]
v = np.zeros(inm.sum())
for j in range(len(L)):
    v += w[j] * back(EF[j] - coef(ssum(np.real(EF[j] * OuF.conj()))) * EuF) * back(OF[j] - coef(ssum(np.real(OF[j] * EuF.conj()))) * OuF)
gw = np.minimum(np.arange(n), n - np.arange(n)); Rw = np.sqrt(gw[:, None, None] ** 2 + gw[None, :, None] ** 2 + gw[None, None, :] ** 2)
def conv(x, k): return np.fft.irfftn(np.fft.rfftn(x) * np.fft.rfftn(k / k.sum()), s=(n, n, n))
def full(x): y = np.zeros((n, n, n)); y[inm] = x; return y
v = conv(full(v), np.exp(-Rw ** 2 / 8.))[inm]; thr = np.quantile(-v, 1 - a.alpha); sig = v > thr
vs = np.where(sig, v, 0.); o = np.argsort(-vs); c = np.cumsum(vs[o]) / max(vs.sum(), 1e-300); keep = np.zeros(len(v), bool); keep[o[:np.searchsorted(c, a.mass) + 1]] = True; keep &= sig
ball = (Rw <= 3).astype(float); mk = (conv(full(keep.astype(float)), ball) > 0.5 / ball.sum()).astype(float); mk = np.clip(conv(mk, np.exp(-Rw ** 2 / 8.)), 0, 1) * M
wr(a.out, mk); print('focus mask %s: %d members, %.1f %% of the sphere (%d significant voxels, threshold %.3g)' % (a.out, len(L), 100 * mk.sum() / M.sum(), sig.sum(), thr))
