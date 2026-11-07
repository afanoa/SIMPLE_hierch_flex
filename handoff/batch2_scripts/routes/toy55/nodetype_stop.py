#!/usr/bin/env python3
"""Stop test of one node from its cheap children (cheapkids output): R = cross-half reproducibility of the children's per-shell-
regressed deviations in a band the assignment never saw (default 15-11 A; the flex latent assigning the children uses <= 15 A), with
its noise: R recomputed in 8 Fourier-space sectors (signs of kx, ky and low / high |kz|), SE = sd(R_sector) / sqrt(8), z = R / SE.
Expected under no resolvable heterogeneity: R ~ 0, |z| small. Prints one line per band.
usage: nodetype_stop.py KIDSDIR [--bands 15:11,15:8,11:8] [--mskdiam 220] [--tag NAME]"""
import argparse, glob, os, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('kids'); ap.add_argument('--bands', default='15:11,15:8,11:8'); ap.add_argument('--mskdiam', type=float, default=220.); ap.add_argument('--tag', default='')
a = ap.parse_args()
def rd(p):
    h = np.fromfile(p, np.int32, 256); n = int(h[0]); return np.fromfile(p, np.float32, offset=1024 + int(h[23]), count=n ** 3).reshape(n, n, n).astype(np.float64), float(np.fromfile(p, np.float32, 13)[10]) / n
lab = np.loadtxt(os.path.join(a.kids, 'children.txt')).astype(int); RR = sorted(glob.glob(os.path.join(a.kids, 'rec', '*_reconstruct3D')))[-1]; K = lab.max()
w = np.array([(lab == k).sum() for k in range(1, K + 1)], float); w /= w.sum()
v0, smpd = rd(os.path.join(RR, 'recvol_state01_even_unfil.mrc')); n = v0.shape[0]; c = n // 2
g = np.arange(n) - c; R3 = np.sqrt(g[:, None, None] ** 2 + g[None, :, None] ** 2 + g[None, None, :] ** 2); r0 = a.mskdiam / 2 / smpd
M = np.where(R3 <= r0, 1., np.where(R3 <= r0 + 6, 0.5 + 0.5 * np.cos(np.pi * (R3 - r0) / 6), 0.))
f = np.fft.fftfreq(n); kx, ky, kz = np.meshgrid(f, f, np.fft.rfftfreq(n), indexing='ij'); kr = np.sqrt(kx ** 2 + ky ** 2 + kz ** 2); res = np.where(kr > 0, smpd / np.maximum(kr, 1e-9), 1e9)
keep = res >= 7.9; sh = np.round(kr * n).astype(int)[keep]; nsh = sh.max() + 1; wb = np.where((kz == 0) | (np.abs(kz) == 0.5), 1., 2.)[keep]; resk = res[keep]
sec = ((kx[keep] >= 0).astype(int) + 2 * (ky[keep] >= 0) + 4 * (np.abs(kz[keep]) >= np.abs(kz[keep]).max() / 2)).astype(int)
Fm = lambda x: np.fft.rfftn(M * x)[keep]
E = np.array([Fm(rd(os.path.join(RR, 'recvol_state%02d_even_unfil.mrc' % k))[0]) for k in range(1, K + 1)]); O = np.array([Fm(rd(os.path.join(RR, 'recvol_state%02d_odd_unfil.mrc' % k))[0]) for k in range(1, K + 1)])
Em, Om = w @ E, w @ O; ss = lambda z: np.bincount(sh, weights=wb * z, minlength=nsh); den = ss(np.real(Em * Om.conj())); ok = np.abs(den) > 1e-9 * np.abs(den).max(); den = np.where(ok, den, 1.)
dE = np.array([e - np.where(ok, ss(np.real(e * Om.conj())) / den, 1.)[sh] * Em for e in E]); dO = np.array([o - np.where(ok, ss(np.real(o * Em.conj())) / den, 1.)[sh] * Om for o in O])
x = (w[:, None] * np.real(dE * dO.conj())).sum(0) * wb; pe = (w[:, None] * np.abs(dE) ** 2).sum(0) * wb; po = (w[:, None] * np.abs(dO) ** 2).sum(0) * wb
tag = (a.tag + ' ') if a.tag else ''
for b in a.bands.split(','):
    lo, hi = [float(t) for t in b.split(':')]; m = (resk <= lo) & (resk >= hi)
    Rv = x[m].sum() / np.sqrt(pe[m].sum() * po[m].sum()); Rs = np.array([x[m & (sec == s)].sum() / np.sqrt(pe[m & (sec == s)].sum() * po[m & (sec == s)].sum()) for s in range(8)])
    se = Rs.std(ddof=1) / np.sqrt(8); print('%sband %g-%g A: R %.4f  SE %.4f  z %.1f   (sector R %s)' % (tag, lo, hi, Rv, se, Rv / se, ' '.join('%.3f' % r for r in Rs)))
