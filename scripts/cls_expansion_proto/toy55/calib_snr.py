#!/usr/bin/env python3
"""Band correlation of simulated particles with the CTF-applied projection of the simulation volume at the true poses,
the same statistic as real_props.py on the real data. usage: calib_snr.py STK ORIS VOL [--smpd 1.464]"""
import sys, argparse, numpy as np
from scipy.ndimage import map_coordinates
ap = argparse.ArgumentParser(); ap.add_argument('stk'); ap.add_argument('oris'); ap.add_argument('vol'); ap.add_argument('--smpd', type=float, default=1.464); ap.add_argument('--nmax', type=int, default=1500); ap.add_argument('--abs', action='store_true', help='|CTF| model (phase-flipped stacks); default signed CTF (ctf=yes stacks)')
a = ap.parse_args(); KV, CS, FRACA = 300., 2.7, 0.1
def read_mrc(fn):
    h = np.fromfile(fn, np.int32, 256); nx, ny, nz = h[:3]; return np.memmap(fn, np.float32, 'r', offset=1024 + int(h[23]), shape=(int(nz), int(ny), int(nx)))
def rotmat(ang, axis):
    r = np.radians(ang); c, s = np.cos(r), np.sin(r); R = np.eye(3)
    if axis == 3: R[0, 0], R[0, 1], R[1, 0], R[1, 1] = c, -s, s, c
    else:         R[0, 0], R[0, 2], R[2, 0], R[2, 2] = c, -s, s, c
    return R
def euler2m(e): return rotmat(e[2], 3) @ rotmat(e[1], 2) @ rotmat(e[0], 3)
def ctf(kx, ky, dfx, dfy, angast, box, smpd):
    kv = KV * 1e3; wl = 12.26 / np.sqrt(kv + 0.9784 * kv * kv * 1e-6) / smpd; cs = CS * 1e7 / smpd
    dfxp = dfx * 1e4 / smpd; dfyp = dfy * 1e4 / smpd; ang = np.arctan2(ky, kx); df = 0.5 * (dfxp + dfyp + np.cos(2. * (ang - np.radians(angast))) * (dfxp - dfyp)); k2 = (kx * kx + ky * ky) / float(box * box)
    return np.sin(np.pi * wl * k2 * (df - 0.5 * wl * wl * k2 * cs) + np.arctan(FRACA / np.sqrt(1. - FRACA * FRACA)))
stk = read_mrc(a.stk); V = np.array(read_mrc(a.vol), float); BOX = V.shape[0]; VC = np.fft.fftshift(np.fft.fftn(np.fft.ifftshift(V)))
oris = [dict(t.split('=', 1) for t in l.split() if '=' in l) for l in open(a.oris) if l.strip()]
ks = np.arange(-BOX // 2 + 1, BOX // 2); KX, KY = np.meshgrid(ks, ks, indexing='xy'); RR = np.hypot(KX, KY); hm = (RR <= BOX // 2 - 1) & ((KX > 0) | ((KX == 0) & (KY >= 0)))
kx, ky, shell = KX[hm].astype(float), KY[hm].astype(float), np.rint(RR[hm]).astype(int); ix = (kx + BOX // 2).astype(int); iy = (ky + BOX // 2).astype(int); wq = np.where(kx > 0, 2., 1.)
BANDS = [(117, 30), (30, 20), (20, 12), (12, 8), (8, 6)]; BM = [wq * ((shell >= BOX * a.smpd / lo) & (shell <= BOX * a.smpd / hi)) for lo, hi in BANDS]
yy, xx = np.mgrid[:BOX, :BOX]; rr = np.hypot(xx - BOX / 2, yy - BOX / 2); RM = BOX / 2 - 1; MASK = 0.5 * (1 + np.cos(np.pi * np.clip((rr - (RM - 6)) / 6., 0, 1)))
res = []
for conv in (0, 1):
    for ss in (1., -1.):
        cc = []
        for i in range(min(a.nmax, len(oris))):
            d = oris[i]; p = [float(d.get(k, 0)) for k in ('e1', 'e2', 'e3', 'x', 'y')]
            im = np.array(stk[i], float); im = (im - im[rr > RM].mean()) * MASK; F = np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(im)))[iy, ix]
            R = euler2m([-p[0], p[1], -p[2]]); M = R if conv == 0 else R.T; c = BOX // 2; loc = np.stack([kx, ky, np.zeros_like(kx)], 1) @ M
            P = map_coordinates(VC.real, [loc[:, 2] + c, loc[:, 1] + c, loc[:, 0] + c], order=1, mode='constant', cval=0.) + 1j * map_coordinates(VC.imag, [loc[:, 2] + c, loc[:, 1] + c, loc[:, 0] + c], order=1, mode='constant', cval=0.)
            P = P * np.exp(-2j * np.pi * ss * (kx * p[3] + ky * p[4]) / BOX) * (np.abs if a.abs else (lambda v: v))(ctf(kx, ky, float(d['dfx']), float(d['dfy']), float(d['angast']), BOX, a.smpd))
            cc.append([float(np.real(np.sum(W * np.conj(F) * P)) / max(np.sqrt(np.real(np.sum(W * np.abs(F) ** 2)) * np.real(np.sum(W * np.abs(P) ** 2))), 1e-30)) for W in BM])
        cc = np.array(cc); res.append((np.median(cc[:, 0]), conv, ss, np.median(cc, 0)))
best = max(res, key=lambda t: abs(t[0])); sgn = 1 if best[0] > 0 else -1
print('BANDCORR ' + ' '.join('%.3f' % (sgn * v) for v in best[3]) + '   (conv %d shift %+.0f, n %d)' % (best[1], best[2], min(a.nmax, len(oris))))
