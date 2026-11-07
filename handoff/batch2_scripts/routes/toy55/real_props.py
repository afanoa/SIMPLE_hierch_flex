#!/usr/bin/env python3
"""Properties of the real RNase P particles the toy must match, from the dumped class images (box 80 @ 2.928 A) and
meta: per-particle correlation with the CTF-applied consensus projection per band (-> band SNR = r^2/(1-r^2)),
defocus distribution, view distribution (e2 tilt histogram, class populations), refined-pose spread per 2D class.
usage: real_props.py DATA_DIR CONS80 OUT [--nptcl 3000]"""
import sys, os, glob, argparse, numpy as np
from scipy.ndimage import map_coordinates
ap = argparse.ArgumentParser(); ap.add_argument('data'); ap.add_argument('cons'); ap.add_argument('out'); ap.add_argument('--nptcl', type=int, default=3000)
a = ap.parse_args(); rng = np.random.default_rng(0); os.makedirs(a.out, exist_ok=True)
KV, CS, FRACA = 300., 2.7, 0.1
def read_mrc(fn):
    with open(fn, 'rb') as f: h = np.frombuffer(f.read(1024), np.int32)
    return np.array(np.memmap(fn, np.float32, 'r', offset=1024 + int(h[23]), shape=(int(h[2]), int(h[1]), int(h[0]))), dtype=np.float64)
def fftc(x): return np.fft.fftshift(np.fft.fftn(np.fft.ifftshift(x)))
def rotmat(ang, axis):
    r = np.radians(ang); c, s = np.cos(r), np.sin(r); R = np.eye(3)
    if axis == 3: R[0, 0], R[0, 1], R[1, 0], R[1, 1] = c, -s, s, c
    else:         R[0, 0], R[0, 2], R[2, 0], R[2, 2] = c, -s, s, c
    return R
def euler2m(e): return rotmat(e[2], 3) @ rotmat(e[1], 2) @ rotmat(e[0], 3)
def ctf_eval(kx, ky, dfx, dfy, angast, box, smpd):
    kv = KV * 1e3; wl = 12.26 / np.sqrt(kv + 0.9784 * kv * kv * 1e-6) / smpd; cs = CS * 1e7 / smpd
    dfxp = dfx * 1e4 / smpd; dfyp = dfy * 1e4 / smpd; ang = np.arctan2(ky, kx)
    df = 0.5 * (dfxp + dfyp + np.cos(2. * (ang - np.radians(angast))) * (dfxp - dfyp)); k2 = (kx * kx + ky * ky) / float(box * box)
    return np.sin(np.pi * wl * k2 * (df - 0.5 * wl * wl * k2 * cs) + np.arctan(FRACA / np.sqrt(1. - FRACA * FRACA)))
cons = read_mrc(a.cons); BOX = cons.shape[0]; VC = fftc(cons)
metas = sorted(glob.glob(os.path.join(a.data, 'cls*_meta.npz'))); SMPD = float(np.load(metas[0])['smpd'])
ks = np.arange(-BOX // 2 + 1, BOX // 2); KX, KY = np.meshgrid(ks, ks, indexing='xy'); RR = np.hypot(KX, KY)
hm = (RR <= BOX // 2 - 1) & ((KX > 0) | ((KX == 0) & (KY >= 0))); kx = KX[hm].astype(float); ky = KY[hm].astype(float); shell = np.rint(RR[hm]).astype(int)
ix = (kx + BOX // 2).astype(int); iy = (ky + BOX // 2).astype(int); wq = np.where(kx > 0, 2., 1.)
BANDS = [(117, 30), (30, 20), (20, 12), (12, 8), (8, 6)]
BM = [wq * ((shell >= BOX * SMPD / lo) & (shell <= BOX * SMPD / hi)) for lo, hi in BANDS]
yy, xx = np.mgrid[:BOX, :BOX]; rr = np.hypot(xx - BOX / 2, yy - BOX / 2); RM = BOX / 2 - 1; MASK = 0.5 * (1 + np.cos(np.pi * np.clip((rr - (RM - 6)) / 6., 0, 1)))
CONV, SS, FS = 0, -1., 0.5
def proj(p):
    R = euler2m([-p[0], p[1], -p[2]]); c = BOX // 2; loc = np.stack([kx, ky, np.zeros_like(kx)], 1) @ (R if CONV == 0 else R.T)
    co = [loc[:, 2] + c, loc[:, 1] + c, loc[:, 0] + c]
    P = map_coordinates(VC.real, co, order=1, mode='constant', cval=0.) + 1j * map_coordinates(VC.imag, co, order=1, mode='constant', cval=0.)
    return P * np.exp(-2j * np.pi * (kx * SS * FS * p[3] + ky * SS * FS * p[4]) / BOX)
def wdot(x, y, w): return float(np.real(np.sum(w * np.conj(x) * y)))
# sample particles across classes
rows = []
per = max(1, a.nptcl // len(metas))
for mf in metas:
    M = np.load(mf); c = int(os.path.basename(mf)[3:6]); imgs = np.load(mf.replace('_meta.npz', '_imgs.npy'), mmap_mode='r')
    refined = ~((M['e3'] == 0) & (M['x'] == 0) & (M['y'] == 0)); idx = np.where(refined)[0]
    if len(idx) == 0: continue
    for i in rng.choice(idx, min(per, len(idx)), replace=False):
        im = np.array(imgs[i], float); im = (im - im[rr > RM].mean()) * MASK; F = fftc(im)[iy, ix]
        p = [M[k][i] for k in ('e1', 'e2', 'e3', 'x', 'y')]; C = np.abs(ctf_eval(kx, ky, M['dfx'][i], M['dfy'][i], M['angast'][i], BOX, SMPD))
        P = C * proj(p)
        r = [wdot(F, P, W) / max(np.sqrt(wdot(F, F, W) * wdot(P, P, W)), 1e-30) for W in BM]
        rows.append([c, M['dfx'][i], M['dfy'][i], M['e1'][i], M['e2'][i]] + r)
R = np.array(rows); np.savez(os.path.join(a.out, 'real_props.npz'), rows=R)
print('%d particles from %d classes' % (len(R), len(np.unique(R[:, 0]))))
print('defocus (um): dfx mean %.2f sd %.2f, range p5-p95 %.2f-%.2f; astig |dfx-dfy| median %.3f' % (R[:, 1].mean(), R[:, 1].std(), *np.percentile(R[:, 1], (5, 95)), np.median(np.abs(R[:, 1] - R[:, 2]))))
for b, (lo, hi) in enumerate(BANDS):
    r = R[:, 5 + b]; r2 = np.clip(np.median(r), 0, 0.999) ** 2
    print('band %3d-%2d A: per-particle corr with the CTF-applied consensus projection median %.3f (p25 %.3f p75 %.3f) -> band SNR %.4f' % (lo, hi, np.median(r), *np.percentile(r, (25, 75)), r2 / (1 - r2)))
e2 = R[:, 4]; h, _ = np.histogram(e2, bins=np.arange(0, 181, 15)); print('tilt (e2) histogram, 15-deg bins from 0: %s' % h.tolist())
cls = R[:, 0].astype(int); print('particles per sampled class: median %d' % np.median(np.bincount(cls)[np.bincount(cls) > 0]))
