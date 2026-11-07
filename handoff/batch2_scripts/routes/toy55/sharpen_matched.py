#!/usr/bin/env python3
"""Sharpened copies of the ten matched-conformer maps: soft mask, B-factor (default -400 A^2), cosine low-pass at each
map's own FSC 0.143 resolution; conformer-named files + README with the per-map numbers.
usage: sharpen_matched.py MAPS_DIR OUT_DIR [--bfac -400]"""
import sys, glob, argparse, numpy as np
sys.path.insert(0, '/mnt/beegfs/elmlund/afan/routes/toy55'); from model_class import read_mrc, write_mrc, Shells, fftc
from scipy.ndimage import gaussian_filter
ap = argparse.ArgumentParser(); ap.add_argument('maps'); ap.add_argument('out'); ap.add_argument('--bfac', type=float, default=-400.); ap.add_argument('--rmask', type=float, default=116., help='soft mask radius in A (the closed arm reaches 114 A; 100 clipped its tip)'); a = ap.parse_args()
import os; os.makedirs(a.out, exist_ok=True); W = '/mnt/beegfs/elmlund/afan/routes/toy55'; SMPD = 1.464; n = 160; ids = [1, 10, 18, 27, 28, 37, 46, 48, 49, 51]
S = Shells(n, SMPD); msk = S.soft_mask(a.rmask / SMPD) if a.rmask > 0 else np.ones((n, n, n)); sh = S.sh; zz, yy, xx = np.indices((n, n, n)); rr = np.hypot(np.hypot(xx - n / 2, yy - n / 2), zz - n / 2); s_inv = rr / (n * SMPD)   # 1/A
cons = read_mrc(W + '/maps160/toy_consensus.mrc'); C = {i: read_mrc(glob.glob(W + '/maps160/docked/conf%02d_*.mrc' % i)[0]) for i in ids}
g = lambda v: gaussian_filter(v, 2.5); L18, Lc = g(C[18]), g(cons); dens18 = L18 > L18.mean() + 2.0 * L18.std(); diff = L18 - Lc
arm = (diff > diff.mean() + 3 * diff.std()) & dens18; core = dens18 & (np.abs(diff) < 0.5 * diff.std()); ref = L18[arm].mean() / L18[core].mean(); T = {i: fftc((C[i] - cons) * msk) for i in ids}
def bc(x, y, lo, hi):
    m = S.band(lo, hi); return float(np.real(np.sum(np.conj(x[m]) * y[m])) / max(np.sqrt(np.sum(np.abs(x[m]) ** 2) * np.sum(np.abs(y[m]) ** 2)), 1e-30))
def fsc_res(x, y):
    a_, b_ = fftc(x * msk), fftc(y * msk); num = np.bincount(sh.ravel(), np.real(np.conj(a_) * b_).ravel()); d1 = np.bincount(sh.ravel(), np.abs(a_.ravel()) ** 2); d2 = np.bincount(sh.ravel(), np.abs(b_.ravel()) ** 2)
    f = num / np.maximum(np.sqrt(d1 * d2), 1e-30); s = 1
    while s < 80 and f[s] > 0.143: s += 1
    return n * SMPD / s
V = [read_mrc('%s/recvol_state%02d.mrc' % (a.maps, k + 1)) for k in range(10)]; Vm = np.mean(V, 0); st = np.loadtxt(a.maps + '/../states.txt', dtype=int)
lines = ['# RNase P, real data: one full-box map per deposited conformer (family 1 of the klem tree, particles of the cells matching that conformer, refined poses)',
         '# conformer  particles  own_corr_60-20A  FSC0.143_A  closed_arm_occupancy  raw_file  sharpened_file (B %+.0f A^2, low-pass at the FSC resolution, soft mask radius %.0f A)' % (a.bfac, a.rmask)]
for k, i in enumerate(ids):
    try: res = fsc_res(read_mrc('%s/recvol_state%02d_even.mrc' % (a.maps, k + 1)), read_mrc('%s/recvol_state%02d_odd.mrc' % (a.maps, k + 1)))
    except Exception: res = 8.
    own = bc(fftc((V[k] - Vm) * msk), T[i], 60, 20); occ = (g(V[k])[arm].mean() / g(V[k])[core].mean()) / ref
    F = np.fft.fftshift(np.fft.fftn(np.fft.ifftshift(V[k] * msk))); F *= np.exp(-a.bfac * s_inv ** 2 / 4.)
    fc = 1. / res; w = 0.5 * (1 + np.cos(np.pi * np.clip((s_inv - 0.85 * fc) / (0.3 * fc), 0, 1))); F *= w
    vs = np.real(np.fft.fftshift(np.fft.ifftn(np.fft.ifftshift(F)))); raw = 'conf%02d_ours_raw.mrc' % i; shp = 'conf%02d_ours_sharp.mrc' % i
    write_mrc('%s/%s' % (a.out, raw), V[k].astype(np.float32), SMPD); write_mrc('%s/%s' % (a.out, shp), vs.astype(np.float32), SMPD)
    lines.append('%9d %9d %+15.2f %11.1f %20.2f  %s  %s' % (i, int((st == k + 1).sum()), own, res, occ, raw, shp))
open(a.out + '/README.txt', 'w').write('\n'.join(lines) + '\n'); print('\n'.join(lines))
