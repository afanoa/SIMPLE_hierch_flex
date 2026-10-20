#!/usr/bin/env python3
"""Identification readout on real RNase P: for every deposited conformer, the best-matching state of a klem run (cells
on the coarse grid vs the docked deposited conformers Fourier-cropped to that grid, deviation correlation at 60-20 A),
how many states have it as their best match, and the mean best-match correlation. The aim: every conformer gets a state.
usage: deposited_match.py CELLS_DIR DOCKED_DIR CONSENSUS [--lo 60 --hi 20]"""
import sys, glob, re, argparse, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('cells'); ap.add_argument('docked'); ap.add_argument('cons'); ap.add_argument('--lo', type=float, default=60.); ap.add_argument('--hi', type=float, default=20.); ap.add_argument('--save', default='', help='write per-state best conformer and corr (state best corr) plus the full matrix'); a = ap.parse_args(); SMPD = 1.464
def read_mrc(p):
    h = np.fromfile(p, dtype=np.int32, count=256); nx, ny, nz = h[0], h[1], h[2]; return np.fromfile(p, dtype=np.float32, offset=1024 + h[23]).reshape(nz, ny, nx).astype(float)
fs = sorted(glob.glob(a.cells + '/cell_*.mrc')) + sorted(glob.glob(a.cells + '/state_*.mrc')); V = np.array([read_mrc(f) for f in fs]); K = len(V); nc = V.shape[1]; Dv = V - (V.mean(0) if K > 1 else read_mrc(a.cons))   # one state: deviation from the consensus, not from itself
dep = sorted(glob.glob(a.docked + '/conf*_m160.mrc')); ids = [int(re.findall(r'conf(\d+)', f.split('/')[-1])[0]) for f in dep]; cons = read_mrc(a.cons); n = cons.shape[0]
zz, yy, xx = np.indices((n, n, n)); r = np.hypot(np.hypot(xx - n / 2, yy - n / 2), zz - n / 2); msk = 1. / (1. + np.exp((r - 100 / SMPD) / 2.))
def crop(v):
    Fv = np.fft.fftshift(np.fft.fftn(np.fft.ifftshift(v * msk))); c = n // 2; h = nc // 2; sub = Fv[c - h:c - h + nc, c - h:c - h + nc, c - h:c - h + nc]; return np.real(np.fft.fftshift(np.fft.ifftn(np.fft.ifftshift(sub))))
Dc = np.array([crop(read_mrc(f) - cons) for f in dep]); smpd_c = n * SMPD / nc; zc, yc, xc = np.indices((nc, nc, nc)); rc = np.hypot(np.hypot(xc - nc / 2, yc - nc / 2), zc - nc / 2); band = (rc >= nc * smpd_c / a.lo) & (rc <= nc * smpd_c / a.hi)
Fc = lambda v: np.fft.fftshift(np.fft.fftn(np.fft.ifftshift(v)))[band]
X = np.array([Fc(d) for d in Dv]); Y = np.array([Fc(d) for d in Dc]); X /= np.maximum(np.sqrt(np.sum(np.abs(X) ** 2, 1))[:, None], 1e-30); Y /= np.maximum(np.sqrt(np.sum(np.abs(Y) ** 2, 1))[:, None], 1e-30); S = np.real(X @ np.conj(Y).T)
best_of_state = S.argmax(1); print('identification vs the %d deposited conformers (deviation corr at %g-%g A on the coarse grid), %d states:' % (len(ids), a.lo, a.hi, K))
print('  conformer  best state corr  states whose best match it is  (coverage: a conformer with >= 1 state and corr >= 0.3 counts as identified)')
nid = 0
for j, i in enumerate(ids):
    nb = int(np.sum(best_of_state == j)); bc = S[:, j].max(); ok = nb >= 1 and bc >= 0.3; nid += ok
    print('  %5d      %+.2f            %2d   %s' % (i, bc, nb, 'identified' if ok else '-'))
print('  identified %d / %d; mean best-match corr over states %.2f' % (nid, len(ids), np.mean(S.max(1))))
if a.save:
    with open(a.save, 'w') as f:
        f.write('# state best_conformer best_corr | corr per conformer %s\n' % ' '.join(str(i) for i in ids))
        for k in range(K): f.write('%d %d %.4f %s\n' % (k + 1, ids[best_of_state[k]], S[k].max(), ' '.join('%.3f' % v for v in S[k])))
    print('  per-state matches -> %s' % a.save)
