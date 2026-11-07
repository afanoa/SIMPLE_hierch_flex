#!/usr/bin/env python3
"""Closed-arm occupancy of every cell map of a klem run, on the run's coarse grid: the deposited closed conformer (18)
and the consensus are Fourier-cropped to the cell grid, the arm = voxels where 18 exceeds the consensus by > 3 sd inside
18's density, the core = dense voxels where they agree; occupancy = (cell[arm] / cell[core]) relative to the same ratio
in 18. Also the deviation correlation with 18 at 60-20 A. usage: closed_arm_cells.py KLEM_DIR [--top 5]"""
import sys, glob, re, argparse, numpy as np
from scipy.ndimage import gaussian_filter
ap = argparse.ArgumentParser(); ap.add_argument('klem'); ap.add_argument('--top', type=int, default=5); a = ap.parse_args()
W = '/mnt/beegfs/elmlund/afan/routes/toy55'; SMPD = 1.464
def read_mrc(p):
    h = np.fromfile(p, dtype=np.int32, count=256); nx, ny, nz = h[0], h[1], h[2]; return np.fromfile(p, dtype=np.float32, offset=1024 + h[23]).reshape(nz, ny, nx).astype(float)
fs = sorted(glob.glob(a.klem + '/cells/cell_*.mrc')) + sorted(glob.glob(a.klem + '/cells/state_*.mrc')); V = [read_mrc(f) for f in fs]; nc = V[0].shape[0]
cons = read_mrc(W + '/maps160/toy_consensus.mrc'); c18 = read_mrc(glob.glob(W + '/maps160/docked/conf18_*.mrc')[0]); n = cons.shape[0]
zz, yy, xx = np.indices((n, n, n)); r = np.hypot(np.hypot(xx - n / 2, yy - n / 2), zz - n / 2); msk = 1. / (1. + np.exp((r - 100 / SMPD) / 2.))
def crop(v):
    Fv = np.fft.fftshift(np.fft.fftn(np.fft.ifftshift(v * msk))); c = n // 2; h = nc // 2; sub = Fv[c - h:c - h + nc, c - h:c - h + nc, c - h:c - h + nc]; return np.real(np.fft.fftshift(np.fft.ifftn(np.fft.ifftshift(sub)))) * (nc / n) ** 3
smpd_c = n * SMPD / nc; sig = 2.5 * SMPD / smpd_c; g = lambda v: gaussian_filter(v, sig)
L18, Lc = g(crop(c18)), g(crop(cons)); dens18 = L18 > L18.mean() + 2.0 * L18.std(); diff = L18 - Lc
arm = (diff > diff.mean() + 3 * diff.std()) & dens18; core = dens18 & (np.abs(diff) < 0.5 * diff.std()); ref = L18[arm].mean() / L18[core].mean()
zc, yc, xc = np.indices((nc, nc, nc)); rc = np.hypot(np.hypot(xc - nc / 2, yc - nc / 2), zc - nc / 2); band = (rc >= nc * smpd_c / 60.) & (rc <= nc * smpd_c / 20.)
Fc = lambda v: np.fft.fftshift(np.fft.fftn(np.fft.ifftshift(v)))[band]; T18 = Fc(crop(c18) - crop(cons)); T18 /= np.sqrt(np.sum(np.abs(T18) ** 2))
zc2, yc2, xc2 = np.indices((nc, nc, nc)); rA = np.hypot(np.hypot(xc2 - nc / 2, yc2 - nc / 2), zc2 - nc / 2) * smpd_c; tip = arm & (rA > 95); base = arm & (rA < 80); rtip = L18[tip].mean() / L18[core].mean(); rbase = L18[base].mean() / L18[core].mean()
Vm = np.mean(V, 0); occ = []; cc = []; tb = []
for v in V:
    Lv = g(v); occ.append((Lv[arm].mean() / Lv[core].mean()) / ref); tb.append(((Lv[tip].mean() / Lv[core].mean()) / rtip, (Lv[base].mean() / Lv[core].mean()) / rbase)); X = Fc(v - Vm); cc.append(float(np.real(np.sum(np.conj(X) * T18)) / max(np.sqrt(np.sum(np.abs(X) ** 2)), 1e-30)))
occ = np.array(occ); cc = np.array(cc); order = np.argsort(-occ)
print('closed-arm occupancy per cell (coarse grid %d^3, %.1f A/px; arm %d voxels, core %d; consensus itself %.2f): median %.2f, max %.2f' % (nc, smpd_c, arm.sum(), core.sum(), (Lc[arm].mean() / Lc[core].mean()) / ref, np.median(occ), occ.max()))
for k in order[:a.top]: print('  %s: occupancy %.2f (tip r>95 A %.2f, base r<80 A %.2f) | corr with 18 %+.2f' % (fs[k].split('/')[-1], occ[k], tb[k][0], tb[k][1], cc[k]))
np.savetxt(a.klem + '/closed_arm_cells.txt', np.c_[np.arange(1, len(V) + 1), occ, cc], fmt='%d %.4f %.4f', header='cell closed_arm_occupancy corr18')
