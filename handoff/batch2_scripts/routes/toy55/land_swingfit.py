#!/usr/bin/env python3
"""land_swingfit.py: position of maps on a PHYSICAL 1-D arm coordinate (the rigid swing of make_swing_maps.py: the closed map's mobile domain
rotated about its hinge, step k = tip travel). Every map is low-passed to LP and standardised as arm_readout.py; inside the swept-volume mask
(voxels where any swing map differs from step 0 by > THR of the largest difference, dilated) the map is fitted by NNLS as sum_k q_k S_k + offset;
reported: q-weighted mean step, q spread (sd, A), the two largest components, relative residual, and the fit's R^2 over the mask.
usage: land_swingfit.py SWINGDIR OUT.txt MAP|DIR [...] [--lp 10 --thr 0.15 --dil 2]"""
import sys, os, glob, argparse, numpy as np
from scipy.optimize import nnls
from scipy.ndimage import binary_dilation
ap = argparse.ArgumentParser(); ap.add_argument('swing'); ap.add_argument('out'); ap.add_argument('maps', nargs='+'); ap.add_argument('--lp', type=float, default=10.)
ap.add_argument('--thr', type=float, default=0.15); ap.add_argument('--dil', type=int, default=2); ap.add_argument('--smpd', type=float, default=1.464); a = ap.parse_args()
def rd(p):
    h = np.fromfile(p, np.int32, 256); n = int(h[0]); return np.fromfile(p, np.float32, offset=1024 + int(h[23]), count=n ** 3).reshape(n, n, n).astype(np.float64)
tr = [l.split() for l in open(a.swing + '/truth_swing.txt') if l.strip() and not l.startswith('#')]; tips = np.array([float(l[3]) for l in tr])
SF = [glob.glob('%s/conf%s_*.mrc' % (a.swing, l[1]))[0] for l in tr]
n = rd(SF[0]).shape[0]; f = np.fft.fftfreq(n, d=a.smpd); kx, ky, kz = np.meshgrid(f, f, f, indexing='ij'); G = np.exp(-0.5 * (kx ** 2 + ky ** 2 + kz ** 2) * (a.lp / 1.2) ** 2)
lp = lambda v: np.real(np.fft.ifftn(np.fft.fftn(v) * G))
S = np.array([lp(rd(p)) for p in SF]); ref = S.mean(0); supp = ref > 0.2 * ref.max(); outside = ref < 0.02 * ref.max()
def std(v): v = v - v[outside].mean(); return v / max(np.sqrt(np.mean(v[supp] ** 2)), 1e-12)
S = np.array([std(s) for s in S]); dev = np.abs(S - S[0]).max(0); M = binary_dilation(dev > a.thr * dev.max(), iterations=a.dil)
A = np.c_[S[:, M].T, np.ones(M.sum()), -np.ones(M.sum())]
print('swing basis: %d steps, tip travel %.0f..%.0f A; swept-volume mask %d voxels (%.1f %% of the support)' % (len(SF), tips[0], tips[-1], M.sum(), 100. * M.sum() / supp.sum()))
files = []
for p in a.maps: files += sorted(glob.glob(p + '/*.mrc')) if os.path.isdir(p) else [p]
rows = []
for fn in files:
    v = std(lp(rd(fn)))[M]; x, rn = nnls(A, v); q = x[:len(SF)]; tot = q.sum()
    if tot <= 0: rows.append((fn, np.nan, np.nan, '-', np.nan, np.nan)); continue
    q = q / tot; mt = q @ tips; sd = np.sqrt(q @ (tips - mt) ** 2); o = np.argsort(-q)[:2]
    fit = A @ x; r2 = 1. - np.sum((v - fit) ** 2) / np.sum((v - v.mean()) ** 2)
    rows.append((fn, mt, sd, '%.0f:%.2f,%.0f:%.2f' % (tips[o[0]], q[o[0]], tips[o[1]], q[o[1]]), rn / np.linalg.norm(v), r2))
with open(a.out, 'w') as fo:
    fo.write('# map mean_tip_travel_A sd_A top2(travel:weight) rel_residual R2 (NNLS on the rigid swing %s, lp %.0f)\n' % (a.swing, a.lp))
    for r in rows: fo.write('%s %.2f %.2f %s %.3f %.3f\n' % r)
for r in rows: print('%-62s travel %6.2f sd %5.2f  %-18s resid %.3f R2 %.3f' % (os.path.basename(r[0]), r[1], r[2], r[3], r[4], r[5]))
