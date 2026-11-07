#!/usr/bin/env python3
"""land_confmix.py: composition of maps in terms of a reference set (e.g. the deposited conformers). Maps and references are low-passed to
LP and standardised as arm_readout.py; inside the consensus support (dilated) every map is fitted by NNLS as sum_j q_j R_j + offset. Each
reference carries its closed-site occupancy from an arm_readout table -> per map: top weights, the weighted closed-site occupancy, and class
shares (closed >= 0.8, intermediate 0.45..0.8, open < 0.45). With --counts (particles per map) the class populations of the pooled set.
usage: land_confmix.py CONS ARMTABLE OUT.txt --refs F [F ...] --maps F [F ...] [--counts n1,n2,.. --lp 10 --dil 3]"""
import os, argparse, numpy as np
from scipy.optimize import nnls
from scipy.ndimage import binary_dilation
ap = argparse.ArgumentParser(); ap.add_argument('cons'); ap.add_argument('arm'); ap.add_argument('out'); ap.add_argument('--refs', nargs='+', required=True); ap.add_argument('--maps', nargs='+', required=True)
ap.add_argument('--counts', default=''); ap.add_argument('--lp', type=float, default=10.); ap.add_argument('--dil', type=int, default=3); ap.add_argument('--smpd', type=float, default=1.464); a = ap.parse_args()
def rd(p):
    h = np.fromfile(p, np.int32, 256); n = int(h[0]); return np.fromfile(p, np.float32, offset=1024 + int(h[23]), count=n ** 3).reshape(n, n, n).astype(np.float64)
c0 = rd(a.cons); n = c0.shape[0]; f = np.fft.fftfreq(n, d=a.smpd); kx, ky, kz = np.meshgrid(f, f, f, indexing='ij'); G = np.exp(-0.5 * (kx ** 2 + ky ** 2 + kz ** 2) * (a.lp / 1.2) ** 2)
lp = lambda v: np.real(np.fft.ifftn(np.fft.fftn(v) * G))
cons = lp(c0); supp = cons > 0.2 * cons.max(); outside = cons < 0.02 * cons.max(); M = binary_dilation(supp, iterations=a.dil)
def std(v): v = v - v[outside].mean(); return v / max(np.sqrt(np.mean(v[supp] ** 2)), 1e-12)
arm = {}
for l in open(a.arm):
    if l.startswith('#') or not l.strip(): continue
    w = l.split(); arm[os.path.basename(w[0])] = (float(w[1]), float(w[2]))
R = np.array([std(lp(rd(p)))[M] for p in a.refs]); RN = [os.path.basename(p).replace('_m160.mrc', '') for p in a.refs]
rc = np.array([arm.get(os.path.basename(p), (np.nan, np.nan))[0] for p in a.refs])
cls = np.where(rc >= 0.8, 0, np.where(rc >= 0.45, 1, 2)); A = np.c_[R.T, np.ones(M.sum()), -np.ones(M.sum())]
print('%d references (closed %d, intermediate %d, open %d by their closed-site occupancy), mask %d voxels' % (len(RN), (cls == 0).sum(), (cls == 1).sum(), (cls == 2).sum(), M.sum()))
cnt = np.array([float(x) for x in a.counts.split(',')]) if a.counts else None; rows = []; S = []
for i, p in enumerate(a.maps):
    v = std(lp(rd(p)))[M]; x, rn = nnls(A, v); q = x[:len(RN)]; q = q / max(q.sum(), 1e-12); o = np.argsort(-q)[:4]; sh = np.array([q[cls == k].sum() for k in range(3)])
    r2 = 1. - np.sum((v - A @ x) ** 2) / np.sum((v - v.mean()) ** 2); S.append(sh)
    rows.append('%s weighted_closed %.3f shares closed %.3f intermediate %.3f open %.3f R2 %.3f top %s' % (os.path.basename(p), q @ np.nan_to_num(rc), sh[0], sh[1], sh[2], r2,
                ','.join('%s:%.2f' % (RN[k], q[k]) for k in o)))
    print(rows[-1])
if cnt is not None and len(cnt) == len(S):
    S = np.array(S); pop = (cnt[:, None] * S).sum(0) / cnt.sum(); rows.append('# POOLED (%d particles): closed %.3f intermediate %.3f open %.3f' % (cnt.sum(), *pop)); print(rows[-1])
open(a.out, 'w').write('\n'.join(rows) + '\n')
