#!/usr/bin/env python3
"""diag_shells.py: per-shell diagnostics of the port on one node (same args as ppinfo_site.py plus --vol): cells hit, mean den,
eps, map power and even/odd FSC per shell; the ring-wise agreement of the raw planes with the slices of a saved map."""
import sys, os, argparse, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ppinfo_site as S
from mrcio import read_mrc, fcrop3
ap = argparse.ArgumentParser()
for k in ('node', 'dump', 'mask', 'out', 'site', 'cache', 'vol'): ap.add_argument('--' + k, default='')
ap.add_argument('--dirs', default='1,3'); ap.add_argument('--box', type=int, default=48); ap.add_argument('--emlim', type=float, default=15.)
ap.add_argument('--mskdiam', type=float, default=225.456); ap.add_argument('--poses', default='input'); ap.add_argument('--maskmode', default='soft')
ap.add_argument('--cachebox', type=int, default=88); ap.add_argument('--contrast', type=int, default=1); ap.add_argument('--shiftsign', type=float, default=1.)
ap.add_argument('--labels', default='em'); ap.add_argument('--truth', default=''); ap.add_argument('--rare', type=int, default=0); ap.add_argument('--nproc', type=int, default=1); ap.add_argument('--latent', default=''); ap.add_argument('--geomcheck', default='')
a = ap.parse_args()
S.setup(a)
G = S.G; B, nh, nyq = G['B'], G['nh'], G['nyq']
v, smpd = read_mrc(a.vol)
if v.shape[0] != B: v = fcrop3(v, B)
V = np.fft.rfftn(np.fft.ifftshift(v)).astype(np.complex64).ravel()
ia = np.where(G['eo'] == 0)[0]; ib = np.where(G['eo'] == 1)[0]; sub = np.arange(len(G['sw']))
epsa = S.density(ia, sub)
na, da = S.nn_acc(ia, np.ones(len(ia), np.float32), sub); nb, db = S.nn_acc(ib, np.ones(len(ib), np.float32), sub)
ma = S.nn_norm(na, da, epsa, len(ia)); mb = S.nn_norm(nb, db, epsa, len(ib))
cs = G['cell_shell']
print('shell ncells nhit  mean_den   eps*n     |ma|^2    |V|^2   fsc(ma,mb) fsc(ma,V) fsc(num_a,num_b) fsc(num_a,V)')
for s in range(1, nyq + 1):
    m = cs == s; hit = m & (da > 0)
    def fsc(x, y): return float((x[m] * np.conj(y[m])).real.sum() / max(np.sqrt((np.abs(x[m]) ** 2).sum() * (np.abs(y[m]) ** 2).sum()), 1e-30))
    print('%2d %6d %6d %10.3e %10.3e %10.3e %10.3e %8.3f %8.3f %8.3f %8.3f' % (s, m.sum(), hit.sum(), da[hit].mean() if hit.any() else 0, epsa[s] * len(ia),
          (np.abs(ma[m]) ** 2).mean(), (np.abs(V[m]) ** 2).mean(), fsc(ma, mb), fsc(ma, V), fsc(na, nb), fsc(na, V)))
# ring-wise agreement of the raw planes with the map's slices: sum_i Re(X_i conj(CTF_i V[cell])) / norms per ring
print('ring  corr(CTF X, CTF^2 slice)  mean|CTF X|^2  mean CTF^2')
P, C, idx, cnj, sw, shell = G['P'], G['C'], G['idx'], G['cnj'], G['sw'], G['shell']
nsh = shell.max() + 1; xy = np.zeros(nsh); xx = np.zeros(nsh); yy = np.zeros(nsh); cc = np.zeros(nsh); cnt = np.zeros(nsh)
for b0 in range(0, 12000, 3000):
    pi = np.arange(b0, b0 + 3000); ii = idx[pi]; ok = ii < G['NC']; ii = np.where(ok, ii, 0)
    pv = V[ii]; pv = np.where(cnj[pi], np.conj(pv), pv) * C[pi] * ok   # CTF^2 x slice against P = CTF x X: the sign follows the data's CTF lobe
    dv = P[pi] * ok
    key = np.broadcast_to(shell[None, :], dv.shape).ravel()
    xy += np.bincount(key, (dv.real * pv.real + dv.imag * pv.imag).ravel(), nsh); xx += np.bincount(key, (np.abs(dv) ** 2).ravel(), nsh)
    yy += np.bincount(key, (np.abs(pv) ** 2).ravel(), nsh); cc += np.bincount(key, (C[pi] * ok).ravel(), nsh); cnt += np.bincount(key, ok.ravel().astype(float), nsh)
for s in range(1, nsh):
    print('%2d  %8.4f  %10.3e  %7.3f' % (s + 1, xy[s] / max(np.sqrt(xx[s] * yy[s]), 1e-30), xx[s] / max(cnt[s], 1), cc[s] / max(cnt[s], 1)))
