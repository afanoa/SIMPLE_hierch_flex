#!/usr/bin/env python3
"""lsamp_viewcheck.py: does a landscape coordinate follow the viewing direction or defocus (a nuisance) rather than the structure? Per run:
12 quantile bins of the full coordinate (coords.npz, first column) against 40 view cells (the scoring's Fibonacci cells of e1/e2) and
against defocus tertiles: Cramer's V of each table (0 = independent; for scale, V of random 12-bin labels is printed too).
usage: lsamp_viewcheck.py NAME:RUNDIR:POSES [...]"""
import sys, numpy as np
def kv(l): return dict(t.split('=', 1) for t in l.split() if '=' in t)
def cramer(x, y):
    _, xi = np.unique(x, return_inverse=True); _, yi = np.unique(y, return_inverse=True); T = np.zeros((xi.max() + 1, yi.max() + 1)); np.add.at(T, (xi, yi), 1)
    E = T.sum(1, keepdims=True) * T.sum(0, keepdims=True) / T.sum(); chi2 = ((T - E) ** 2 / np.maximum(E, 1e-12)).sum(); return np.sqrt(chi2 / (T.sum() * (min(T.shape) - 1)))
P = {}; rng = np.random.default_rng(0)
for spec in sys.argv[1:]:
    nm, rd, ps = spec.split(':'); C = np.load(rd + '/coords.npz'); rows = C['rows']
    if ps not in P: P[ps] = np.array([[float(kv(l).get(k, 0)) for k in ('e1', 'e2', 'dfx', 'dfy')] for l in open(ps)])
    pp = P[ps][rows]; ph, th = np.radians(pp[:, 0]), np.radians(pp[:, 1]); v = np.c_[np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)]
    i = np.arange(40) + 0.5; zz = 1 - 2 * i / 40; rr = np.sqrt(np.maximum(1 - zz * zz, 0)); gg = np.pi * (3 - np.sqrt(5)) * i
    vb = np.argmax(v @ np.c_[rr * np.cos(gg), rr * np.sin(gg), zz].T, 1); dfm = 0.5 * (pp[:, 2] + pp[:, 3]); dq = np.searchsorted(np.quantile(dfm, [1 / 3., 2 / 3.]), dfm)
    t = C['tfull'][:, 0]; b = np.searchsorted(np.quantile(t, np.linspace(0, 1, 13)[1:-1]), t); rb = rng.permutation(b)
    print('%-22s view cells V %.3f | defocus tertiles V %.3f | random labels: view V %.3f, defocus V %.3f' % (nm, cramer(b, vb), cramer(b, dq), cramer(rb, vb), cramer(rb, dq)), flush=True)
