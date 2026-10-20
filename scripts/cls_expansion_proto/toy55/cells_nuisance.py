#!/usr/bin/env python3
"""Is a board's partition a nuisance partition? Per argmax cell: particles, mean pose error (true vs refined direction,
deg), mean image correlation, junk fraction; the between-cell spread of the cell means as a z-score against random
cells of the same sizes, and the correlation of the two latent coordinates with pose error / corr.
usage: cells_nuisance.py KLEM_LATENT TRUE_ORIS REFINED_ORIS GT_PARTICLES"""
import sys, numpy as np
lat, truef, reff, gtf = sys.argv[1:5]
def oris(fn, keys):
    out = []
    for l in open(fn):
        kv = dict(t.split('=', 1) for t in l.split() if '=' in t); out.append([float(kv.get(k, 0)) for k in keys])
    return np.array(out)
T = oris(truef, ('e1', 'e2')); Rr = oris(reff, ('e1', 'e2', 'corr')); conf = np.loadtxt(gtf, comments='#')[:, 1].astype(int)
def nvec(e1, e2): a, b = np.radians(e1), np.radians(e2); return np.c_[np.cos(a) * np.sin(b), np.sin(a) * np.sin(b), np.cos(b)]
ang = np.degrees(np.arccos(np.clip(np.sum(nvec(T[:, 0], T[:, 1]) * nvec(Rr[:, 0], Rr[:, 1]), 1), -1, 1))); corr = Rr[:, 2]; junk = (conf == 0).astype(float)
L = np.loadtxt(lat, comments='#'); pind = L[:, 0].astype(int) - 1; z1, z2, am = L[:, 1], L[:, 2], L[:, 3].astype(int)
a, c, j = ang[pind], corr[pind], junk[pind]; cells = np.unique(am); rng = np.random.default_rng(0)
def spread(v, lab):
    m = np.array([v[lab == k].mean() for k in cells]); return m.std()
print('%d particles in %d cells; whole subset: pose error %.2f deg, corr %.3f, junk %.3f' % (len(pind), len(cells), a.mean(), c.mean(), j.mean()))
for name, v in (('pose error', a), ('corr', c), ('junk', j)):
    s = spread(v, am); null = np.array([spread(v, rng.permutation(am)) for _ in range(200)])
    r1, r2 = np.corrcoef(z1, v)[0, 1], np.corrcoef(z2, v)[0, 1]
    print('  %-10s between-cell std of cell means %.4f, random cells %.4f +- %.4f -> z %+6.1f; corr with latent z1 %+.3f z2 %+.3f' % (name, s, null.mean(), null.std(), (s - null.mean()) / max(null.std(), 1e-12), r1, r2))
print('cells (n, pose error, corr, junk):', ' '.join('%d:(%d,%.1f,%.3f,%.2f)' % (k, (am == k).sum(), a[am == k].mean(), c[am == k].mean(), j[am == k].mean()) for k in cells))
