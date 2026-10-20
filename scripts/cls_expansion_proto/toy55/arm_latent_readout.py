#!/usr/bin/env python3
"""How much of each arm's angles a klem board's 2-D latent explains: R^2 of a quadratic fit of every arm's theta and phi
on the posterior means (z1, z2), non-junk particles only (junk = conformer 0), plus the free-component share.
usage: arm_latent_readout.py KLEM_LATENT GT_PARTICLES TRUTH [TAG]"""
import sys, numpy as np
lat, gtf, trf = sys.argv[1:4]; tag = sys.argv[4] if len(sys.argv) > 4 else lat
conf = {int(w[0]): int(w[1]) for w in (l.split() for l in open(gtf)) if w and not w[0].startswith('#')}
tr = {int(w[0]): np.array([float(x) for x in w[1:]]) for w in (l.split() for l in open(trf)) if w and not w[0].startswith('#')}
L = np.loadtxt(lat, comments='#'); pind = L[:, 0].astype(int); z1, z2, gm = L[:, 1], L[:, 2], L[:, 5]
c = np.array([conf.get(p, 0) for p in pind]); ok = (c > 0) & np.array([cc in tr for cc in c])
T = np.array([tr[cc] for cc in c[ok]]); X = np.c_[np.ones(ok.sum()), z1[ok], z2[ok], z1[ok] ** 2, z2[ok] ** 2, z1[ok] * z2[ok]]
r2 = []
for j in range(T.shape[1]):
    y = T[:, j]; b = np.linalg.lstsq(X, y, rcond=None)[0]; r2.append(1. - np.sum((y - X @ b) ** 2) / np.sum((y - y.mean()) ** 2))
r2 = np.array(r2).reshape(-1, 2)
print('%-28s %6d non-junk; R^2 theta/phi per arm: %s; best arm %d (%.3f/%.3f); grid mass < 0.5 for %.3f' % (tag, ok.sum(), ' '.join('%.2f/%.2f' % tuple(v) for v in r2), np.argmax(r2.sum(1)) + 1, *r2[np.argmax(r2.sum(1))], np.mean(gm[ok] < 0.5)))
