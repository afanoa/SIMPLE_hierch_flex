#!/usr/bin/env python3
"""Is the sub-class step needed to pick the motion direction? Compares, per particle, |corr| with the truth variable of: the most
split-half-reliable direction of the sub-class means (as used), the leading PCs of the bin-centred per-particle latent (no sub-classes),
and each raw latent component. usage: subcls_dircheck.py EMBED_DIR NROWS GT [--codes swing:200:230]"""
import argparse, os, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('d'); ap.add_argument('nrows', type=int); ap.add_argument('gt'); ap.add_argument('--codes', default='swing:200:230'); a = ap.parse_args()
E = np.load(os.path.join(a.d, 'subcls_embed.npz')); p, sid, sbin, Z, V = E['p'], E['sid'], E['sbin'], E['Z'], E['bincentred_V']
binn = sbin[sid]; ub, bi = np.unique(binn, return_inverse=True); nb = len(ub)
bm = np.stack([np.bincount(bi, Z[:, j], minlength=nb) for j in range(Z.shape[1])], 1) / np.bincount(bi, minlength=nb)[:, None]; Zc = Z - bm[bi]; Zc -= Zc.mean(0)
G = np.loadtxt(a.gt, comments='#', ndmin=2); conf_of = np.zeros(a.nrows + 1, int); conf_of[G[:, 0].astype(int)] = G[:, 1].astype(int); conf = conf_of[p]
nm, lo, hi = a.codes.split(':'); v = np.where((conf >= int(lo)) & (conf <= int(hi)), conf - int(lo), np.nan); dv = np.isfinite(v)
c = lambda y: abs(np.corrcoef(y[dv], v[dv])[0, 1])
w, U = np.linalg.eigh(np.cov(Zc.T)); U = U[:, ::-1]; w = w[::-1]
print('sub-class reliability direction 1: %.3f | 2: %.3f' % (c(Zc @ V[:, 0]), c(Zc @ V[:, 1])))
print('per-particle PCs of the bin-centred latent (variance share): ' + ' '.join('PC%d %.3f (%.2f)' % (k + 1, c(Zc @ U[:, k]), w[k] / w.sum()) for k in range(5)))
print('raw components: ' + ' '.join('z%d %.3f' % (k + 1, c(Z[:, k])) for k in range(Z.shape[1])))
