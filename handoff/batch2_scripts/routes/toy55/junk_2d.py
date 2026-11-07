#!/usr/bin/env python3
"""usage: junk_2d.py SCORES TRUTH [--tcol junk --jtype COL --guard COLS]
Does a 2D cut on (full-band amplitude, high-band amplitude) beat the 1D amplitude cut? Two-component full-covariance Gaussian mixture,
junk = the component with the lower full-band mean, posterior >= 0.5; precision / recall / false-junk / guards next to the 1D mixture."""
import argparse, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('scores'); ap.add_argument('truth'); ap.add_argument('--tcol', default='junk'); ap.add_argument('--jtype', default=''); ap.add_argument('--guard', default='')
a = ap.parse_args(); S = np.loadtxt(a.scores); node = np.isfinite(S[:, 1]); X = S[node][:, [1, 2]]
hd = [l.split()[1:] for l in open(a.truth) if l.startswith('#')][0]; T = np.genfromtxt(a.truth, comments='#'); tr = T[:, 0].astype(int) - 1
def col(cn): v = np.full(len(S), np.nan); v[tr] = T[:, hd.index(cn)]; return v[node]
tj = col(a.tcol) > 0
def gmm(X, iters=300):
    d = X.shape[1]; lo, hi = np.percentile(X[:, 0], [5, 50]); mu = np.array([X[np.abs(X[:, 0] - lo) < 0.05].mean(0), X[np.abs(X[:, 0] - hi) < 0.05].mean(0)])
    C = np.array([np.cov(X.T).reshape(d, d) * 0.5] * 2); pi = np.array([0.2, 0.8])
    for _ in range(iters):
        L = np.stack([pi[k] * np.exp(-0.5 * np.einsum('ni,ij,nj->n', X - mu[k], np.linalg.inv(C[k]), X - mu[k])) / np.sqrt(np.linalg.det(C[k])) for k in range(2)], 1)
        R = L / np.maximum(L.sum(1, keepdims=True), 1e-300); w = R.sum(0); pi = w / len(X); mu = (R.T @ X) / w[:, None]
        C = np.array([((X - mu[k]).T * R[:, k]) @ (X - mu[k]) / w[k] + 1e-9 * np.eye(d) for k in range(2)])
    j = int(np.argmin(mu[:, 0])); return R[:, j], mu, C, pi, j
for nm, Xs in (('1D amplitude', X[:, :1]), ('2D amplitude + high-band', X)):
    p, mu, C, pi, j = gmm(Xs); J = p >= 0.5
    print('%s: junk component mean %s weight %.3f | genuine %s | flagged %d (%.4f)' % (nm, np.round(mu[j], 2), pi[j], np.round(mu[1 - j], 2), J.sum(), J.mean()))
    print('  precision %.3f recall %.3f; false-junk on genuine %.4f' % ((J & tj).sum() / max(J.sum(), 1), (J & tj).sum() / max(tj.sum(), 1), (J & ~tj).sum() / max((~tj).sum(), 1))
          + ''.join('; recall junk type %g %.3f' % (v, (J & tj & (col(a.jtype) == v)).sum() / max((tj & (col(a.jtype) == v)).sum(), 1)) for v in (np.unique(col(a.jtype)[tj]) if a.jtype else [])))
    for g in [x for x in a.guard.split(',') if x]:
        gc = col(g); print('  guard %s: ' % g + ', '.join('= %g false-junk %.4f' % (v, (J & ~tj & (gc == v)).sum() / max((~tj & (gc == v)).sum(), 1)) for v in np.unique(gc[~tj])))
