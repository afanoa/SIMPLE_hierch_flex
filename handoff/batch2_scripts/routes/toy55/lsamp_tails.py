#!/usr/bin/env python3
"""lsamp_tails.py: truth diagnostic of a score space (lsamp_score2 output): per generalised-eigen axis, the closed share and the top conformers
in the particles beyond the 2 % / 98 % quantiles. usage: lsamp_tails.py SCORES GT [EXCLUDE_ROWS]"""
import sys, numpy as np
from scipy.linalg import eigh
d = sys.argv[1]; M = np.load(d + '/meta.npz'); S = np.load(d + '/scores_pA.npy') + np.load(d + '/scores_pB.npy'); own = M['own']; eo = M['eo']
rows = np.where((own >= 0) & np.isfinite(S).all(1))[0]
if len(sys.argv) > 3: drop = np.loadtxt(sys.argv[3]) > 0; rows = rows[~drop[rows]]
lam_i = M['ppnoise'][rows]; X = S[rows] / lam_i[:, None]; X -= X.mean(1, keepdims=True); K = X.shape[1]
G = (M['gramA0'] + M['gramB0'] + M['gramA1'] + M['gramB1']) / 2 * np.mean(1. / lam_i)
w, U = np.linalg.eigh(np.eye(K) - 1. / K); B = U[:, w > 0.5]; lam, W = eigh(B.T @ np.cov(X.T) @ B, B.T @ G @ B); o = np.argsort(-lam); V = B @ W[:, o]; T = X @ V
g = np.loadtxt(sys.argv[2], comments='#'); conf = np.zeros(len(own), int); clo = np.zeros(len(own)); conf[g[:, 0].astype(int) - 1] = g[:, 1]; clo[g[:, 0].astype(int) - 1] = g[:, 2]
cf, cl = conf[rows], clo[rows]; print('%d particles, closed share %.3f; best leaf (own) closed share %.3f' % (len(rows), cl.mean(), max(cl[own[rows] == u].mean() for u in np.unique(own[rows]))))
for j in range(min(12, K - 1)):
    t = T[:, j]; lo, hi = np.percentile(t, [2, 98]); out = []
    for nm, m in (('low 2 %', t < lo), ('high 2 %', t > hi)):
        cc = np.bincount(cf[m]); top = ' '.join('%d:%.2f' % (k, cc[k] / m.sum()) for k in np.argsort(-cc)[:3] if cc[k] > 0); out.append('%s closed %.3f (%s)' % (nm, cl[m].mean(), top))
    print('axis %2d lam-1 %6.2f | %s' % (j + 1, lam[j] - 1, ' | '.join(out)))
