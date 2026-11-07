#!/usr/bin/env python3
"""usage: pick_score_any.py OUT COLUMN TRUTH TREE_DIR   (pick_score.py for any column of a sites truth, e.g. B_blob)
pick_score.py's ranking of the population-wide kept sites' pole-1 posteriors (level >= 90 % of rows) and the root flex / fflex z
components (as percentile ranks) by abs(Spearman) with COLUMN > 0 over the genuine level rows. Writes OUT (one value per row,
< 0 outside the level), oriented so that a high score means COLUMN is set, and prints the ranking."""
import sys, numpy as np
from scipy.stats import spearmanr
OUT, TC, TF, D = sys.argv[1:5]
L = np.loadtxt(D + '/hetero_tree_labels.txt', dtype=int); NM = len(L)
hd = [l.split()[1:] for l in open(TF) if l.startswith('#')][0]; M = np.genfromtxt(TF, comments='#'); tr = M[:, 0].astype(int) - 1
def col(c): v = np.zeros(NM); v[tr] = M[:, hd.index(c)]; return v
G = ~(col('junk') > 0); y = (col(TC) > 0).astype(float); cands = []
for w in [l.split() for l in open(D + '/factorial/table.txt') if l.strip() and not l.startswith('#')]:
    Q = np.loadtxt(w[3], ndmin=2); Q = Q.T if Q.shape[1] > Q.shape[0] else Q; p = Q[:, 0]
    if (p >= 0).mean() >= 0.9: cands.append(('site ' + w[0], p.copy()))
for lat in ('flex', 'fflex'):
    b = open('%s/node_r/%s/flex_pca_embedding.bin' % (D, lat), 'rb').read(); ver, n, nc, box = (int(v) for v in np.frombuffer(b[8:24], '<i4'))
    pi = np.frombuffer(b[28:28 + 4 * n], '<i4').astype(int) - 1; z = np.frombuffer(b[28 + 4 * n:28 + 4 * n + 8 * n * nc], '<f8').reshape(nc, n).T
    for j in range(nc): r = np.full(NM, -1.); r[pi] = (np.argsort(np.argsort(z[:, j])) + 0.5) / n; cands.append(('%s z%d' % (lat, j + 1), r))
res = []
for nm, s in cands: m = (s >= 0) & G; rho = spearmanr(s[m], y[m])[0]; res.append((abs(rho), rho, nm, s))
res.sort(key=lambda x: -x[0])
for a, r, nm, _ in res[:8]: print('   %-14s abs rho %.3f (%+.3f)' % (nm, a, r))
a, r, nm, s = res[0]; s = s.copy(); m = s >= 0
if r < 0: s[m] = 1. - s[m]
np.savetxt(OUT, s, fmt='%.6f'); print('picked %s for %s (abs rho %.3f, %s): %s' % (nm, TC, a, 'flipped' if r < 0 else 'as is', OUT))
