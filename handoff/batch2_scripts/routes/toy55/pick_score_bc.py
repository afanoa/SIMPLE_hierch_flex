#!/usr/bin/env python3
"""usage: pick_score_bc.py OUT TARGET KIND TRUTH TREE_DIR   (pick_score.py on brightness-corrected candidates)
pick_score.py's candidates (the population-wide kept sites' pole-1 posteriors and the root flex / fflex z) with the per-particle
contrast regressed out: a posterior as its logit, a z as is, is replaced by its residual on [1, c, c^2] over the level rows, with
c the root's posterior-mean contrast (TREE_DIR/node_r/site/site_contrast.txt). Picks by abs(Spearman) with the truth over the
genuine level rows as pick_score.py does, and prints the truth's own Spearman with c and, per candidate, the Spearman with the
truth and with c before and after. Writes OUT as percentile ranks of the residual (high = more TARGET; < 0 outside the level)."""
import sys, numpy as np
from scipy.stats import spearmanr
OUT, TGT, KIND, TF, D = sys.argv[1:6]
L = np.loadtxt(D + '/hetero_tree_labels.txt', dtype=int); NM = len(L)
if KIND == 'cahra':
    sys.path.insert(0, '/mnt/beegfs/elmlund/afan/cahra_13991'); import cahra_score as cs
    y = (np.array(cs.read_star(TF)['conformation'], int)[:NM] == 1).astype(float); G = np.ones(NM, bool)
else:
    hd = [l.split()[1:] for l in open(TF) if l.startswith('#')][0]; M = np.genfromtxt(TF, comments='#'); tr = M[:, 0].astype(int) - 1
    def col(c): v = np.zeros(NM); v[tr] = M[:, hd.index(c)]; return v
    G = ~(col('junk') > 0); y = (col({'A': 'A_domain', 'C': 'C_closed'}[TGT]) > 0).astype(float)
C = np.loadtxt(D + '/node_r/site/site_contrast.txt', comments='#', ndmin=2); c = np.full(NM, np.nan); c[C[:, 0].astype(int) - 1] = C[:, 1]
cands = []
for w in [l.split() for l in open(D + '/factorial/table.txt') if l.strip() and not l.startswith('#')]:
    Q = np.loadtxt(w[3], ndmin=2); Q = Q.T if Q.shape[1] > Q.shape[0] else Q; p = Q[:, 0]
    if (p >= 0).mean() >= 0.9: q = np.clip(p, 1e-6, 1 - 1e-6); cands.append(('site ' + w[0], np.where(p >= 0, np.log(q / (1 - q)), np.nan)))
for lat in ('flex', 'fflex'):
    b = open('%s/node_r/%s/flex_pca_embedding.bin' % (D, lat), 'rb').read(); ver, n, nc, box = (int(v) for v in np.frombuffer(b[8:24], '<i4'))
    pi = np.frombuffer(b[28:28 + 4 * n], '<i4').astype(int) - 1; z = np.frombuffer(b[28 + 4 * n:28 + 4 * n + 8 * n * nc], '<f8').reshape(nc, n).T
    for j in range(nc): x = np.full(NM, np.nan); x[pi] = z[:, j]; cands.append(('%s z%d' % (lat, j + 1), x))
ok = ~np.isnan(c); print('contrast over %d rows: mean %.3f sd %.3f; truth %s vs contrast Spearman %+.3f over the genuine rows' % (ok.sum(), np.nanmean(c), np.nanstd(c), TGT, spearmanr(c[ok & G], y[ok & G])[0]))
res = []
for nm, x in cands:
    m = ok & ~np.isnan(x); A = np.column_stack([np.ones(m.sum()), c[m], c[m] ** 2]); r = np.full(NM, np.nan); r[m] = x[m] - A @ np.linalg.lstsq(A, x[m], rcond=None)[0]; g = m & G
    res.append((abs(spearmanr(r[g], y[g])[0]), spearmanr(r[g], y[g])[0], nm, r, m, spearmanr(x[g], y[g])[0], spearmanr(x[m], c[m])[0], spearmanr(r[m], c[m])[0]))
res.sort(key=lambda t: -t[0])
for a, r1, nm, _, _, r0, c0, c1 in res[:10]: print('   %-14s abs rho %.3f (%+.3f; before %+.3f) | Spearman with contrast before %+.3f, after %+.3f' % (nm, a, r1, r0, c0, c1))
a, r1, nm, r, m = res[0][:5]; s = np.full(NM, -1.); s[m] = (np.argsort(np.argsort(r[m])) + 0.5) / m.sum()
if r1 < 0: s[m] = 1. - s[m]
np.savetxt(OUT, s, fmt='%.6f'); print('picked %s for %s (abs rho %.3f, %s): %s' % (nm, TGT, a, 'flipped' if r1 < 0 else 'as is', OUT))
