#!/usr/bin/env python3
"""usage: pick_igm.py OUTDIR TREE_DIR TRUTH
For each arm angle of the IgM truth, the population-wide kept site posterior that tracks it best: every posterior column of the
kept sites (TREE_DIR/factorial/table.txt) whose level covers >= 90 % of the rows, by abs(Spearman) over the genuine level rows.
Writes OUTDIR/<site>_p<column>.txt (one value per row, < 0 outside the level, oriented so that a high score means a high angle)
once per picked column, and prints one line per angle: 'PICK angle site column rho file', with ' same' when an earlier angle picked
the column, then the next-best column of each angle."""
import sys, numpy as np
from scipy.stats import spearmanr
OUT, D, TF = sys.argv[1:4]
hd = [l.split()[1:] for l in open(TF) if l.startswith('#')][0]; T = np.genfromtxt(TF, comments='#'); tr = T[:, 0].astype(int) - 1
L = np.loadtxt(D + '/hetero_tree_labels.txt', dtype=int); NM = len(L); assert tr.max() + 1 == NM, 'truth rows %d != tree rows %d' % (tr.max() + 1, NM)
def col(c): v = np.full(NM, np.nan); v[tr] = T[:, hd.index(c)]; return v
junk = col('junk') > 0; ang = hd[3:]; cands = []
for w in [l.split() for l in open(D + '/factorial/table.txt') if l.strip() and not l.startswith('#')]:
    Q = np.loadtxt(w[3], ndmin=2); Q = Q if Q.shape[0] == NM else Q.T
    for j in range(Q.shape[1]):
        if (Q[:, j] >= 0).mean() >= 0.9: cands.append((w[0], j + 1, Q[:, j].copy()))
print('%d population-wide posterior columns over %d rows' % (len(cands), NM)); done = {}
for a in ang:
    x = col(a); rk = []
    for s, j, p in cands:
        m = (p >= 0) & ~junk & ~np.isnan(x); rk.append((abs(spearmanr(p[m], x[m])[0]), spearmanr(p[m], x[m])[0], s, j, p))
    rk.sort(key=lambda t: -t[0]); _, r, s, j, p = rk[0]; f = '%s/%s_p%d.txt' % (OUT, s, j); same = (s, j) in done
    if not same:
        q = p.copy(); m = q >= 0
        if r < 0: q[m] = 1. - q[m]
        np.savetxt(f, q, fmt='%.6f'); done[(s, j)] = a
    print('PICK %s %s %d %+.3f %s%s' % (a, s, j, r, f, ' same' if same else ''))
    print('   next %s p%d %+.3f, %s p%d %+.3f' % (rk[1][2], rk[1][3], rk[1][1], rk[2][2], rk[2][3], rk[2][1]))
