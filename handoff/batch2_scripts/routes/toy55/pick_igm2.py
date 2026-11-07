#!/usr/bin/env python3
"""usage: pick_igm2.py OUTDIR TREE_DIR TRUTH [MIN_ROWS]   (pick_igm.py with each site at its own level)
For each arm angle of the IgM truth, the kept site posterior column (TREE_DIR/factorial/table.txt; a level of >= MIN_ROWS rows,
default 20000) that tracks it best by abs(Spearman) over the genuine rows of its level. Writes OUTDIR/<site>_p<column>.txt (one value
per row, < 0 outside the level, high = high angle) once per picked column and prints one line per angle, 'PICK angle site column rho
file level rows' (' same' when an earlier angle picked the column), then the next two candidates."""
import sys, numpy as np
from scipy.stats import spearmanr
OUT, D, TF = sys.argv[1:4]; MINR = int(sys.argv[4]) if len(sys.argv) > 4 else 20000
hd = [l.split()[1:] for l in open(TF) if l.startswith('#')][0]; T = np.genfromtxt(TF, comments='#'); tr = T[:, 0].astype(int) - 1
L = np.loadtxt(D + '/hetero_tree_labels.txt', dtype=int); NM = len(L); assert tr.max() + 1 == NM, 'truth rows %d != tree rows %d' % (tr.max() + 1, NM)
def col(c): v = np.full(NM, np.nan); v[tr] = T[:, hd.index(c)]; return v
junk = col('junk') > 0; ang = hd[3:]; cands = []
for w in [l.split() for l in open(D + '/factorial/table.txt') if l.strip() and not l.startswith('#')]:
    Q = np.loadtxt(w[3], ndmin=2); Q = Q if Q.shape[0] == NM else Q.T
    for j in range(Q.shape[1]):
        if (Q[:, j] >= 0).sum() >= MINR: cands.append((w[0], j + 1, w[1], Q[:, j].copy()))
print('%d posterior columns on levels of >= %d rows (of %d)' % (len(cands), MINR, NM)); done = set()
for a in ang:
    x = col(a); rk = []
    for s, j, lv, p in cands:
        m = (p >= 0) & ~junk & ~np.isnan(x); r = spearmanr(p[m], x[m])[0]; rk.append((abs(r), r, s, j, lv, p))
    rk.sort(key=lambda t: -t[0]); _, r, s, j, lv, p = rk[0]; f = '%s/%s_p%d.txt' % (OUT, s, j); same = (s, j) in done
    if not same:
        q = p.copy(); m = q >= 0
        if r < 0: q[m] = 1. - q[m]
        np.savetxt(f, q, fmt='%.6f'); done.add((s, j))
    print('PICK %s %s %d %+.3f %s %s %d%s' % (a, s, j, r, f, lv, (p >= 0).sum(), ' same' if same else ''))
    print('   next %s p%d %+.3f (%s), %s p%d %+.3f (%s)' % (rk[1][2], rk[1][3], rk[1][1], rk[1][4], rk[2][2], rk[2][3], rk[2][1], rk[2][4]))
