#!/usr/bin/env python3
"""usage: score_port_sites.py (--out DIR | --table SITES) --truth T [--tcols C_closed,A_domain --tree_labels F --tops 2500,10000]
--table: one site per line, 'name level labels post rows' (per master row: labels 0 / 1 / 2 = neither / pole 1 / pole 2; post = one column
(P of pole 2) or two (P of pole 1, P of pole 2), < 0 outside the level; rows = 0/1 per master row, or the node's 1-based master-row list).
Per kept site of a tree_sites.py run (OUT/sites_final.txt) and each of its poles (lo / hi at the pooled level), against truth: particles,
junk share, and per binary truth column the precision (share of the pole's genuine particles in that class), recall (share of the
level's genuine particles of that class in the pole), the pole posterior's AUC and top-N purity. Reference: the tree's own leaves and its
junk removal (tree_labels: leaf per master row, 0 = removed)."""
import argparse, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('--out', default=''); ap.add_argument('--truth', required=True); ap.add_argument('--tcols', default='C_closed,A_domain')
ap.add_argument('--tree_labels', default=''); ap.add_argument('--tops', default='2500,10000'); ap.add_argument('--table', default=''); a = ap.parse_args()
hd = [l.split()[1:] for l in open(a.truth) if l.startswith('#')][0]; T = np.genfromtxt(a.truth, comments='#'); tr = T[:, 0].astype(int) - 1; NM = tr.max() + 1
def col(c): v = np.zeros(NM); v[tr] = T[:, hd.index(c)]; return v
junk = col('junk') > 0; TC = {c: (col(c) > 0) & ~junk for c in a.tcols.split(',') if c}; tops = [int(x) for x in a.tops.split(',') if x]
def auc(s, y): rk = np.argsort(np.argsort(s)) + 1.; n1 = y.sum(); n0 = len(y) - n1; return (rk[y].sum() - n1 * (n1 + 1) / 2) / max(n1 * n0, 1)
def cls(m, rows, tag):
    g = rows & ~junk; s = '  %s: %d particles, junk %.3f' % (tag, m.sum(), (m & junk).sum() / max(m.sum(), 1))
    for c, y in TC.items(): s += ' | %s precision %.3f recall %.3f (base %.3f)' % (c, (m & y).sum() / max((m & ~junk).sum(), 1), (m & y).sum() / max((g & y).sum(), 1), (g & y).sum() / max(g.sum(), 1))
    return s
if a.tree_labels:
    TL = np.loadtxt(a.tree_labels, dtype=int); allr = np.ones(NM, bool); rem = TL == 0
    print('tree: junk removal precision %.3f recall %.3f (%d removed); leaves:' % ((rem & junk).sum() / max(rem.sum(), 1), (rem & junk).sum() / max(junk.sum(), 1), rem.sum()))
    for l in np.unique(TL[TL > 0]): print(cls(TL == l, allr, 'leaf %d' % l))
def rowmask(f):
    r = np.loadtxt(f, dtype=int).ravel()
    if len(r) == NM and r.max() <= 1: return r > 0
    m = np.zeros(NM, bool); m[r - 1] = True; return m
if a.table:
    SITES = []
    for name, lev, labf, postf, rowf in [l.split()[:5] for l in open(a.table) if l.strip() and not l.startswith('#')]:
        Q = np.loadtxt(postf, ndmin=2); Q = Q.T if Q.shape[0] != NM else Q
        SITES.append((name, lev, '-', '-', np.loadtxt(labf, dtype=int), rowmask(rowf), {1: Q[:, 0], 2: Q[:, -1]} if Q.shape[1] == 2 else {2: Q[:, 0]}))
else:
    SITES = [(name, lev, home + ', port candidate ' + rep, repro, np.loadtxt(labf, dtype=int), np.loadtxt('%s/view/rows_%s.txt' % (a.out, lev)) > 0,
              {k: np.loadtxt('%s/fits/%s@%s/%s/post_master.txt' % (a.out, name, lev, pole)) for k, pole in ((1, 'lo'), (2, 'hi'))})
             for name, home, lev, rep, repro, nlo, nhi, labf, mask in [l.split() for l in open(a.out + '/sites_final.txt') if l.strip() and not l.startswith('#')]]
for name, lev, home, repro, L, rows, POST in SITES:
    print('site %s (home %s, level %s, repro %s): %d particles at the level, junk %.3f' % (name, home, lev, repro, rows.sum(), (rows & junk).sum() / max(rows.sum(), 1)))
    for k, pole in ((1, 'lo'), (2, 'hi')):
        s = cls(L == k, rows, '%s pole' % pole)
        if k not in POST: print(s); continue
        P = POST[k]; ok = rows & ~junk & np.isfinite(P) & (P >= 0)
        for c, y in TC.items():
            o = np.argsort(-np.where(ok, P, -np.inf)); s += '\n      %s by the %s posterior: AUC %.3f, top-N purity %s' % (c, pole, auc(P[ok], y[ok]), ' '.join('%d: %.3f' % (N, y[o[:N]].mean()) for N in tops))
        print(s)
    if a.tree_labels:
        for l in np.unique(TL[TL > 0]):
            for k in (1, 2): print(cls((L == k) & (TL == l), rows, 'leaf %d x %s pole' % (l, ('lo', 'hi')[k - 1])))
