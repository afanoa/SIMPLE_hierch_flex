#!/usr/bin/env python3
"""usage: track_node.py REF_TREE REF_NODE TREE TRUTH
Where the rows of REF_TREE's node REF_NODE (the leaves on paths through it) went in TREE: the three TREE nodes or leaves that overlap
them best (Jaccard). For the best one, its split: each child's size, share of the reference rows, and truth shares (junk, A, B, C
among genuine). A K >= 2 split whose children differ by < 0.1 in every share is the same null-type split; a leaf means no split."""
import sys, numpy as np
RD, RN, D, TF = sys.argv[1:5]
hd = [l.split()[1:] for l in open(TF) if l.startswith('#')][0]; T = np.genfromtxt(TF, comments='#'); tr = T[:, 0].astype(int) - 1
def path(nm):
    node, _, rest = nm.partition('_'); d = node[1:]; p = ['r' + d[:i] for i in range(len(d) + 1)]
    return p + [nm] if rest.startswith('s') else p
def tree(d):
    L = np.loadtxt(d + '/hetero_tree_labels.txt', dtype=int); names = [l.split()[0] for l in open(d + '/hetero_tree_leaves.txt') if l.strip()]
    sub, kids = {}, {}
    for k in [k for k in np.unique(L) if k > 0]:
        p = path(names[k - 1])
        for i, x in enumerate(p):
            sub.setdefault(x, []).append(k)
            if i + 1 < len(p): kids.setdefault(x, set()).add(p[i + 1])
    return L, sub, kids
Lr, subr, _ = tree(RD); L, sub, kids = tree(D); NM = len(L)
def col(c): v = np.full(NM, np.nan); v[tr] = T[:, hd.index(c)]; return v
junk = col('junk') > 0; fac = {'junk': junk, 'A': (col('A_domain') > 0) & ~junk, 'B': (col('B_blob') > 0) & ~junk, 'C': (col('C_closed') > 0) & ~junk}
ref = np.isin(Lr, subr[RN]); print('== %s %s: %d rows (junk %.3f, A %.3f, C %.3f among genuine) -> %s' % (RD, RN, ref.sum(), junk[ref].mean(), fac['A'][ref & ~junk].mean(), fac['C'][ref & ~junk].mean(), D))
jac = {x: (np.isin(L, ks) & ref).sum() / max((np.isin(L, ks) | ref).sum(), 1) for x, ks in sub.items()}
for x in sorted(jac, key=jac.get, reverse=True)[:3]:
    m = np.isin(L, sub[x]); print('   %-12s Jaccard %.3f: %d rows, %.3f of them reference rows, %.3f of the reference' % (x, jac[x], m.sum(), (m & ref).sum() / m.sum(), (m & ref).sum() / ref.sum()))
x = 'r'; placed = (L > 0) & ref; print('   %.3f of the reference rows are in leaves' % (placed.sum() / ref.sum()))
while True:
    ch = sorted(kids.get(x, [])); rx = placed & np.isin(L, sub[x]); shr = {c: (rx & np.isin(L, sub[c])).sum() / max(rx.sum(), 1) for c in ch}
    if len(ch) >= 2 and max(shr.values()) < 0.8: break
    if not ch: print('VERDICT %s: the reference rows stay together (>= 80 %%) down to leaf %s -> the null split is gone' % (RN, x)); sys.exit()
    x = max(shr, key=shr.get)
parts = {c: rx & np.isin(L, sub[c]) for c in ch if shr[c] >= 0.05}; sh = {}
for c, m in parts.items():
    g = m & ~junk; sh[c] = {f: (v[m] if f == 'junk' else v[g]).mean() for f, v in fac.items()}
    print('   %s child %-12s %.3f of the reference rows (%d) | %s' % (x, c, shr[c], m.sum(), ' '.join('%s %.3f' % (f, sh[c][f]) for f in fac)))
big = [c for c in parts if parts[c].sum() >= 200]; dmax = max(np.ptp([sh[c][f] for c in big]) for f in fac) if len(big) > 1 else 0.
print('VERDICT %s: node %s divides the reference rows K%d, max truth share difference between the parts %.3f -> %s' % (RN, x, len(parts), dmax,
      'the same null-type split' if dmax < 0.1 else 'a division with truth support' if dmax >= 0.2 else 'a weak division'))
