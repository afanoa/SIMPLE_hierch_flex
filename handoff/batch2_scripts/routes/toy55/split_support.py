#!/usr/bin/env python3
"""usage: split_support.py KIND TRUTH NAME1 TREE1 NAME2 TREE2
Seed recurrence of splits, like bootstrap support. Every inner node of a tree is a split whose children are particle sets (the union
of the final leaves under each child). A split X of one tree is compared with every split Y of the other tree:
- on the particles the two parents share, agreement = max over pairings of X's children with distinct children of Y of the mean
  Jaccard over X's children (an unmatched child scores 0);
- weighted by the parents' overlap |S1 & S2| / |S1|.
Support(X) = the best weighted agreement over Y. Both directions are reported, so every split of both trees gets one support.
Truth class, cut at >= 0.2 supported, >= 0.1 weak, else none (children of >= 200 rows):
- sites (simsites truth): validate_tree3.py's rule, the max share difference of junk, each junk group, and A, B, C among genuine;
- conf (a STAR truth, conformation 1 = open): the open-share difference.
Leaves: each leaf's best Jaccard with a leaf of the other tree. Prints SPLIT and LEAF rows."""
import sys, itertools, numpy as np
KIND, TF, N1, D1, N2, D2 = sys.argv[1:7]
if KIND == 'conf':
    sys.path.insert(0, '/mnt/beegfs/elmlund/afan/cahra_13991'); import cahra_score as cs
    op = np.array(cs.read_star(TF)['conformation'], int) == 1; NM = len(op); fac = {'open': op}
else:
    hd = [l.split()[1:] for l in open(TF) if l.startswith('#')][0]; T = np.genfromtxt(TF, comments='#'); tr = T[:, 0].astype(int) - 1
    NM = int(tr.max()) + 1
    def col(c): v = np.zeros(NM); v[tr] = T[:, hd.index(c)]; return v
    junk = col('junk') > 0; jt = col('jtype').astype(int); fac = {'junk': junk}
    for q in sorted(set(jt[junk])): fac['j%d' % q] = junk & (jt == q)
    for c in ('A_domain', 'B_blob', 'C_closed'): fac[c[0]] = (col(c) > 0, ~junk)
def path(nm):
    node, _, rest = nm.partition('_'); d = node[1:]; p = ['r' + d[:i] for i in range(len(d) + 1)]
    return p + [nm] if rest.startswith('s') else p
def tree(D):
    L = np.loadtxt(D + '/hetero_tree_labels.txt', dtype=int)[:NM]; names = [l.split()[0] for l in open(D + '/hetero_tree_leaves.txt') if l.strip()]
    kids = {}
    for k in [k for k in np.unique(L) if k > 0]:
        p = path(names[k - 1])
        for i, x in enumerate(p[:-1]): kids.setdefault(x, {}).setdefault(p[i + 1], []).append(k)
    spl = {x: [np.isin(L, ks) for ks in ch.values()] for x, ch in kids.items() if len(ch) >= 2}
    lvs = {names[k - 1]: L == k for k in np.unique(L) if k > 0}
    return spl, lvs
def truth(ch):
    big = [m for m in ch if m.sum() >= 200]; best = ('-', 0.)
    for f, v in fac.items():
        m_, g = (v, None) if not isinstance(v, tuple) else v
        sh = [m_[m].mean() if g is None else m_[m & g].mean() for m in big if (m.sum() if g is None else (m & g).sum()) >= 200]
        if len(sh) > 1 and np.ptp(sh) > best[1]: best = (f, float(np.ptp(sh)))
    return ('supported' if best[1] >= 0.2 else 'weak' if best[1] >= 0.1 else 'none'), best
def jac(a, b): u = (a | b).sum(); return (a & b).sum() / u if u else 0.
def pairing(J):
    K, M = J.shape
    if K <= M: return max(sum(J[i, q[i]] for i in range(K)) for q in itertools.permutations(range(M), K)) / K
    return max(sum(J[q[j], j] for j in range(M)) for q in itertools.permutations(range(K), M)) / K
def support(ch, other):
    S1 = np.any(ch, axis=0); best = (0., '-')
    for y, dh in other.items():
        S2 = np.any(dh, axis=0); sh = S1 & S2; w = sh.sum() / S1.sum()
        if w == 0: continue
        agr = pairing(np.array([[jac(c & sh, d & sh) for d in dh] for c in ch]))
        if w * agr > best[0]: best = (w * agr, y)
    return best
(s1, l1), (s2, l2) = tree(D1), tree(D2)
for nA, sA, lA, nB, sB, lB in ((N1, s1, l1, N2, s2, l2), (N2, s2, l2, N1, s1, l1)):
    for x in sorted(sA, key=lambda s: (len(s), s)):
        ch = sA[x]; tc, (f, d) = truth(ch); sp, y = support(ch, sB)
        print('SPLIT %-10s %-9s n %6d depth %d K%d %-9s (%s %.2f) support %.3f (best %s)' % (nA, x, np.any(ch, axis=0).sum(), len(x) - 1, len(ch), tc, f, d, sp, y))
    for nm, m in sorted(lA.items()):
        bj, bn = max((jac(m, mm), nn) for nn, mm in lB.items())
        print('LEAF  %-10s %-14s n %6d best Jaccard %.3f (%s)' % (nA, nm, m.sum(), bj, bn))
