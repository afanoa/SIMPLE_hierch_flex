#!/usr/bin/env python3
"""usage: validate_conf.py NAME TREE_DIR TRUTH_STAR [--table factorial/table.txt] [--rep TREE_DIR2]
validate_tree.py's columns for a tree on a two-conformation set with STAR truth (CAHRA 13991: conformation 1 = open, no junk):
rows in no leaf (all genuine, as the set has no junk); leaf purity = size-weighted majority-conformation share of the leaf members
(the junk-inclusive purity is the same number); the tree's most open leaf (precision / recall); the open-majority leaves pooled
(precision / recall, binary accuracy over all rows, NMI(leaf; conformation) over the leaf members); the per-site layer's best open
pole (labels precision / recall, posterior AUC) from TABLE; with --rep, the replicate NMI of the two trees' leaves and the Jaccard
of their open-majority pools."""
import os, sys, numpy as np
sys.path.insert(0, '/mnt/beegfs/elmlund/afan/cahra_13991'); import cahra_score as cs
NAME, D, GT = sys.argv[1:4]; opt = dict(zip(sys.argv[4::2], sys.argv[5::2])); TAB = opt.get('--table', 'factorial/table.txt')
op = np.array(cs.read_star(GT)['conformation'], int) == 1; NM = len(op)
def tree(d):
    L = np.loadtxt(d + '/hetero_tree_labels.txt', dtype=int)[:NM]; ks = [k for k in np.unique(L) if k > 0]
    pool = np.isin(L, [k for k in ks if op[L == k].mean() > 0.5]); return L, ks, pool
L, leaves, pool = tree(D); inl = L > 0
names = [l.split()[0] for l in open(D + '/hetero_tree_leaves.txt') if l.strip()]
def auc(s, y): rk = np.argsort(np.argsort(s)) + 1.; n1 = y.sum(); n0 = len(y) - n1; return (rk[y].sum() - n1 * (n1 + 1) / 2) / max(n1 * n0, 1)
print('== %s: %s, %d leaves, %d rows in no leaf (truth open %d / closed %d)' % (NAME, D, len(leaves), (~inl).sum(), op.sum(), (~op).sum()))
lf = {k: (op[L == k].mean(), (op & (L == k)).sum() / op.sum(), (L == k).sum()) for k in leaves}
for k in leaves: print('   leaf %2d %-14s %6d, open %.3f (recall %.3f)' % (k, names[k - 1] if k <= len(names) else '?', lf[k][2], lf[k][0], lf[k][1]))
pur = sum(max(op[L == k].sum(), (~op[L == k]).sum()) for k in leaves) / max(inl.sum(), 1); kb = max(lf, key=lambda k: lf[k][0])
pp, pr = (pool & op).sum() / max(pool.sum(), 1), (pool & op).sum() / op.sum(); acc = ((pool & op) | (~pool & ~op)).mean()
best = None; p = TAB if os.path.isabs(TAB) else os.path.join(D, TAB)
for w in ([l.split() for l in open(p) if l.strip() and not l.startswith('#')] if os.path.exists(p) else []):
    lab = np.loadtxt(w[2], dtype=int); Q = np.loadtxt(w[3], ndmin=2); Q = Q if Q.shape[0] == NM else Q.T; r = np.loadtxt(w[4], dtype=int).ravel()
    rows = (r > 0) if (len(r) == NM and r.max() <= 1) else np.isin(np.arange(NM), r - 1)
    for j in range(Q.shape[1]):
        m = (lab == j + 1) & rows; P = Q[:, j]; ok = rows & (P >= 0)
        if not m.any() or op[ok].sum() == 0: continue
        a_ = auc(P[ok], op[ok]); c_ = ((m & op).sum() / m.sum(), (m & op).sum() / (rows & op).sum(), max(a_, 1 - a_), '%s pole %d' % (w[0], j + 1), m.sum())
        if best is None or c_[0] > best[0]: best = c_
site = 'site layer best open pole %s: %d labelled, precision %.3f recall %.3f, posterior AUC %.3f' % (best[3], best[4], best[0], best[1], best[2]) if best else 'site layer: no table'
rep = ''
if '--rep' in opt:
    L2, _, pool2 = tree(opt['--rep']); b = inl & (L2 > 0)
    rep = ' | replicate NMI %.3f, open-pool Jaccard %.3f' % (cs.nmi(list(L[b]), list(L2[b])), (pool & pool2).sum() / max((pool | pool2).sum(), 1))
print('SUMMARY %-12s junk n/a (no junk; %d rows in no leaf) | junk in leaves n/a | %d leaves, purity %.3f (incl. junk %.3f) | A-pure share n/a | tree open leaf %.3f/%.3f | open-majority pool %.3f/%.3f, accuracy %.3f, NMI %.3f | %s%s'
      % (NAME, (~inl).sum(), len(leaves), pur, pur, lf[kb][0], lf[kb][1], pp, pr, acc, cs.nmi(list(L[inl]), list(op[inl].astype(int))), site, rep))
