#!/usr/bin/env python3
"""Cross-seed reproducibility of two hetero_tree runs on the same input rows: NMI of the leaves, of the units and of the site products
(site_product_readout.py labels_out), and per root-level site of run A the best |posterior correlation| with any root-level site of run B.
usage: seed_repro.py RUN_A RUN_B PRODUCT_A PRODUCT_B"""
import sys, os, numpy as np
A, B, pa, pb = sys.argv[1:5]
def nmi(a, b):
    a = np.unique(a, return_inverse=True)[1]; b = np.unique(b, return_inverse=True)[1]; p = np.zeros((a.max() + 1, b.max() + 1)); np.add.at(p, (a, b), 1); p /= p.sum()
    x, y = p.sum(1), p.sum(0); nz = p > 0; mi = (p[nz] * np.log(p[nz] / np.outer(x, y)[nz])).sum()
    return mi / max(np.sqrt(-(x[x > 0] * np.log(x[x > 0])).sum() * -(y[y > 0] * np.log(y[y > 0])).sum()), 1e-12)
def ld(R, f): return np.loadtxt(R + '/1_hetero_tree/' + f, dtype=int, ndmin=1)
la, lb = ld(A, 'hetero_tree_labels.txt'), ld(B, 'hetero_tree_labels.txt'); ok = (la > 0) & (lb > 0)
print('leaves %d vs %d: NMI %.3f' % (la.max(), lb.max(), nmi(la[ok], lb[ok])))
if os.path.exists(A + '/1_hetero_tree/hetero_tree_merged.txt'):
    ua, ub = ld(A, 'hetero_tree_merged.txt'), ld(B, 'hetero_tree_merged.txt'); ok2 = (ua > 0) & (ub > 0); print('units %d vs %d: NMI %.3f' % (ua.max(), ub.max(), nmi(ua[ok2], ub[ok2])))
qa, qb = np.loadtxt(pa, dtype=int), np.loadtxt(pb, dtype=int); print('site products %d vs %d combinations: NMI %.3f' % (len(np.unique(qa)), len(np.unique(qb)), nmi(qa, qb)))
def sites(R):
    T = R + '/1_hetero_tree'; out = []
    for l in open(T + '/factorial/table.txt'):
        if l.startswith('#') or l.split()[1] != 'r': continue
        f = l.split()[3]; f = f if os.path.exists(f) else T + '/factorial/post/' + os.path.basename(f); out.append((l.split()[0], np.loadtxt(f, ndmin=2)[:, 0]))
    return out
sa, sb = sites(A), sites(B)
for n1, p1 in sa:
    best = max((abs(np.corrcoef(p1, p2)[0, 1]), n2) for n2, p2 in sb)
    print('  site %-6s best match %-6s |r| %.2f' % (n1, best[1], best[0]))
