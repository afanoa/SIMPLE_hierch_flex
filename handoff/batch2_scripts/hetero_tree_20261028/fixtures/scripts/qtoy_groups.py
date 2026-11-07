#!/usr/bin/env python3
"""qtoy split groups against the simulation truth: per group of a node's split (split/split_children.txt) its size, junk share,
the closed / A-domain / B-blob / holo shares, mean particle quality (snrq) and leading structure combos, plus the normalised mutual
information of the groups with the structure combo and with the quality quintile (genuine particles), so a group that differs by
quality and not by structure shows as a quality axis. usage: qtoy_groups.py TREE_DIR NODE [TRUTH]"""
import sys, numpy as np
from pathlib import Path
T, node = Path(sys.argv[1]), sys.argv[2]
truth = Path(sys.argv[3] if len(sys.argv) > 3 else '/mnt/beegfs/elmlund/afan/routes/toy55/simqtoy2/truth_qtoy.txt')
COLS = 'row combo junk A_domain B_blob C_closed jtype conformer s1 closed_arm holo bfac snrq'.split()
t = np.loadtxt(truth, ndmin=2); C = {c: t[:, COLS.index(c)] for c in COLS}; tix = {int(r): i for i, r in enumerate(C['row'])}
d = T / ('node_' + node); rows = np.loadtxt(d / 'rows.txt', dtype=int, ndmin=1)
ch = np.loadtxt(d / 'split' / 'split_children.txt', dtype=int, ndmin=1)
if len(ch) != len(rows): ch = ch[rows - 1]
ti = np.array([tix[r] for r in rows]); g = ch
def nmi(a, b):
    a = np.unique(a, return_inverse=True)[1]; b = np.unique(b, return_inverse=True)[1]; n = len(a)
    pab = np.zeros((a.max() + 1, b.max() + 1)); np.add.at(pab, (a, b), 1); pab /= n
    pa, pb = pab.sum(1), pab.sum(0); nz = pab > 0
    mi = (pab[nz] * np.log(pab[nz] / (pa[:, None] * pb[None, :])[nz])).sum()
    ha, hb = -(pa[pa > 0] * np.log(pa[pa > 0])).sum(), -(pb[pb > 0] * np.log(pb[pb > 0])).sum()
    return mi / max(np.sqrt(ha * hb), 1e-12)
gen = (C['junk'][ti] == 0) & (g > 0)
q = np.digitize(C['snrq'][ti], np.quantile(C['snrq'][ti][gen], [0.2, 0.4, 0.6, 0.8]))
print('node %s: %d rows, %d groups; NMI(groups, structure combo) %.3f, NMI(groups, quality quintile) %.3f (genuine particles)' %
      (node, len(rows), len(np.unique(g[g > 0])), nmi(g[gen], C['combo'][ti][gen]), nmi(g[gen], q[gen])))
print('%5s %7s %6s %7s %7s %7s %6s %8s  %s' % ('group', 'rows', 'junk', 'closed', 'A_dom', 'B_blob', 'holo', 'snrq', 'top combos (share)'))
for k in sorted(np.unique(g)):
    m = g == k; tt = ti[m]; cb, cn = np.unique(C['combo'][tt], return_counts=True); o = np.argsort(-cn)[:3]
    print('%5d %7d %6.3f %7.3f %7.3f %7.3f %6.3f %8.3f  %s' % (k, m.sum(), C['junk'][tt].mean(), C['C_closed'][tt].mean(), C['A_domain'][tt].mean(),
          C['B_blob'][tt].mean(), C['holo'][tt].mean(), C['snrq'][tt].mean(), ' '.join('%d (%.2f)' % (cb[j], cn[j] / m.sum()) for j in o)))
