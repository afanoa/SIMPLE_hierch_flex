#!/usr/bin/env python3
"""qtoy final leaves against the simulation truth (the leaf form of qtoy_groups.py): per tree (hetero_tree_labels.txt, one leaf per
project row) each leaf's size, junk share, closed / A-domain / B-blob / holo shares, mean quality (snrq) and top structure combos;
NMI of the leaves with the structure combo and with the quality quintile (genuine particles); per genuine combo the best leaf's F1.
usage: qtoy_leaves.py NAME=TREE_DIR [NAME=TREE_DIR ..] [truth=TRUTH]"""
import sys, numpy as np
args = dict(a.split('=', 1) for a in sys.argv[1:])
truth = args.pop('truth', '/mnt/beegfs/elmlund/afan/routes/toy55/simqtoy2/truth_qtoy.txt')
COLS = 'row combo junk A_domain B_blob C_closed jtype conformer s1 closed_arm holo bfac snrq'.split()
t = np.loadtxt(truth, ndmin=2); C = {c: t[:, COLS.index(c)] for c in COLS}; rows = C['row'].astype(int)
def nmi(a, b):
    a = np.unique(a, return_inverse=True)[1]; b = np.unique(b, return_inverse=True)[1]; n = len(a)
    pab = np.zeros((a.max() + 1, b.max() + 1)); np.add.at(pab, (a, b), 1); pab /= n
    pa, pb = pab.sum(1), pab.sum(0); nz = pab > 0
    mi = (pab[nz] * np.log(pab[nz] / (pa[:, None] * pb[None, :])[nz])).sum()
    ha, hb = -(pa[pa > 0] * np.log(pa[pa > 0])).sum(), -(pb[pb > 0] * np.log(pb[pb > 0])).sum()
    return mi / max(np.sqrt(ha * hb), 1e-12)
for name, d in args.items():
    lab = np.loadtxt(d + '/hetero_tree_labels.txt', dtype=int, ndmin=1)[rows - 1]
    gen = (C['junk'] == 0) & (lab > 0)
    q = np.digitize(C['snrq'], np.quantile(C['snrq'][gen], [0.2, 0.4, 0.6, 0.8]))
    ul = np.unique(lab[lab > 0])
    print('== %s: %d leaves over %d particles; NMI(leaves, structure combo) %.3f, NMI(leaves, quality quintile) %.3f (genuine)' %
          (name, len(ul), (lab > 0).sum(), nmi(lab[gen], C['combo'][gen]), nmi(lab[gen], q[gen])))
    print('%5s %7s %6s %7s %7s %7s %6s %7s  %s' % ('leaf', 'rows', 'junk', 'closed', 'A_dom', 'B_blob', 'holo', 'snrq', 'top combos (share)'))
    for k in sorted(ul, key=lambda k: -(lab == k).sum()):
        m = lab == k; cb, cn = np.unique(C['combo'][m], return_counts=True); o = np.argsort(-cn)[:3]
        print('%5d %7d %6.3f %7.3f %7.3f %7.3f %6.3f %7.3f  %s' % (k, m.sum(), C['junk'][m].mean(), C['C_closed'][m].mean(),
              C['A_domain'][m].mean(), C['B_blob'][m].mean(), C['holo'][m].mean(), C['snrq'][m].mean(),
              ' '.join('%d (%.2f)' % (cb[j], cn[j] / m.sum()) for j in o)))
    out = []
    for c in np.unique(C['combo'][C['junk'] == 0]):
        mc = (C['combo'] == c) & (C['junk'] == 0); best = (0., 0)
        for k in ul:
            ml = lab == k; tp = (mc & ml).sum()
            if tp: f = 2 * tp / (mc.sum() + ml.sum()); best = max(best, (f, k))
        out.append('%d:%.2f(l%d)' % (c, best[0], best[1]))
    print('best-leaf F1 per combo: ' + ' '.join(out))
