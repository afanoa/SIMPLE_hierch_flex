#!/usr/bin/env python3
"""qflex_leaves.py: the positive-control toy's final leaves against its truth: three discrete kinds (arm = the 30-map s1 continuum,
closed, holo) plus junk. Per leaf its size, kind shares and the arm particles' s1 (mean, p10-p90); the NMI of the leaves with the kind
(genuine particles); the best leaf's F1 for closed and for holo; and how the tree cut the arm: the leaves holding >= 5 % of the arm
particles, ordered by s1, with their s1 ranges (a cut along a continuum gives contiguous, overlapping ranges).
usage: qflex_leaves.py TREE_DIR [TRUTH]"""
import sys
import numpy as np

d = sys.argv[1]
truth = sys.argv[2] if len(sys.argv) > 2 else '/mnt/beegfs/elmlund/afan/routes/toy55/simqtoyflex/truth_qtoy.txt'
COLS = 'row combo junk A_domain B_blob C_closed jtype conformer s1 closed_arm holo bfac snrq'.split()
t = np.loadtxt(truth, ndmin=2)
C = {c: t[:, COLS.index(c)] for c in COLS}
rows = C['row'].astype(int)
kind = np.where(C['junk'] > 0, 0, np.where(C['C_closed'] > 0, 2, np.where(C['holo'] > 0, 3, 1)))
NAMES = {1: 'arm', 2: 'closed', 3: 'holo'}
lab = np.loadtxt(d + '/hetero_tree_labels.txt', dtype=int, ndmin=1)[rows - 1]


def nmi(a, b):
    a = np.unique(a, return_inverse=True)[1]
    b = np.unique(b, return_inverse=True)[1]
    P = np.zeros((a.max() + 1, b.max() + 1))
    np.add.at(P, (a, b), 1)
    P /= len(a)
    pa, pb = P.sum(1), P.sum(0)
    nz = P > 0
    mi = (P[nz] * np.log(P[nz] / np.outer(pa, pb)[nz])).sum()
    ha = -(pa[pa > 0] * np.log(pa[pa > 0])).sum()
    hb = -(pb[pb > 0] * np.log(pb[pb > 0])).sum()
    return mi / max(np.sqrt(ha * hb), 1e-12)


gen = (kind > 0) & (lab > 0)
ul = np.unique(lab[lab > 0])
arm = gen & (kind == 1)
print('== %s: %d leaves over %d particles; NMI(leaves, kind) %.3f; NMI(leaves, s1 decile | arm) %.3f'
      % (d, len(ul), (lab > 0).sum(), nmi(lab[gen], kind[gen]),
         nmi(lab[arm], np.digitize(C['s1'][arm], np.quantile(C['s1'][arm], np.linspace(.1, .9, 9))))))
print('%5s %7s %6s %6s %7s %6s   %s' % ('leaf', 'rows', 'junk', 'arm', 'closed', 'holo', 'arm s1 mean [p10, p90]'))
for k in sorted(ul, key=lambda k: -(lab == k).sum()):
    m = lab == k
    ma = m & arm
    s1 = ('%6.2f [%6.2f, %6.2f]' % (C['s1'][ma].mean(), *np.percentile(C['s1'][ma], [10, 90]))) if ma.sum() > 20 else '-'
    print('%5d %7d %6.3f %6.3f %7.3f %6.3f   %s' % (k, m.sum(), (kind[m] == 0).mean(), (kind[m] == 1).mean(),
                                                 (kind[m] == 2).mean(), (kind[m] == 3).mean(), s1))
for kd in (2, 3):
    mk = gen & (kind == kd)
    best = max((2 * (mk & (lab == k)).sum() / (mk.sum() + (lab == k).sum()), k) for k in ul)
    print('%s: %d particles, best leaf %d F1 %.3f (recall %.3f, purity %.3f)' % (NAMES[kd], mk.sum(), best[1], best[0],
          (mk & (lab == best[1])).sum() / max(mk.sum(), 1), (mk & (lab == best[1])).sum() / max((lab == best[1]).sum(), 1)))
big = [k for k in ul if (arm & (lab == k)).sum() >= 0.05 * arm.sum()]
big.sort(key=lambda k: C['s1'][arm & (lab == k)].mean())
print('arm cut into %d leaves (>= 5 %% of the arm): %s' % (len(big), ', '.join(
    'l%d s1 [%.1f, %.1f] (%d)' % (k, *np.percentile(C['s1'][arm & (lab == k)], [10, 90]), (arm & (lab == k)).sum()) for k in big)))
