#!/usr/bin/env python3
"""flexsub_cut.py: a node's stage-B states and split children against the positive-control toy's kinds (qflex_leaves.py's readout
for one node): per group the size, kind shares and the arm particles' s1 (mean, p10-p90); NMI with the kind (genuine particles) and
with the s1 decile (arm particles); the best group's F1 for closed and holo.
usage: flexsub_cut.py NODE TRUTH NAME=DIR [NAME=DIR ...]   (DIR holds b_ptcl3D.txt and split_children.txt)"""
import os, sys
import numpy as np

node, truth = sys.argv[1], sys.argv[2]
COLS = open(truth).readline().lstrip('#').split()
t = np.loadtxt(truth, ndmin=2)
rows = np.loadtxt(os.path.join(node, 'rows.txt'), dtype=int, ndmin=1)
C = {c: t[rows - 1, COLS.index(c)] for c in COLS}
kind = np.where(C['junk'] > 0, 0, np.where(C['C_closed'] > 0, 2, np.where(C['holo'] > 0, 3, 1)))
NAMES = {1: 'arm', 2: 'closed', 3: 'holo'}


def nmi(a, b):
    a = np.unique(a, return_inverse=True)[1]; b = np.unique(b, return_inverse=True)[1]
    P = np.zeros((a.max() + 1, b.max() + 1)); np.add.at(P, (a, b), 1); P /= len(a)
    pa, pb = P.sum(1), P.sum(0); nz = P > 0
    mi = (P[nz] * np.log(P[nz] / np.outer(pa, pb)[nz])).sum()
    ha = -(pa[pa > 0] * np.log(pa[pa > 0])).sum(); hb = -(pb[pb > 0] * np.log(pb[pb > 0])).sum()
    return mi / max(np.sqrt(ha * hb), 1e-12)


def states(fn):
    st = []
    for ln in open(fn):
        if 'e1=' in ln:
            d = dict(t.split('=', 1) for t in ln.split() if '=' in t)
            st.append(int(float(d.get('state', 0))))
    st = np.array(st)
    return st[rows - 1] if len(st) > len(rows) else st


def show(name, lab):
    gen = (kind > 0) & (lab > 0)
    arm = gen & (kind == 1)
    s1 = C['s1']
    print('== %s: %d groups; NMI(groups, kind) %.3f; NMI(groups, s1 decile | arm) %.3f'
          % (name, len(np.unique(lab[lab > 0])), nmi(lab[gen], kind[gen]),
             nmi(lab[arm], np.digitize(s1[arm], np.quantile(s1[arm], np.linspace(.1, .9, 9))))))
    for g in sorted(np.unique(lab[lab > 0]), key=lambda g: -(lab == g).sum()):
        m = lab == g
        ma = m & arm
        r = ('%6.2f [%6.2f, %6.2f]' % (s1[ma].mean(), *np.percentile(s1[ma], [10, 90]))) if ma.sum() > 20 else '-'
        print('   group %d: %6d  junk %.3f arm %.3f closed %.3f holo %.3f   arm s1 %s'
              % (g, m.sum(), (kind[m] == 0).mean(), (kind[m] == 1).mean(), (kind[m] == 2).mean(), (kind[m] == 3).mean(), r))
    for kd in (2, 3):
        mk = gen & (kind == kd)
        best = max((2 * (mk & (lab == g)).sum() / (mk.sum() + (lab == g).sum()), g) for g in np.unique(lab[lab > 0]))
        print('   %s: best group %d F1 %.3f (recall %.3f, purity %.3f)' % (NAMES[kd], best[1], best[0],
              (mk & (lab == best[1])).sum() / max(mk.sum(), 1), (mk & (lab == best[1])).sum() / max((lab == best[1]).sum(), 1)))


for arg in sys.argv[3:]:
    name, d = arg.split('=', 1)
    show(name + ' stage B', states(os.path.join(d, 'b_ptcl3D.txt')))
    sc = os.path.join(d, 'split_children.txt')
    if os.path.exists(sc):
        ch = np.loadtxt(sc, dtype=int, ndmin=1)
        ch = ch[rows - 1] if len(ch) > len(rows) else ch
        show(name + ' split children', ch)
