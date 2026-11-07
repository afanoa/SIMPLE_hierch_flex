#!/usr/bin/env python3
"""Truth-free junk flag for a split tree node (the root-child rule): a child is junk when the node's own board put most
of its particles on the board's free (off-lattice) components in BOTH halves (share above THRESHOLD, default 0.8;
real junk children sit at 0.99-1.00, a junk-free root left 0.37-0.47 of a real child there). With GT files the truth
junk share of every child is printed beside it. usage: junk_children.py NODE_A NODE_B NAME [THRESHOLD [GT_A GT_B]]
(last output line = the junk child numbers, space separated, empty when none)"""
import sys, numpy as np
na, nb, name = sys.argv[1:4]; thr = float(sys.argv[4]) if len(sys.argv) > 4 else 0.8; gts = sys.argv[5:7] if len(sys.argv) > 6 else None
fr, tj = {}, {}
for h, N in enumerate((na, nb)):
    ch = np.loadtxt(N + '/children.txt', dtype=int); cfg = dict(l.split()[:2] for l in open(N + '/board/board.txt') if l.strip()); ng = int(cfg['n1']) * int(cfg['n2'])
    L = np.loadtxt(N + '/board/klem_latent.txt', comments='#', ndmin=2); free = np.zeros(len(ch) + 1, bool); free[L[:, 0].astype(int)] = L[:, 3] > ng
    conf = np.loadtxt(gts[h], comments='#')[:, 1].astype(int)[:len(ch)] if gts else None
    for c in sorted(set(ch[ch > 0])):
        rows = np.where(ch == c)[0]; fr.setdefault(c, [np.nan, np.nan])[h] = free[rows + 1].mean()
        if conf is not None: tj.setdefault(c, [np.nan, np.nan])[h] = np.mean(conf[rows] == 0)
junk = [c for c, v in sorted(fr.items()) if min(v) > thr]
for c, v in sorted(fr.items()):
    print('   child %s.%d: on the free components A %.2f B %.2f%s%s' % (name, c, v[0], v[1], '; truth junk share A %.2f B %.2f' % tuple(tj[c]) if c in tj else '', ' -> JUNK leaf, not split' if c in junk else ''))
print(' '.join(str(c) for c in junk))
