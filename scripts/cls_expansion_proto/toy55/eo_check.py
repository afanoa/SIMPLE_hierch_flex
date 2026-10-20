#!/usr/bin/env python3
"""Is the even/odd split of a refined project exchangeable? Per family table: eo=1 fraction by conformer, in the first
vs second half of the rows, and by junk. usage: eo_check.py GT_PARTICLES TABLE[,TABLE] SUBSET[,SUBSET]"""
import sys, numpy as np
G = np.loadtxt(sys.argv[1], comments='#'); conf = G[:, 1].astype(int)
for tab, sub in zip(sys.argv[2].split(','), sys.argv[3].split(',')):
    eo = np.array([float(dict(t.split('=', 1) for t in l.split() if '=' in t).get('eo', -1)) for l in open(tab) if l.strip()])
    rows = np.loadtxt(sub).astype(int) - 1; rows = rows[rows < len(eo)]; e = eo[rows]; c = conf[rows]; h = rows < np.median(rows)
    print('%s: %d rows, eo values %s; eo=1 fraction: first half of rows %.3f, second half %.3f; by conformer %s' % (tab.split('/')[-2], len(rows), np.unique(e), e[h].mean(), e[~h].mean(), ' '.join('%d:%.3f(n=%d)' % (k, e[c == k].mean(), (c == k).sum()) for k in np.unique(c))))
