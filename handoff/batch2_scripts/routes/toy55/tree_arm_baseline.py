#!/usr/bin/env python3
"""usage: tree_arm_baseline.py TRUTH LABELS: a tree's leaves (per master row, 0 = removed) against continuous truth (arm angles):
eta^2 of every truth column over the leaves on genuine particles, and the tree's junk removal precision / recall."""
import sys, numpy as np
TF, LF = sys.argv[1:3]; hd = [l.split()[1:] for l in open(TF) if l.startswith('#')][0]; T = np.genfromtxt(TF, comments='#'); tr = T[:, 0].astype(int) - 1; NM = tr.max() + 1
L = np.loadtxt(LF, dtype=int); junk = np.zeros(NM, bool); junk[tr] = T[:, hd.index('junk')] > 0; rem = L == 0
print('tree: %d leaves; junk removal precision %.3f recall %.3f (%d removed)' % (len(np.unique(L[L > 0])), (rem & junk).sum() / max(rem.sum(), 1), (rem & junk).sum() / max(junk.sum(), 1), rem.sum()))
g = (L > 0) & ~junk
def eta2(lab, y): mu = y.mean(); return sum(((lab == k).sum() * (y[lab == k].mean() - mu) ** 2) for k in np.unique(lab)) / max(((y - mu) ** 2).sum(), 1e-30)
cols = [c for c in hd[3:] if c not in ('junk',)]
out = []
for c in cols:
    y = np.full(NM, np.nan); y[tr] = T[:, hd.index(c)]; m = g & np.isfinite(y); out.append('%s %.3f' % (c, eta2(L[m], y[m])))
print('eta^2 of each arm angle over the tree leaves (genuine particles): ' + ', '.join(out))
