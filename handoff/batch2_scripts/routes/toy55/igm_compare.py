#!/usr/bin/env python3
"""IgM arms: per arm angle, how well each representation orders the particles (genuine rows):
 leaves  corr(leaf mean of the angle, angle) = sqrt(eta^2), the best any leaf-constant score can do;
 factors |corr| of the particle's leaf code on the best leaf-map component (still leaf-constant: a re-assignment would go beyond it);
 sites   |Spearman| of the best root-level site pole posterior (site layer, every particle scored), and that site's name.
usage: igm_compare.py TREE_DIR LF_DIR TRUTH"""
import sys, os, numpy as np
from scipy.stats import spearmanr
T, lf, tr = sys.argv[1:4]
cols = open(tr).readline().lstrip('#').split(); t = np.loadtxt(tr, comments='#', ndmin=2); C = {c: t[:, i] for i, c in enumerate(cols)}
rows = C['row'].astype(int); gen = C['junk'] == 0; ang = [c for c in cols if c.startswith('a')]
lab = np.loadtxt(T + '/hetero_tree_labels.txt', dtype=int, ndmin=1)[rows - 1]
leafnames = [l.split()[0] for l in open(T + '/hetero_tree_leaves.txt') if l.strip()]
L = [l.split() for l in open(lf + '/codes.txt') if not l.startswith('#')]; r = (len(L[0]) - 2) // 2
code = {leafnames.index(l[0]) + 1: np.array([float(x) for x in l[2:2 + r]]) for l in L}
g = gen & (lab > 0); P = np.array([code[k] for k in lab[g]])
sites = []
for line in open(T + '/factorial/table.txt'):
    if line.startswith('#'): continue
    nm, lev, lfp, pfp = line.split()[:4]
    if lev != 'r': continue
    pp = np.loadtxt(os.path.join(T, 'factorial', 'post', os.path.basename(pfp)), ndmin=2)
    sites.append((nm, pp[rows - 1]))
print('%-6s %7s %7s %7s  %s' % ('angle', 'leaves', 'factors', 'sites', 'best site (pole)'))
for a in ang:
    y = C[a][g]; u, inv = np.unique(lab[g], return_inverse=True); mu = np.bincount(inv, weights=y) / np.bincount(inv); lv = np.corrcoef(mu[inv], y)[0, 1]
    fc = max(abs(np.corrcoef(P[:, j], y)[0, 1]) for j in range(r))
    best = (0, '')
    for nm, pp in sites:
        for j in range(pp.shape[1]):
            v = pp[g, j]; ok = v >= 0
            if ok.sum() > 1000:
                rho = abs(spearmanr(v[ok], y[ok])[0])
                if rho > best[0]: best = (rho, '%s p%d' % (nm, j + 1))
    print('%-6s %7.3f %7.3f %7.3f  %s' % (a, lv, fc, best[0], best[1]))
