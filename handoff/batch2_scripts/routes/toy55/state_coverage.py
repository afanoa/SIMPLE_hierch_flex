#!/usr/bin/env python3
"""Does the latent give EVERY state a home? Per true conformer: its best cell's purity, the share of its particles in
cells where it is the majority, the number of cells it owns (majority), recall at those cells; plus the Hungarian
one-cell-per-conformer assignment. Cells = rounded posterior means on the n1 x n2 grid (argmax cells optional).
usage: state_coverage.py klem_latent.txt gt_particles.txt [--n1 7 --n2 7] [--argmax] [--tag T]"""
import sys, argparse, numpy as np
from scipy.optimize import linear_sum_assignment
ap = argparse.ArgumentParser(); ap.add_argument('latent'); ap.add_argument('gt'); ap.add_argument('--n1', type=int, default=7); ap.add_argument('--n2', type=int, default=7); ap.add_argument('--argmax', action='store_true'); ap.add_argument('--tag', default='')
a = ap.parse_args(); gt = {}
for l in open(a.gt):
    if l.strip() and not l.startswith('#'): w = l.split(); gt[int(w[0])] = int(w[1])
L = np.loadtxt(a.latent, comments='#'); pind = L[:, 0].astype(int); z1, z2, am, gm = L[:, 1], L[:, 2], L[:, 3].astype(int), L[:, 5]
conf = np.array([gt.get(p, 0) for p in pind]); on = gm > 0.5; nc = a.n1 * a.n2
cell = np.minimum(am - 1, nc - 1) if a.argmax else (np.clip(np.rint(z1), 0, a.n1 - 1).astype(int) * a.n2 + np.clip(np.rint(z2), 0, a.n2 - 1).astype(int))
m = on & (conf > 0) & ((am - 1) < nc if a.argmax else True); cell, conf = cell[m], conf[m]; cs = np.unique(conf); C = np.zeros((nc, len(cs)))
for j, c in enumerate(cs): C[:, j] = np.bincount(cell[conf == c], minlength=nc)
tot = C.sum(1); P = C / np.maximum(tot[:, None], 1); maj = P.argmax(1); base = C.sum(0) / C.sum()
print('%s%d conformers, %d on-grid particles, %d used cells' % (a.tag + ': ' if a.tag else '', len(cs), int(C.sum()), int((tot > 0).sum())))
print('  conf  base  best-cell purity (x base)  cells owned  recall in owned cells  share of its particles in its owned cells')
owned_any = 0
for j, c in enumerate(cs):
    own = np.where((maj == j) & (tot > 0))[0]; rec = C[own, j].sum() / C[:, j].sum() if len(own) else 0.; owned_any += len(own) > 0
    print('  %3d  %.2f   %.2f (%.1fx)                 %2d          %.2f' % (c, base[j], P[:, j].max(), P[:, j].max() / base[j], len(own), rec))
# one-to-one: Hungarian on the purity matrix (each conformer its own best distinct cell)
r_, c_ = linear_sum_assignment(-P[tot > 0]); cells_ok = np.where(tot > 0)[0]; pur = P[tot > 0][r_, c_]
print('  conformers owning at least one cell (majority): %d of %d; one-cell-per-conformer assignment: purity per conformer %s (mean %.2f, min %.2f)' % (owned_any, len(cs), ' '.join('%d:%.2f' % (cs[k], v) for k, v in zip(c_, pur)), pur.mean(), pur.min()))
