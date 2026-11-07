#!/usr/bin/env python3
"""Deterministic-annealing split criterion for the genetic expansion (Rose 1998): a family of cells splits when it is
unstable, along the direction where it is. The within-family covariance of the cell deviations is estimated CROSS-HALF
(Gram matrix 1/2 (X0 X1^H + X1 X0^H) of the population-weighted, centred even/odd deviations, so noise cannot create a
mode): its leading eigenvalue is the reproducible within-family variance along the principal mode and its eigenvector
the split direction (cells by the sign of their coordinate). The noise along that mode comes from the half DIFFERENCES
(1/4 (X0 - X1)(X0 - X1)^H along the same vector). The critical-temperature rule: split while the reproducible variance
along the principal mode exceeds tau x the noise variance along it. Recursion from one family of all grid cells; the
number of families is decided by the data. Writes family<j>_subset.txt (1-based particle indices) + families.txt.
usage: anneal_split.py KLEM_DIR OUT_DIR [--tau 1] [--min_cells 3] [--min_frac 0.05] [--max_fam 8] [--gt gt_particles.txt]"""
import sys, os, argparse, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('klem'); ap.add_argument('out'); ap.add_argument('--tau', type=float, default=1.); ap.add_argument('--min_cells', type=int, default=3)
ap.add_argument('--min_frac', type=float, default=0.05); ap.add_argument('--max_fam', type=int, default=8); ap.add_argument('--gt', default=''); a = ap.parse_args()
H = np.load(a.klem + '/klem_halves.npz'); D0 = H['D0'].astype(complex); D1 = H['D1'].astype(complex); w = H['w'].astype(float); i1 = H['i1']; i2 = H['i2']; ncell = len(w)
L = np.loadtxt(a.klem + '/klem_latent.txt', comments='#'); pind = L[:, 0].astype(int); am = L[:, 3].astype(int) - 1; gm = L[:, 5]; on = (gm > 0.5) & (am < ncell)
gt = {}
if a.gt:
    for l in open(a.gt):
        if l.strip() and not l.startswith('#'): t = l.split(); gt[int(t[0])] = int(t[1])


def node_stats(cells):
    """leading eigenpair of the population-weighted cross-half within-family covariance (Gram form) and the noise variance along it"""
    ww = w[cells] / max(w[cells].sum(), 1e-30); X0 = D0[cells] - (ww[:, None] * D0[cells]).sum(0); X1 = D1[cells] - (ww[:, None] * D1[cells]).sum(0)
    sw = np.sqrt(ww)[:, None]; X0 = X0 * sw; X1 = X1 * sw
    C = 0.5 * np.real(X0 @ np.conj(X1).T + X1 @ np.conj(X0).T); Nn = 0.25 * np.real((X0 - X1) @ np.conj(X0 - X1).T)
    ev, U = np.linalg.eigh(C); u = U[:, -1]; lam = float(ev[-1]); noise = float(u @ Nn @ u); coord = u / np.maximum(np.sqrt(ww), 1e-12)
    return lam, noise, coord


fams = []; queue = [np.arange(ncell)]
print('annealing split of %s: %d cells, tau %.2f, min %d cells and %.0f%% of the parent mass per child' % (a.klem, ncell, a.tau, a.min_cells, 100 * a.min_frac))
while queue:
    cells = queue.pop(0); A = B = None
    if len(cells) >= 2 * a.min_cells: lam, noise, coord = node_stats(cells)
    else: lam, noise, coord = 0., 1., None
    ratio = lam / max(noise, 1e-30); ok = ratio > a.tau and coord is not None and len(fams) + len(queue) + 1 < a.max_fam
    if ok:
        A = cells[coord > 0]; B = cells[coord <= 0]
        for _ in range(10):   # settle (Rose: the perturbed pair re-assigns): every cell to the nearer sub-family mean, distance cross-half so its own noise cannot pull it
            if len(A) == 0 or len(B) == 0: break
            mu = {}
            for nm, sub in (('A', A), ('B', B)): ws = w[sub] / max(w[sub].sum(), 1e-30); mu[nm] = ((ws[:, None] * D0[sub]).sum(0), (ws[:, None] * D1[sub]).sum(0))
            dA = np.real(np.sum(np.conj(D0[cells] - mu['A'][0]) * (D1[cells] - mu['A'][1]), 1)); dB = np.real(np.sum(np.conj(D0[cells] - mu['B'][0]) * (D1[cells] - mu['B'][1]), 1))
            A2 = cells[dA <= dB]; B2 = cells[dA > dB]
            if len(A2) == len(A) and np.array_equal(np.sort(A2), np.sort(A)): break
            A, B = A2, B2
        if len(A) < a.min_cells or len(B) < a.min_cells or w[A].sum() < a.min_frac * w[cells].sum() or w[B].sum() < a.min_frac * w[cells].sum(): ok = False
    print('  node of %2d cells (mass %.2f): reproducible variance along the principal mode %.3g, noise along it %.3g, ratio %.2f -> %s' % (len(cells), w[cells].sum() / w.sum(), lam, noise, ratio, ('SPLIT %d | %d cells' % (len(A), len(B))) if ok else 'leaf'))
    if ok: queue += [A, B]
    else: fams.append(cells)
os.makedirs(a.out, exist_ok=True); fam_of = np.full(ncell, -1)
for j, cells in enumerate(fams): fam_of[cells] = j
famp = np.where(on, fam_of[np.minimum(am, ncell - 1)], -1)
tot = {}
if gt:
    for p in pind[on]: tot[gt.get(int(p), 0)] = tot.get(gt.get(int(p), 0), 0) + 1
for j, cells in enumerate(fams):
    members = pind[famp == j]; np.savetxt('%s/family%d_subset.txt' % (a.out, j + 1), members, fmt='%d'); comp = ''
    if gt:
        c = np.array([gt.get(int(p), 0) for p in members]); cs, cn = np.unique(c[c > 0], return_counts=True)
        comp = ' | share of each conformer here: ' + ' '.join('%d:%.2f' % (k, nn / max(tot.get(k, 1), 1)) for k, nn in zip(cs, cn) if nn / max(tot.get(k, 1), 1) >= 0.05)
    print('family %d: %d cells [%s], %d particles%s' % (j + 1, len(cells), ' '.join('(%g,%g)' % (i1[k], i2[k]) for k in cells), len(members), comp))
np.savetxt(a.out + '/families.txt', np.c_[np.arange(1, ncell + 1), fam_of + 1], fmt='%d', header='state family'); print('%d families -> %s' % (len(fams), a.out))
