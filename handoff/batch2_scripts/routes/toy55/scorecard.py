#!/usr/bin/env python3
"""Truth scorecard of a labelling (per row, 0 = none) on a synthetic benchmark: number of states, NMI(conformer; state), purity by
conformer and by class (CLASSES: name:lo-hi ranges of conformer ids), and per truth coordinate (CODES as in factor_sample.py:
name:id0:step for one motion per conformer range, grid:id0:kstep:jstep) the share of its variance the states explain (eta2) and the
extent the states reach (lowest / highest state mean among states holding >= 50 % of that motion, >= MINN particles).
usage: scorecard.py LABELS GT [--classes junk:0-0,apo:1-100,...] [--codes swing:200:1] [--minn 200] [--tag NAME]"""
import sys, argparse, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('labels'); ap.add_argument('gt'); ap.add_argument('--classes', default=''); ap.add_argument('--codes', default='')
ap.add_argument('--minn', type=int, default=200); ap.add_argument('--tag', default=''); a = ap.parse_args()
lab = np.loadtxt(a.labels, dtype=int); G = np.loadtxt(a.gt, comments='#'); conf = np.full(len(lab), -1); conf[G[:, 0].astype(int)[G[:, 0] <= len(lab)] - 1] = G[:, 1].astype(int)[G[:, 0] <= len(lab)]
on = (lab > 0) & (conf >= 0); L, C = lab[on], conf[on]; S = np.unique(L)
if len(L) == 0: print('%sno labelled particles' % ((a.tag + ': ') if a.tag else '')); sys.exit(0)
def nmi(x, y):
    _, xi = np.unique(x, return_inverse=True); _, yi = np.unique(y, return_inverse=True); J = np.zeros((xi.max() + 1, yi.max() + 1)); np.add.at(J, (xi, yi), 1); J /= J.sum()
    px, py = J.sum(1), J.sum(0); nz = J > 0; I = (J[nz] * np.log(J[nz] / np.outer(px, py)[nz])).sum(); hx, hy = -(px * np.log(px)).sum(), -(py * np.log(py)).sum(); return I / max(np.sqrt(hx * hy), 1e-12)
def purity(cls): return sum(np.bincount(cls[L == s]).max() for s in S) / len(L)
_, Ci = np.unique(C, return_inverse=True); out = ['%s%d states over %d particles; NMI(conformer; state) %.2f, conformer purity %.2f' % ((a.tag + ': ') if a.tag else '', len(S), len(L), nmi(C, L), purity(Ci))]
if a.classes:
    rng = [(w.split(':')[0], int(w.split(':')[1].split('-')[0]), int(w.split(':')[1].split('-')[1])) for w in a.classes.split(',') if w]
    cl = np.full(len(C), -1)
    for i, (nm, lo, hi) in enumerate(rng): cl[(C >= lo) & (C <= hi)] = i
    out.append('class purity %.2f (%s)' % (purity(np.where(cl >= 0, cl, len(rng))), ', '.join('%s %d' % (nm, (cl == i).sum()) for i, (nm, _, _) in enumerate(rng))))
tv = {}
for c in [w for w in a.codes.split(',') if w]:
    p = c.split(':')
    if p[0] == 'grid':
        i0, ks, js = int(p[1]), float(p[2]), float(p[3]); m = C >= i0; tv['grid_k'] = np.where(m, ((C - i0) // 10) * ks, np.nan); tv['grid_j'] = np.where(m, ((C - i0) % 10) * js, np.nan)
    else:
        nm, i0, st = p[0], int(p[1]), float(p[2]); hi = min([int(q.split(':')[1]) for q in a.codes.split(',') if q and not q.startswith('grid:') and int(q.split(':')[1]) > i0] + [10 ** 9])
        tv[nm] = np.where((C >= i0) & (C < hi), (C - i0) * st, np.nan)
for nm, x in tv.items():
    m = np.isfinite(x)
    if m.sum() < 50: continue
    xs, ls = x[m], L[m]; mu = xs.mean(); bet = sum((ls == s).sum() * (xs[ls == s].mean() - mu) ** 2 for s in np.unique(ls)); eta = bet / max(((xs - mu) ** 2).sum(), 1e-12)
    ms = [(xs[ls == s].mean(), (ls == s).sum()) for s in np.unique(ls) if (ls == s).sum() >= a.minn and (ls == s).sum() >= 0.5 * (L == s).sum()]
    ext = ('%.1f .. %.1f of %.1f .. %.1f in %d states' % (min(v for v, _ in ms), max(v for v, _ in ms), np.nanmin(xs), np.nanmax(xs), len(ms))) if ms else 'no state holds it'
    out.append('%s: eta2 %.2f, extent %s' % (nm, eta, ext))
print(' | '.join(out))
