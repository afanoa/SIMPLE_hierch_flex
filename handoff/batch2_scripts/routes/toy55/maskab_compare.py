#!/usr/bin/env python3
"""Run-to-run spread vs mask effect for the focused refinement at one node: three labellings (original current-mask run, current-mask
replicate, children-mask run) -> per run the conformer NMI / purity and swing eta2 (truth), and pairwise agreement (NMI, share of
particles in the same state; the runs share their starting refs so state numbers correspond).
usage: maskab_compare.py GT name=labels name=labels name=labels"""
import sys, numpy as np
gt = sys.argv[1]; runs = [a.split('=', 1) for a in sys.argv[2:]]
L = {nm: np.loadtxt(fn).astype(int) for nm, fn in runs}; n = len(next(iter(L.values())))
G = np.loadtxt(gt, comments='#', ndmin=2); conf = np.zeros(n, int); conf[G[:, 0].astype(int) - 1] = G[:, 1].astype(int)
def nmi(x, y):
    _, xi = np.unique(x, return_inverse=True); _, yi = np.unique(y, return_inverse=True); J = np.zeros((xi.max() + 1, yi.max() + 1)); np.add.at(J, (xi, yi), 1); J /= J.sum()
    px, py = J.sum(1), J.sum(0); nz = J > 0; I = (J[nz] * np.log(J[nz] / np.outer(px, py)[nz])).sum(); h = lambda p: -(p[p > 0] * np.log(p[p > 0])).sum(); return I / max(np.sqrt(h(px) * h(py)), 1e-12)
def eta2(v, lab):
    u, inv = np.unique(lab, return_inverse=True); m = np.bincount(inv, v) / np.bincount(inv); return 1 - ((v - m[inv]) ** 2).sum() / ((v - v.mean()) ** 2).sum()
for nm, lab in L.items():
    m = lab > 0; c = conf[m]; l = lab[m]; pur = sum(np.bincount(c[l == s]).max() for s in np.unique(l)) / m.sum(); sw = m & (conf >= 200) & (conf <= 230)
    print('%-28s conformer NMI %.3f purity %.3f | swing eta2 %.3f (%d swing particles)' % (nm, nmi(c, l), pur, eta2(conf[sw] - 200., lab[sw]) if sw.sum() > 100 else np.nan, sw.sum()))
names = list(L)
for i in range(len(names)):
    for j in range(i + 1, len(names)):
        a, b = L[names[i]], L[names[j]]; m = (a > 0) & (b > 0)
        print('agreement %-26s vs %-26s NMI %.3f, same state %.1f %%' % (names[i], names[j], nmi(a[m], b[m]), 100 * (a[m] == b[m]).mean()))
