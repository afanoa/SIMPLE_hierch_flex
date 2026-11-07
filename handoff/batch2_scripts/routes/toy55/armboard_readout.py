#!/usr/bin/env python3
"""Readout of one arm-focused board on a hinged-arm synthetic set: |corr| of the lattice coordinates (z1, z2) with every
arm angle, and the per-cell (argmax) std ratio of the FOCUSED arm's theta / phi (cells with >= MINC particles, conformer
particles only), median and best, plus the junk share captured by the free components.
usage: armboard_readout.py BOARD_DIR GT_PARTICLES TRUTH ARM [MINC 200] [NCELL 21 = lattice cells; higher = free components]"""
import sys, numpy as np
bd, gtf, trf, arm = sys.argv[1:5]; arm = int(arm); minc = int(sys.argv[5]) if len(sys.argv) > 5 else 200; ncell = int(sys.argv[6]) if len(sys.argv) > 6 else 21
G = np.loadtxt(gtf, comments='#'); conf = np.zeros(int(G[:, 0].max()) + 1, int); conf[G[:, 0].astype(int)] = G[:, 1].astype(int)
ang = {int(w[0]): np.array([float(x) for x in w[1:]]) for w in (l.split() for l in open(trf)) if w and not w[0].startswith('#')}
narm = len(next(iter(ang.values()))) // 2
L = np.loadtxt(bd + '/klem_latent.txt', comments='#'); p = L[:, 0].astype(int); z = L[:, 1:3]; am = L[:, 3].astype(int); c = conf[p]
ok = c > 0; A = np.array([ang[k] for k in c[ok]]); sall = A.std(0)
cc = [max(abs(np.corrcoef(z[ok, j], A[:, k])[0, 1]) for j in range(2)) for k in range(2 * narm)]
own = [2 * (arm - 1), 2 * (arm - 1) + 1]; R = []
for m in np.unique(am[ok & (am <= ncell)]):
    s = am[ok] == m
    if s.sum() >= minc: R.append((A[s][:, own].std(0) / sall[own]))
R = np.array(R) if R else np.ones((1, 2)); junk = (~ok).sum(); fj = 1 - ok[am > ncell].mean() if (am > ncell).any() else float('nan')
print('arm %d board %s: |corr(z, angle)| own arm theta %.2f phi %.2f; best other %.2f  (all %s)' % (arm, bd.rstrip('/').split('/')[-2], cc[own[0]], cc[own[1]], max(v for k, v in enumerate(cc) if k not in own), ' '.join('%.2f' % v for v in cc)))
print('   cells >= %d ptcls: %d; own-arm std ratio theta median %.2f best %.2f, phi median %.2f best %.2f; junk %d particles, free-component purity %.2f' % (minc, len(R), np.median(R[:, 0]), R[:, 0].min(), np.median(R[:, 1]), R[:, 1].min(), junk, fj))
