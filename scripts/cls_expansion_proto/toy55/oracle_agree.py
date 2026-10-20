#!/usr/bin/env python3
"""Oracle-start agreement of a klem board in the BOARD's own lattice numbering (argmax m -> i1 = (m-1)//n2, i2 = (m-1)%n2;
the oracle's init_labels use g1 + n1*g2 + 1, a different order). Non-junk particles: exact cell, within one lattice step,
free-component share, and the correlation of the posterior-mean coordinates z1/z2 with the start g1/g2.
usage: oracle_agree.py ORACLE_DIR [ORACLE_DIR ...]"""
import sys, os, numpy as np
for d in sys.argv[1:]:
    lat = d + '/board/klem_latent.txt'
    if not os.path.isfile(lat): continue
    cfg = dict(l.split(None, 1) for l in open(d + '/board/board.txt') if l.strip() and len(l.split()) > 1)
    n1, n2 = int(cfg['n1']), int(cfg['n2']); ng = n1 * n2
    P = np.loadtxt(d + '/init_positions.txt', dtype=int); L = np.loadtxt(d + '/init_labels.txt', dtype=int)
    g1 = np.full(len(L), -1); g2 = np.full(len(L), -1); g1[P[:, 0] - 1] = P[:, 1]; g2[P[:, 0] - 1] = P[:, 2]
    B = np.loadtxt(lat, comments='#'); pind = B[:, 0].astype(int) - 1; z1, z2, am = B[:, 1], B[:, 2], B[:, 3].astype(int)
    ok = L[pind] > 0; a1, a2 = (am - 1) // n2, (am - 1) % n2; grid = am <= ng
    s1, s2 = g1[pind], g2[pind]; m = ok & grid
    ex = np.mean((a1[m] == s1[m]) & (a2[m] == s2[m])); nb = np.mean(np.maximum(np.abs(a1[m] - s1[m]), np.abs(a2[m] - s2[m])) <= 1)
    c1 = np.corrcoef(z1[ok], s1[ok])[0, 1]; c2 = np.corrcoef(z2[ok], s2[ok])[0, 1]
    print('%-34s exact cell %.3f (chance %.3f), within 1 step %.3f, free %.3f; corr(z1, g1) %.3f corr(z2, g2) %.3f' % (os.path.basename(d.rstrip('/')), ex, 1. / ng, nb, np.mean(~grid[ok]), c1, c2))
