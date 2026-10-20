#!/usr/bin/env python3
"""Bins along the main axis of a klem latent -> selection file for a SIMPLE project (one int per project row; 0 = not in
a bin). The main axis = first principal direction of the on-grid posterior means. usage:
klem_bins.py KLEMDIR META OUT_STATES --nb 7 [--minmass 0.5]"""
import sys, argparse, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('klemdir'); ap.add_argument('meta'); ap.add_argument('out'); ap.add_argument('--nb', type=int, default=7); ap.add_argument('--minmass', type=float, default=0.5); ap.add_argument('--pc', type=int, default=1, help='principal direction of the posterior means to bin along (1 = main axis, 2 = second axis)'); ap.add_argument('--nptcl', type=int, default=0, help='project rows (toy: instead of META)')
a = ap.parse_args(); N = a.nptcl if a.nptcl else len(np.load(a.meta)['flags'])
L = np.loadtxt(a.klemdir + '/klem_latent.txt', comments='#'); pind = L[:, 0].astype(int) - 1; z = (L[:, [1, 2, 7]] if L.shape[1] >= 8 else L[:, 1:3]); gm = L[:, 5]; ok = gm >= a.minmass
zc = z[ok] - z[ok].mean(0); u = np.linalg.svd(zc, full_matrices=False)[2][a.pc - 1]; t = zc @ u
edges = np.quantile(t, np.linspace(0, 1, a.nb + 1)); b = np.clip(np.searchsorted(edges, t, side='right') - 1, 0, a.nb - 1) + 1
states = np.zeros(N, int); states[pind[ok]] = b; np.savetxt(a.out, states, fmt='%d')
print('axis PC%d direction %s;' % (a.pc, np.array2string(u, precision=2)), end=' '); print(' %d particles binned into %d bins of %s; %d left out (grid mass < %.1f)' % (ok.sum(), a.nb, np.bincount(b)[1:], (~ok).sum(), a.minmass))
