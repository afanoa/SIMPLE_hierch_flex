#!/usr/bin/env python3
"""Seeds for the minority test on a synthetic set with truth: per particle a noisy closed score s = d' x closed + N(0, 1) (AUC of s for
closed vs the rest = Phi(d'/sqrt 2); d' 1.47 -> 0.85, the per-particle separability measured on real RNase P), cut into NBINS equal-width
bins between the 0.5 and 99.5 % quantiles (the outer bins take the tails), the same edges in both halves; writes one label (1..NBINS) per
half row and prints each bin's size and closed share.
usage: noisy_seeds.py GT_A GT_B OUT_A OUT_B [--dprime 1.47 --nbins 12 --seed 0 --closed_col 2]"""
import argparse, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('gta'); ap.add_argument('gtb'); ap.add_argument('outa'); ap.add_argument('outb')
ap.add_argument('--dprime', type=float, default=1.47); ap.add_argument('--nbins', type=int, default=12); ap.add_argument('--seed', type=int, default=0); ap.add_argument('--closed_col', type=int, default=2)
a = ap.parse_args(); rng = np.random.default_rng(a.seed); S, CL = [], []
for f in (a.gta, a.gtb):
    G = np.loadtxt(f, comments='#'); n = int(G[:, 0].max()); cl = np.zeros(n, bool); cl[G[:, 0].astype(int) - 1] = G[:, a.closed_col] > 0
    S.append(a.dprime * cl + rng.standard_normal(n)); CL.append(cl)
lo, hi = np.quantile(np.r_[S[0], S[1]], [0.005, 0.995]); edges = np.linspace(lo, hi, a.nbins + 1)
for s, cl, out, t in zip(S, CL, (a.outa, a.outb), 'AB'):
    lab = np.clip(np.searchsorted(edges, s, side='right'), 1, a.nbins); np.savetxt(out, lab, fmt='%d')
    print('half %s: %d particles, %d closed; bins (size, closed share): %s' % (t, len(s), cl.sum(), '  '.join('%d %.2f' % ((lab == b).sum(), cl[lab == b].mean() if (lab == b).any() else 0.) for b in range(1, a.nbins + 1))))
from math import erf, sqrt
print("d' %.2f: AUC %.3f" % (a.dprime, 0.5 * (1 + erf(a.dprime / 2.))))
