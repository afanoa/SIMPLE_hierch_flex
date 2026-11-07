#!/usr/bin/env python3
"""land_latsep.py: does a flex latent carry the truth at all? For two truth groups (simsites: arm closed vs open via combos.txt column 4;
swing30: t <= T0 vs t >= T1), the Fisher discriminant separation d' = sqrt(dmu^T W^-1 dmu) (W = pooled within-group covariance) of the
posterior-mean latents (flex_pca_coordinates.txt) and of the raw MAP latents (flex_pca_embedding.bin), the per-dimension |dmu| / sd, and
the held-out (eo) classification accuracy of the linear discriminant.  usage: land_latsep.py --emb EMB --gt GT (--combos C | --t 5:25)"""
import os, argparse, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('--emb', required=True); ap.add_argument('--gt', required=True); ap.add_argument('--combos', default='')
ap.add_argument('--t', default='5:25'); ap.add_argument('--codes', default='200:230'); a = ap.parse_args()
b = open(a.emb, 'rb').read(); n, D = [int(x) for x in np.frombuffer(b, np.int32, 2, 12)]; pind = np.frombuffer(b, np.int32, n, 28); zr = np.frombuffer(b, np.float64, n * D, 28 + 4 * n).reshape(D, n).T
C = np.loadtxt(os.path.join(os.path.dirname(a.emb), 'flex_pca_coordinates.txt'), comments='#'); eo = C[:, 1].astype(int); zp = C[:, 5:]
conf = np.loadtxt(a.gt, comments='#')[:, 1].astype(int)[pind - 1]
if a.combos:
    cmb = {int(l.split()[0]): int(l.split()[3]) for l in open(a.combos) if l.strip() and not l.startswith('#')}; g = np.array([cmb.get(c, -1) if c > 0 else -1 for c in conf])
    lab = 'arm closed (1) vs open (0), non-junk'
else:
    i0, i1 = [int(x) for x in a.codes.split(':')]; t0, t1 = [int(x) for x in a.t.split(':')]; tt = conf - i0; g = np.where((tt >= 0) & (tt <= t0), 0, np.where((tt >= t1) & (tt <= i1 - i0), 1, -1))
    lab = 'swing t <= %d (0) vs t >= %d (1)' % (t0, t1)
print('%s: %s; group sizes %d / %d' % (a.emb, lab, (g == 0).sum(), (g == 1).sum()))
for nm, Z in (('posterior-mean', zp), ('raw MAP', zr)):
    ok = g >= 0; Zg, gg, ee = Z[ok], g[ok], eo[ok]; m0, m1 = Zg[gg == 0].mean(0), Zg[gg == 1].mean(0); W = 0.5 * (np.cov(Zg[gg == 0].T) + np.cov(Zg[gg == 1].T)); dm = m1 - m0
    dp = float(np.sqrt(dm @ np.linalg.solve(W, dm))); per = np.abs(dm) / np.sqrt(np.diag(W))
    acc = []
    for h in (0, 1):   # train on one half, test on the other (balanced accuracy at the prior-free midpoint)
        tr, te = ee == h, ee != h; a0, a1 = Zg[tr & (gg == 0)].mean(0), Zg[tr & (gg == 1)].mean(0); Wt = 0.5 * (np.cov(Zg[tr & (gg == 0)].T) + np.cov(Zg[tr & (gg == 1)].T))
        w = np.linalg.solve(Wt, a1 - a0); s = Zg[te] @ w; thr = 0.5 * (a0 + a1) @ w; pr = (s > thr).astype(int)
        acc.append(0.5 * (((pr == 0) & (gg[te] == 0)).sum() / max(1, (gg[te] == 0).sum()) + ((pr == 1) & (gg[te] == 1)).sum() / max(1, (gg[te] == 1).sum())))
    print('  %-14s: d\' %.2f; |dmu|/sd per dim %s; held-out balanced accuracy %.3f / %.3f' % (nm, dp, ' '.join('%.2f' % x for x in per), acc[0], acc[1]))
