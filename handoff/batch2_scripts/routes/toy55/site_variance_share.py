#!/usr/bin/env python3
"""usage: site_variance_share.py MASKDIR VARMAP [VARMAP ..]: share of each variance map's positive mass inside the true site masks (MASKDIR/
mask_*.mrc > 0.5) and each site's peak / the global peak; plus the noise-free frequency-weighted variance of MASKDIR's combination maps
(combos.txt: id flags.. frequency; conf<id>_sites.mrc) as the infinite-data ceiling."""
import sys, glob, os, numpy as np
def rd(p): h = np.fromfile(p, np.int32, 256); n = int(h[0]); return np.fromfile(p, np.float32, offset=1024 + int(h[23]), count=n ** 3).reshape(n, n, n).astype(np.float64)
MD = sys.argv[1]; MS = {os.path.basename(p)[5:-4]: rd(p) > 0.5 for p in sorted(glob.glob(MD + '/mask_*.mrc'))}
def report(tag, v):
    vp = np.clip(v, 0, None); tot = vp.sum(); g = v.max()
    print('%s: ' % tag + ', '.join('%s share %.3f peak/max %.2f' % (k, vp[m].sum() / tot, v[m].max() / g) for k, m in MS.items()) + ' (core cut 0.20 of max)')
for p in sys.argv[2:]: report(p.split('1_hetero_tree/')[-1], rd(p))
C = [l.split() for l in open(MD + '/combos.txt') if not l.startswith('#')]; w = np.array([float(c[-1]) for c in C]); V = [rd('%s/conf%02d_sites.mrc' % (MD, int(c[0]))) for c in C]
mean = sum(wi * x for wi, x in zip(w, V)) / w.sum(); report('noise-free truth (8 combination maps, frequency-weighted)', sum(wi * (x - mean) ** 2 for wi, x in zip(w, V)) / w.sum())
