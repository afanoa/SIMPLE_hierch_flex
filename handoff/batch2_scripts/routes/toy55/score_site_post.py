#!/usr/bin/env python3
"""Per-site posterior vs every truth variable (continuous or discrete): for each site dir (site_size.py output: post_master.txt, top_star.txt)
the |corr| of P(class 2) with every true site variable over the non-junk rows in the site EM, the best-matching variable, the eta2 of the
binary class over it, and the oracle ceiling when given (ORACLE: 'name value' lines). usage: score_site_post.py TRUTH SITEDIR [SITEDIR ...]"""
import sys, numpy as np
tr = sys.argv[1]; names = [l.split()[1:] for l in open(tr) if l.startswith('#')][0][3:]; T = np.genfromtxt(tr, comments='#'); X = T[:, 3:]; ok0 = T[:, 2] == 0
for d in sys.argv[2:]:
    p = np.loadtxt(d + '/post_master.txt'); top = np.loadtxt(d + '/top_star.txt', dtype=int) > 0; ok = ok0 & (p >= 0)
    cs = [abs(np.corrcoef(p[ok], X[ok, j])[0, 1]) if np.nanstd(X[ok, j]) > 0 else 0. for j in range(X.shape[1])]; j = int(np.nanargmax(cs))
    v = X[ok, j]; t = top[ok]; e2 = ((t.sum() * (v[t].mean() - v.mean()) ** 2 + (~t).sum() * (v[~t].mean() - v.mean()) ** 2) / (len(v) * v.var()))
    print('%s: best variable %s |corr(P, x)| %.3f (next %s), class-2 rows %d, eta2 of the class over it %.3f' % (d, names[j], cs[j], ', '.join('%s %.2f' % (names[k], cs[k]) for k in np.argsort(cs)[::-1][1:4]), t.sum(), e2))
    import os
    if os.path.exists(d + '/motion_coord.txt'):
        mc = np.loadtxt(d + '/motion_coord.txt'); okm = ok0 & (mc >= 0); cm = [abs(np.corrcoef(mc[okm], X[okm, k])[0, 1]) if np.nanstd(X[okm, k]) > 0 else 0. for k in range(X.shape[1])]; k = int(np.nanargmax(cm))
        print('   motion coordinate: best variable %s |corr| %.3f (next %s)' % (names[k], cm[k], ', '.join('%s %.2f' % (names[q], cm[q]) for q in np.argsort(cm)[::-1][1:4])))
