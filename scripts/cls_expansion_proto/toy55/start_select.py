#!/usr/bin/env python3
"""Multi-start selection for region boards: in every RUN/arm_<r>/ the start board_s<k> with the highest mean max
responsibility (klem_latent column maxr; the board's own confidence, no truth) becomes RUN/arm_<r>/board (symlink).
--gt GT --truth TRUTH: also prints each start's best-arm cross-validated R^2 (theta + phi) so the rule can be judged.
usage: start_select.py RUN [--gt GT_PARTICLES --truth TRUTH]"""
import sys, os, glob, re, argparse, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('run'); ap.add_argument('--gt', default=''); ap.add_argument('--truth', default=''); a = ap.parse_args()
if a.gt:
    G = np.loadtxt(a.gt, comments='#'); conf = np.zeros(int(G[:, 0].max()) + 1, int); conf[G[:, 0].astype(int)] = G[:, 1].astype(int)
    ang = {int(w[0]): np.array([float(x) for x in w[1:]]) for w in (l.split() for l in open(a.truth)) if w and not w[0].startswith('#')}; narm = len(next(iter(ang.values()))) // 2
def r2cv(z, y):
    X = np.c_[np.ones(len(z)), z, z ** 2, z[:, :1] * z[:, 1:2]]; e = s = 0.
    for f in (0, 1):
        tr, te = (np.arange(len(z)) % 2) != f, (np.arange(len(z)) % 2) == f; b = np.linalg.lstsq(X[tr], y[tr], rcond=None)[0]
        e += ((y[te] - X[te] @ b) ** 2).sum(); s += ((y[te] - y[tr].mean()) ** 2).sum()
    return 1 - e / s
for d in sorted(glob.glob(a.run + '/arm_*'), key=lambda p: int(re.search(r'arm_(\d+)$', p).group(1))):
    sc = []
    for b in sorted(glob.glob(d + '/board_s*/klem_latent.txt')):
        L = np.loadtxt(b, comments='#'); conf_mean = L[:, 4].mean(); line = '%s mean maxr %.4f' % (b.split('/')[-2], conf_mean)
        if a.gt:
            c = conf[L[:, 0].astype(int)]; ok = c > 0; A = np.array([ang[k] for k in c[ok]]); z = L[ok, 1:3]
            r = [(r2cv(z, A[:, 2 * k]), r2cv(z, A[:, 2 * k + 1])) for k in range(narm)]; k = int(np.argmax([u + v for u, v in r])); line += '  (truth: arm %d R^2 %.2f/%.2f)' % (k + 1, *r[k])
        sc.append((conf_mean, os.path.dirname(b), line))
    if not sc: print('%s: no start finished' % os.path.basename(d)); continue
    best = max(sc)[1]; lk = d + '/board'
    if os.path.islink(lk) or os.path.exists(lk): os.remove(lk)
    os.symlink(os.path.basename(best), lk); print('%s: %s' % (os.path.basename(d), '; '.join(('*' if s[1] == best else ' ') + s[2] for s in sc)))
