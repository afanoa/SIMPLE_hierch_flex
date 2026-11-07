#!/usr/bin/env python3
"""Truth readout of the factorised (one focused board per data-derived region) route on a hinged-arm set. Per region:
its overlap with every arm's sector mask; its board's |corr(z1/z2, angle)| for every arm angle and the 2-fold
cross-validated R^2 of each arm's (theta, phi) from a quadratic in (z1, z2) (the continuous coordinate as a predictor);
then per arm the best region (coverage: which arms the route reaches and how well). Truth is used only here.
usage: factorised_readout.py RUN GT_PARTICLES TRUTH MASKDIR   (RUN/regions/region_*.mrc, RUN/arm_<region>/board)"""
import sys, os, glob, re, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from model_class import read_mrc
run, gtf, trf, mk = sys.argv[1:5]
G = np.loadtxt(gtf, comments='#'); conf = np.zeros(int(G[:, 0].max()) + 1, int); conf[G[:, 0].astype(int)] = G[:, 1].astype(int)
ang = {int(w[0]): np.array([float(x) for x in w[1:]]) for w in (l.split() for l in open(trf)) if w and not w[0].startswith('#')}
narm = len(next(iter(ang.values()))) // 2; S = [read_mrc('%s/arm_%d.mrc' % (mk, k + 1)) for k in range(narm)]
def r2cv(z, y):
    X = np.c_[np.ones(len(z)), z, z ** 2, z[:, :1] * z[:, 1:2]]; out = 0.; sst = 0.
    for f in (0, 1):
        tr, te = (np.arange(len(z)) % 2) != f, (np.arange(len(z)) % 2) == f
        b = np.linalg.lstsq(X[tr], y[tr], rcond=None)[0]; out += ((y[te] - X[te] @ b) ** 2).sum(); sst += ((y[te] - y[tr].mean()) ** 2).sum()
    return 1 - out / sst
best = {}
for rf in sorted(glob.glob(run + '/regions/region_*.mrc'), key=lambda p: int(re.search(r'region_(\d+)', p).group(1))):
    rn = int(re.search(r'region_(\d+)', rf).group(1)); R = read_mrc(rf); ov = [(R * s).sum() / max(R.sum(), 1e-9) for s in S]
    bd = '%s/arm_%d/board' % (run, rn)
    if not os.path.exists(bd + '/klem_latent.txt'): print('region %d: overlap with arm sectors %s; no board' % (rn, ' '.join('%.2f' % v for v in ov))); continue
    L = np.loadtxt(bd + '/klem_latent.txt', comments='#'); p = L[:, 0].astype(int); z = L[:, 1:3]; c = conf[p]; ok = c > 0
    A = np.array([ang[k] for k in c[ok]]); zz = z[ok]
    cc = [max(abs(np.corrcoef(zz[:, j], A[:, k])[0, 1]) for j in range(2)) for k in range(2 * narm)]
    r2 = [r2cv(zz, A[:, k]) for k in range(2 * narm)]; k = int(np.argmax([r2[2 * a] + r2[2 * a + 1] for a in range(narm)]))
    print('region %d (mask overlap with arm sectors %s): board on %d particles (%d junk left)' % (rn, ' '.join('%.2f' % v for v in ov), len(p), (~ok).sum()))
    print('   |corr(z, angle)| %s' % ' '.join('a%d %.2f/%.2f' % (a + 1, cc[2 * a], cc[2 * a + 1]) for a in range(narm)))
    print('   CV R^2 (theta/phi)  %s  -> tracks arm %d' % (' '.join('a%d %.2f/%.2f' % (a + 1, r2[2 * a], r2[2 * a + 1]) for a in range(narm)), k + 1))
    for a in range(narm):
        if r2[2 * a] + r2[2 * a + 1] > best.get(a, (-9, 0, 0, 0))[0]: best[a] = (r2[2 * a] + r2[2 * a + 1], rn, r2[2 * a], r2[2 * a + 1])
print('== per arm, best region: %s' % '; '.join('arm %d: region %d R^2 theta %.2f phi %.2f' % (a + 1, best[a][1], best[a][2], best[a][3]) if a in best else 'arm %d: none' % (a + 1) for a in range(narm)))
