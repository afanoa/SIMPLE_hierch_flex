#!/usr/bin/env python3
"""usage: motion_bins2.py POST OUT_STATES VARIANT [TRUTH ANGLE]   (motion_bins.py reading the pole-1 column of multi-column POST)
Bins of a site score for a motion series. POST has one value per project row, < 0 outside the site's level.
- VARIANT q5: 5 bins with edges at the quantiles of the level's genuine scores (all level rows when there is no TRUTH).
- VARIANT end10k: 10k-particle end bins + 3 equal middle bins.
Ties are broken at random (seed 7), and junk falls in the same score intervals. Writes OUT_STATES, one state per row: 0 outside the
level, 1..5 = bin from low to high score. Prints per bin the rows and junk share, and with TRUTH the mean / sd of ANGLE
(standardised over all genuine particles)."""
import sys, numpy as np
post, out, var = sys.argv[1:4]; Q = np.loadtxt(post, ndmin=2); Q = Q.T if Q.shape[1] > Q.shape[0] else Q; sc = Q[:, 0].copy(); NM = len(sc); lev = sc >= 0
junk = np.zeros(NM, bool); Z = None
if len(sys.argv) > 5:
    hd = [l.split()[1:] for l in open(sys.argv[4]) if l.startswith('#')][0]; T = np.genfromtxt(sys.argv[4], comments='#'); tr = T[:, 0].astype(int) - 1
    assert tr.max() + 1 == NM, 'truth rows %d != score rows %d' % (tr.max() + 1, NM); junk[tr] = T[:, hd.index('junk')] > 0; x = np.full(NM, np.nan); x[tr] = T[:, hd.index(sys.argv[5])]; G = ~junk & ~np.isnan(x)
    Z = (x - x[G].mean()) / x[G].std()
ref = lev & ~junk; s = sc + np.random.default_rng(7).uniform(0., 1e-9, NM); sg = np.sort(s[ref]); n = len(sg)
if var == 'q5': cuts = [sg[int(n * q) - 1] for q in (0.2, 0.4, 0.6, 0.8)]
elif var == 'end10k':
    assert n >= 25000, 'end10k needs >= 25,000 reference rows'; mid = (n - 20000) // 3
    cuts = [sg[9999], sg[10000 + mid - 1], sg[10000 + 2 * mid - 1], sg[n - 10001]]
else: sys.exit('VARIANT is q5 or end10k')
st = np.where(lev, np.searchsorted(np.array(cuts), s, side='left') + 1, 0); np.savetxt(out, st, fmt='%d')
print('%s: %d rows, %d in the level (%d reference rows for the edges), variant %s' % (post, NM, lev.sum(), n, var))
for k in range(1, 6):
    m = st == k; g = m & ~junk
    print('bin %d: rows %6d, genuine %6d, junk share %.3f, score %.4f-%.4f%s' % (k, m.sum(), g.sum(), (m & junk).sum() / max(m.sum(), 1), sc[m].min(), sc[m].max(),
          ', %s mean %+.2f sd %.2f' % (sys.argv[5], Z[g].mean(), Z[g].std()) if Z is not None else ''))
