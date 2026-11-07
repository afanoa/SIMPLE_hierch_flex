#!/usr/bin/env python3
"""usage: motion_bins3.py SCORE OUT_STATES KIND TRUTH   (KIND cahra | cdisc)
motion_bins2.py's 5 genuine-quantile bins of SCORE (one value per row, < 0 outside the level; ties broken at random, seed 7; junk in
the same score intervals), with the truth composition per bin:
- cahra: the open fraction;
- cdisc: A-present and closed fractions among genuine, and the junk share.
Writes OUT_STATES (0 outside, 1..5 low to high score)."""
import sys, numpy as np
SC, OUT, KIND, TF = sys.argv[1:5]; sc = np.loadtxt(SC); NM = len(sc); lev = sc >= 0
if KIND == 'cahra':
    sys.path.insert(0, '/mnt/beegfs/elmlund/afan/cahra_13991'); import cahra_score as cs
    junk = np.zeros(NM, bool); comp = {'open': (np.array(cs.read_star(TF)['conformation'], int)[:NM] == 1)}
else:
    hd = [l.split()[1:] for l in open(TF) if l.startswith('#')][0]; M = np.genfromtxt(TF, comments='#'); tr = M[:, 0].astype(int) - 1
    assert tr.max() + 1 == NM, 'truth rows %d != score rows %d' % (tr.max() + 1, NM)
    def col(c): v = np.zeros(NM); v[tr] = M[:, hd.index(c)]; return v
    junk = col('junk') > 0; comp = {'A': col('A_domain') > 0, 'closed': col('C_closed') > 0}
ref = lev & ~junk; s = sc + np.random.default_rng(7).uniform(0., 1e-9, NM); sg = np.sort(s[ref]); n = len(sg)
cuts = [sg[int(n * q) - 1] for q in (0.2, 0.4, 0.6, 0.8)]; st = np.where(lev, np.searchsorted(np.array(cuts), s, side='left') + 1, 0); np.savetxt(OUT, st, fmt='%d')
print('%s: %d rows, %d in the level, %d reference rows' % (SC, NM, lev.sum(), n))
for k in range(1, 6):
    m = st == k; g = m & ~junk
    print('bin %d: rows %6d, genuine %6d, junk share %.3f, %s' % (k, m.sum(), g.sum(), (m & junk).sum() / max(m.sum(), 1), ', '.join('%s %.3f' % (c, v[g].mean()) for c, v in comp.items())))
