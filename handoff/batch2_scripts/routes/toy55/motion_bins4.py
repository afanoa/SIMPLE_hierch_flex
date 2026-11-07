#!/usr/bin/env python3
"""usage: motion_bins4.py SCORE OUT_STATES KIND TRUTH EXCLUDE TREE_DIR [ANGLE]   (KIND cdisc | cahra | arms; EXCLUDE tree | truth | none)
5 quantile bins of SCORE (one value per row, < 0 outside the level), cut as a tree node would cut them, with junk out of every bin:
- tree: rows the tree calls junk (no leaf, or a leaf TREE_DIR/hetero_tree_leaves.txt names 'junk') get state 0, and the edges are
  quantiles of the remaining level rows, so no truth enters;
- truth: the truth junk gets state 0 and the edges are the genuine quantiles (an oracle);
- none: every level row is binned.
Ties are broken at random (seed 7). Writes OUT_STATES (0 out, 1..5 low to high score) and prints per bin the rows, the truth junk
share left and the truth composition: cdisc the A-present and closed fractions among genuine, cahra the open fraction, arms the
mean of ANGLE ('ANGLE mean', standardised over all genuine particles) or, without ANGLE, the mean of every arm angle."""
import sys, numpy as np
SC, OUT, KIND, TF, EX, D = sys.argv[1:7]; ANG = sys.argv[7] if len(sys.argv) > 7 else None
sc = np.loadtxt(SC); NM = len(sc); lev = sc >= 0
if KIND == 'cahra':
    sys.path.insert(0, '/mnt/beegfs/elmlund/afan/cahra_13991'); import cahra_score as cs
    junk = np.zeros(NM, bool); comp = {'open': (np.array(cs.read_star(TF)['conformation'], int)[:NM] == 1).astype(float)}
else:
    hd = [l.split()[1:] for l in open(TF) if l.startswith('#')][0]; M = np.genfromtxt(TF, comments='#'); tr = M[:, 0].astype(int) - 1
    assert tr.max() + 1 == NM, 'truth rows %d != score rows %d' % (tr.max() + 1, NM)
    def col(c): v = np.full(NM, np.nan); v[tr] = M[:, hd.index(c)]; return v
    junk = col('junk') > 0; comp = {}
    if KIND == 'cdisc': comp = {'A': (col('A_domain') > 0).astype(float), 'closed': (col('C_closed') > 0).astype(float)}
    else:
        for a in hd[3:]: x = col(a); g = ~junk & ~np.isnan(x); comp[a] = (x - x[g].mean()) / x[g].std()
L = np.loadtxt(D + '/hetero_tree_labels.txt', dtype=int); assert len(L) == NM, 'tree rows %d != score rows %d' % (len(L), NM)
LV = [l.split() for l in open(D + '/hetero_tree_leaves.txt') if l.strip()]; jm = [k + 1 for k, w in enumerate(LV) if w[1] == 'junk']
out = {'tree': (L == 0) | np.isin(L, jm), 'truth': junk, 'none': np.zeros(NM, bool)}[EX]
ref = lev & ~out; s = sc + np.random.default_rng(7).uniform(0., 1e-9, NM); sg = np.sort(s[ref]); n = len(sg)
cuts = [sg[int(n * q) - 1] for q in (0.2, 0.4, 0.6, 0.8)]; st = np.where(ref, np.searchsorted(np.array(cuts), s, side='left') + 1, 0); np.savetxt(OUT, st, fmt='%d')
print('%s: %d rows, %d in the level, %d excluded (%s; truth junk share among them %.3f), %d binned' % (SC, NM, lev.sum(), (lev & out).sum(), EX,
      (lev & out & junk).sum() / max((lev & out).sum(), 1), n))
for k in range(1, 6):
    m = st == k; g = m & ~junk
    if KIND != 'arms': txt = ', '.join('%s %.3f' % (c, v[g].mean()) for c, v in comp.items())
    elif ANG: txt = '%s mean %+.2f sd %.2f' % (ANG, comp[ANG][g].mean(), comp[ANG][g].std())
    else: txt = 'angles ' + ' '.join('%s %+.2f' % (a, np.nanmean(v[g])) for a, v in comp.items())
    print('bin %d: rows %6d, genuine %6d, junk share %.3f, score %.4f-%.4f, %s' % (k, m.sum(), g.sum(), (m & junk).sum() / max(m.sum(), 1), sc[m].min(), sc[m].max(), txt))
