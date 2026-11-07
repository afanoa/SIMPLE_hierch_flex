#!/usr/bin/env python3
"""usage: promote_score.py TRUTH NAME:POST:ANGLE ..
Scores site posteriors (one value per project row, < 0 outside the level) the way paths.py scores kept sites, on genuine particles,
with each angle standardised over all genuine particles:
- abs(Spearman) with the target ANGLE, and the best-matching angle;
- 5 genuine-quantile bins and 10k-end bins along the posterior, with junk counted in the same score intervals. Per variant: the bin
  Spearman, the span (sd), the span / the oracle's span at the end bins' genuine sizes, and the max / min step;
- tails: the top and bottom 10k by the posterior, as a ratio to the oracle of the same genuine count (as in extremes.py)."""
import sys, numpy as np
TF = sys.argv[1]; ENT = [a.split(':') for a in sys.argv[2:]]
hd = [l.split()[1:] for l in open(TF) if l.startswith('#')][0]; T = np.genfromtxt(TF, comments='#'); tr = T[:, 0].astype(int) - 1; NM = tr.max() + 1
ang = hd[3:]; X = np.full((NM, len(ang)), np.nan); X[tr] = T[:, 3:3 + len(ang)]
junk = np.zeros(NM, bool); junk[tr] = T[:, hd.index('junk')] > 0; G = ~junk & ~np.isnan(X[:, 0]); Z = np.nan_to_num((X - X[G].mean(0)) / X[G].std(0))
EH = [np.cumsum(np.sort(Z[G, c])[::-1]) / np.arange(1, G.sum() + 1) for c in range(len(ang))]
EL = [np.cumsum(np.sort(-Z[G, c])[::-1]) / np.arange(1, G.sum() + 1) for c in range(len(ang))]
def rank(v): r = np.empty(len(v)); r[np.argsort(v, kind='stable')] = np.arange(len(v)); return r
def spear(x, y): return float(np.corrcoef(rank(np.asarray(x, float)), rank(np.asarray(y, float)))[0, 1])
for nm, f, a in ENT:
    sc = np.array(open(f).read().split(), float)[:NM]; lev = sc >= 0; gl = lev & G; c = ang.index(a)
    rh = np.array([spear(sc[gl], Z[gl, j]) for j in range(len(ang))]); sg = 1. if rh[c] > 0 else -1.; jb = int(np.argmax(np.abs(rh)))
    s = sc + np.random.default_rng(7).uniform(0., 1e-9, NM); sgs = np.sort(s[gl]); n = len(sgs); line = []
    for var in ('q5', 'end10k'):
        if var == 'q5': cuts = [sgs[int(n * q) - 1] for q in (0.2, 0.4, 0.6, 0.8)]
        else: mid = (n - 20000) // 3; cuts = [sgs[9999], sgs[10000 + mid - 1], sgs[10000 + 2 * mid - 1], sgs[n - 10001]]
        b = np.searchsorted(np.array(cuts), s, side='left'); mu, ng, js = [], [], []
        for k in range(5):
            m = lev & (b == k); g = m & G; mu.append(sg * Z[g, c].mean()); ng.append(g.sum()); js.append((m & junk).sum() / max(m.sum(), 1))
        mu = np.array(mu); span = mu[-1] - mu[0]; osp = (EH[c] if sg > 0 else EL[c])[ng[-1] - 1] + (EL[c] if sg > 0 else EH[c])[ng[0] - 1]; st = np.diff(mu)
        line.append('%s: bins Spearman %.2f, span %.2f, ratio %.2f, steps %.2f-%.2f | means %s | junk %s' % (var, spear(np.arange(5), mu), span, span / osp, st.min(), st.max(),
                    ' '.join('%+.2f' % x for x in mu), ' '.join('%.2f' % x for x in js)))
    o = np.where(lev)[0][np.argsort(s[lev], kind='stable')]; tails = []
    for end, rows in (('high', o[-10000:]), ('low', o[:10000])):
        g = rows[G[rows]]; e = (1 if end == 'high' else -1) * sg * Z[g, c].mean(); want = 'high' if (end == 'high') == (sg > 0) else 'low'
        tails.append('%s-score 10k -> %s end ratio %.2f (genuine %d, junk %.2f)' % (end, want, e / (EH[c] if want == 'high' else EL[c])[len(g) - 1], len(g), 1 - len(g) / len(rows)))
    print('== %s (%s): target %s abs rho %.3f (sign %+d); best angle %s %.3f; level %d rows, %d genuine' % (nm, f.split('/')[-2] + '/' + f.split('/')[-1], a, abs(rh[c]), sg, ang[jb], abs(rh[jb]), lev.sum(), gl.sum()))
    for l in line + tails: print('   ' + l)
