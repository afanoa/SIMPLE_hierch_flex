#!/usr/bin/env python3
"""Is the neighbouring-bin matching right where it can be checked? Per sub-class truth = member mean of the arm angles /
site occupancies, centred within its bin (the contrasts are bin-centred too). For every edge: correlation of the truth
contrasts of the MATCHED sub-class pairs vs all other pairs of the same two bins. usage: propagate_check.py PROP_DIR CLSEXP_DIR GT TRUTH MODE"""
import sys, numpy as np
d, clsexp, gtf, trf, mode = sys.argv[1:6]
conf = {int(w[0]): int(w[1]) for w in (l.split() for l in open(gtf)) if w and not w[0].startswith('#')}
T = {}
for w in (l.split() for l in open(trf)):
    if not w or w[0].startswith('#'): continue
    v = np.array([float(x) for x in w[1:]]); T[int(w[0])] = v[0::4] if mode.startswith('tethered') else v
L = np.loadtxt(clsexp + '/cls_expansion_labels.txt', comments='#').astype(int); sub = L[:, 2]; N = sub.max(); dim = len(next(iter(T.values())))
tp = np.array([T[conf[p]] if conf.get(p) in T else np.full(dim, np.nan) for p in L[:, 0]]); ts = np.full((N, dim), np.nan)
for g in range(1, N + 1):
    m = sub == g
    if m.any(): ts[g - 1] = np.nanmean(tp[m], 0)
E = np.load(d + '/edges.npz'); S = np.loadtxt(d + '/subclass_state.txt', dtype=int); par = None
import glob
M = np.load(glob.glob(d + '/../stitch/avgs_meta.npz')[0]); par = M['parent'].astype(int)
tc = ts.copy()
for b in np.unique(par):
    m = par == b; tc[m] -= np.nanmean(ts[m], 0)
def c(u, v): u, v = u - u.mean(), v - v.mean(); return float((u * v).sum() / (np.sqrt((u * u).sum() * (v * v).sum()) + 1e-30))
matched, other = [], []
for ie, (i, j) in enumerate(E['pairs']):
    mi, mj = np.where(par == i + 1)[0], np.where(par == j + 1)[0]; mt = {(u, v) for (u, v), e in zip(E['matches'], E['edge_of_match']) if e == ie}
    for u in mi:
        for v in mj:
            if np.isfinite(tc[u]).all() and np.isfinite(tc[v]).all() and np.abs(tc[u]).sum() > 0 and np.abs(tc[v]).sum() > 0: (matched if (u, v) in mt else other).append(c(tc[u], tc[v]))
print('truth-contrast agreement of matched sub-class pairs: mean %.3f (n=%d) vs other pairs of the same bins %.3f (n=%d); fraction of matched pairs with agreement > 0.5: %.2f' % (np.mean(matched), len(matched), np.mean(other), len(other), np.mean(np.array(matched) > 0.5)))
hi = E['conf'] > np.median(E['conf']); mh = [c_ for c_, e in zip(matched, [ie for ie in range(len(E['pairs'])) for _ in range(1)]) ] if False else None
