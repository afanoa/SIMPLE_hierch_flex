#!/usr/bin/env python3
"""IgM arms truth vs leaf-map factor codes: per component, the size-weighted correlation of its leaf code with each arm angle's leaf mean
(genuine rows); best arm angle per component, and how many distinct arm angles the components cover. usage: igm_codes_score.py LF_DIR TREE_DIR TRUTH"""
import sys, numpy as np
lf, T, tr = sys.argv[1:4]
cols = open(tr).readline().lstrip('#').split(); t = np.loadtxt(tr, comments='#', ndmin=2); C = {c: t[:, i] for i, c in enumerate(cols)}
rows = C['row'].astype(int); junk = C['junk'] > 0; ang = [c for c in cols if c.startswith('a')]
lab = np.loadtxt(T + '/hetero_tree_labels.txt', dtype=int, ndmin=1)[rows - 1]
L = [l.split() for l in open(lf + '/codes.txt') if not l.startswith('#')]; names = [l[0] for l in L]; r = (len(L[0]) - 2) // 2
code = np.array([[float(x) for x in l[2:2 + r]] for l in L])
leafnames = [l.split()[0] for l in open(T + '/hetero_tree_leaves.txt') if l.strip()]
M = np.full((len(names), len(ang)), np.nan); n = np.zeros(len(names))
for i, nm in enumerate(names):
    k = leafnames.index(nm) + 1; g = (lab == k) & ~junk; n[i] = g.sum()
    if g.sum(): M[i] = [C[a][g].mean() for a in ang]
ok = n > 0; w = n[ok] / n[ok].sum()
def wcorr(x, y): x = x - (w * x).sum(); y = y - (w * y).sum(); return (w * x * y).sum() / np.sqrt((w * x * x).sum() * (w * y * y).sum())
R = np.array([[wcorr(code[ok, j], M[ok, q]) for q in range(len(ang))] for j in range(r)])
best = {}
for j in range(r):
    o = np.argsort(-np.abs(R[j])); q = o[0]; best.setdefault(ang[q], []).append((abs(R[j, q]), j + 1))
    print('c%-2d best %s r %+.2f (second %s %+.2f)' % (j + 1, ang[q], R[j, q], ang[o[1]], R[j, o[1]]))
print('arm angles covered (best |r| >= 0.5): %s' % ' '.join('%s:c%s(%.2f)' % (a, max(v)[1], max(v)[0]) for a, v in sorted(best.items()) if max(v)[0] >= 0.5))
# leaf-level R^2 of each arm angle from all component codes vs the leaves' own eta^2 ceiling at leaf level
for q, a in enumerate(ang):
    y = M[ok, q]; X = np.c_[np.ones(ok.sum()), code[ok]]; b = np.linalg.lstsq(X * np.sqrt(w)[:, None], y * np.sqrt(w), rcond=None)[0]
    res = y - X @ b; r2 = 1 - (w * res ** 2).sum() / (w * (y - (w * y).sum()) ** 2).sum(); print('%s leaf-mean R2 from codes %.2f' % (a, r2), end='; ')
print()
