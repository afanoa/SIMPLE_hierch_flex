#!/usr/bin/env python3
"""What particle property separates two groups (e.g. the two families of a tree)? Univariate AUC of every per-particle
attribute (view, in-plane angle, defocus, astigmatism, pose error, shift error, refinement score, conformer) and the
cross-validated accuracy of a logistic model on all of them. usage: split_discrim.py TRUTH REFINED ROWS1 ROWS2 [GT]  (ROWS files: 1-based rows, or 'LABELS:k' = rows where a label file equals k)"""
import sys, numpy as np
def kv(l): return dict(t.split('=', 1) for t in l.split() if '=' in t)
def load(f, keys): return np.array([[float(kv(l).get(k, 0)) for k in keys] for l in open(f) if l.strip()])
def rot(e1, e2, e3):
    a, b, g = np.radians([e1, e2, e3]); ca, sa, cb, sb, cg, sg = np.cos(a), np.sin(a), np.cos(b), np.sin(b), np.cos(g), np.sin(g)
    return np.array([[cg * cb * ca - sg * sa, cg * cb * sa + sg * ca, -cg * sb], [-sg * cb * ca - cg * sa, -sg * cb * sa + cg * ca, sg * sb], [sb * ca, sb * sa, cb]])
T = load(sys.argv[1], ['e1', 'e2', 'e3', 'x', 'y', 'dfx', 'dfy', 'angast']); Rr = load(sys.argv[2], ['e1', 'e2', 'e3', 'x', 'y', 'corr'])
def rows_of(spec):
    if ':' in spec: f, k = spec.rsplit(':', 1); lab = np.loadtxt(f).astype(int); return np.where(lab == int(k))[0]
    return np.loadtxt(spec).astype(int) - 1
r1 = rows_of(sys.argv[3]); r2 = rows_of(sys.argv[4]); rows = np.r_[r1, r2]; y = np.r_[np.zeros(len(r1)), np.ones(len(r2))]
conf = np.loadtxt(sys.argv[5], comments='#')[:, 1].astype(int) if len(sys.argv) > 5 else None
ok = np.ones(len(rows), bool) if conf is None else conf[rows] > 0
if conf is not None and len(sys.argv) > 6: ok &= conf[rows] == int(sys.argv[6])   # restrict to one conformer
rows, y = rows[ok], y[ok]
F = {}
Rt = np.array([rot(*T[i, :3]) for i in rows]); Rf = np.array([rot(*Rr[i, :3]) for i in rows])
F['view x'], F['view y'], F['view z'] = Rt[:, 2, 0], Rt[:, 2, 1], Rt[:, 2, 2]
F['|view z|'] = np.abs(Rt[:, 2, 2]); F['in-plane e3'] = np.cos(np.radians(T[rows, 2]))
F['defocus'] = 0.5 * (T[rows, 5] + T[rows, 6]); F['astigmatism'] = np.abs(T[rows, 5] - T[rows, 6])
F['pose error'] = np.degrees(np.arccos(np.clip((np.einsum('nij,nij->n', Rt, Rf) - 1) / 2, -1, 1)))
F['direction error'] = np.degrees(np.arccos(np.clip(np.einsum('ni,ni->n', Rt[:, 2], Rf[:, 2]), -1, 1)))
F['shift error'] = np.hypot(T[rows, 3] - Rr[rows, 3], T[rows, 4] - Rr[rows, 4]); F['refinement corr'] = Rr[rows, 5]
from scipy.spatial.transform import Rotation
q = Rotation.from_matrix(np.einsum('nji,njk->nik', Rt, Rf)).as_rotvec(degrees=True)
for k in range(3): F['pose err rot %d' % k] = q[:, k]
sm = np.einsum('nji,nj->ni', Rt, np.c_[Rr[rows, 3] - T[rows, 3], Rr[rows, 4] - T[rows, 4], np.zeros(len(rows))])
for k in range(3): F['shift err map %d' % k] = sm[:, k]
F['|view x|'], F['|view y|'] = np.abs(Rt[:, 2, 0]), np.abs(Rt[:, 2, 1]); F['in-plane error'] = np.abs((T[rows, 2] - Rr[rows, 2] + 180) % 360 - 180)
F['row index'] = rows.astype(float)
if conf is not None: F['conformer'] = conf[rows].astype(float)
def auc(v, y):
    o = np.argsort(v); rk = np.empty(len(v)); rk[o] = np.arange(1, len(v) + 1); n1 = y.sum(); n0 = len(y) - n1
    return (rk[y == 1].sum() - n1 * (n1 + 1) / 2) / max(n0 * n1, 1)
print('%d vs %d particles (non-junk)' % ((y == 0).sum(), (y == 1).sum()))
for k, v in sorted(F.items(), key=lambda kv: -abs(auc(kv[1], y) - 0.5)): print('  %-17s AUC %.3f   group means %.4g vs %.4g' % (k, auc(v, y), v[y == 0].mean(), v[y == 1].mean()))
X = np.column_stack([(v - v.mean()) / max(v.std(), 1e-12) for k, v in F.items() if k != 'conformer']); X = np.c_[np.ones(len(X)), X]
rng = np.random.default_rng(0); fold = rng.integers(0, 5, len(y)); acc = []
for f in range(5):
    tr, te = fold != f, fold == f; w = np.zeros(X.shape[1])
    for _ in range(300): p = 1 / (1 + np.exp(-X[tr] @ w)); w += 0.5 * X[tr].T @ (y[tr] - p) / tr.sum()
    acc.append(np.mean(((X[te] @ w) > 0) == (y[te] == 1)))
print('  logistic model on all attributes (no conformer): 5-fold accuracy %.3f (majority %.3f)' % (np.mean(acc), max(y.mean(), 1 - y.mean())))
