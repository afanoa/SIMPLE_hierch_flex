#!/usr/bin/env python3
"""What tells same-conformer leaves apart on a null set: per state, mean pose error (refined vs truth), shift error,
defocus, refinement score, against random subsets of the same size (z). usage: leaf_diag.py TRUTH.txt REFINED.txt STATES.txt"""
import sys, numpy as np
def kv(line): return dict(t.split('=', 1) for t in line.split() if '=' in t)
def rot(e1, e2, e3):
    a, b, g = np.radians([e1, e2, e3]); ca, sa, cb, sb, cg, sg = np.cos(a), np.sin(a), np.cos(b), np.sin(b), np.cos(g), np.sin(g)
    return np.array([[cg * cb * ca - sg * sa, cg * cb * sa + sg * ca, -cg * sb], [-sg * cb * ca - cg * sa, -sg * cb * sa + cg * ca, sg * sb], [sb * ca, sb * sa, cb]])
def load(f, keys): return np.array([[float(kv(l).get(k, 0)) for k in keys] for l in open(f) if l.strip()])
T = load(sys.argv[1], ['e1', 'e2', 'e3', 'x', 'y', 'dfx', 'state']); R = load(sys.argv[2], ['e1', 'e2', 'e3', 'x', 'y', 'dfx', 'corr'])
st = np.loadtxt(sys.argv[3]).astype(int)[:len(T)]; n = min(len(T), len(R), len(st)); T, R, st = T[:n], R[:n], st[:n]
ang = np.array([np.degrees(np.arccos(np.clip((np.trace(rot(*t[:3]).T @ rot(*r[:3])) - 1) / 2, -1, 1))) for t, r in zip(T, R)])
dirc = np.array([np.degrees(np.arccos(np.clip(rot(*t[:3])[2] @ rot(*r[:3])[2], -1, 1))) for t, r in zip(T, R)])
sh = np.hypot(T[:, 3] - R[:, 3], T[:, 4] - R[:, 4]); feats = {'angerr': ang, 'direrr': dirc, 'shifterr': sh, 'defocus': R[:, 5], 'corr': R[:, 6], 'junk': (T[:, 6] == 0).astype(float)}
rng = np.random.default_rng(0); K = st.max(); print('%5s %7s' % ('state', 'n') + ''.join('%14s' % k for k in feats))
for s in range(1, K + 1):
    m = st == s; line = '%5d %7d' % (s, m.sum())
    for k, v in feats.items():
        rs = np.array([v[rng.choice(n, m.sum(), replace=False)].mean() for _ in range(50)]); z = (v[m].mean() - rs.mean()) / max(rs.std(), 1e-9)
        line += '%8.2f z%+5.1f' % (v[m].mean(), z)
    print(line)
print('whole set means: ' + ' '.join('%s %.2f' % (k, v.mean()) for k, v in feats.items()))
