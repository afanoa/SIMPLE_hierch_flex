#!/usr/bin/env python3
"""Do states differ by a COMMON pose error (a rotated / shifted copy of the map) rather than by conformation? Per state:
the mean rotation vector of the per-particle pose error in both composition orders (R_t^T R_r and R_r R_t^T; the map-frame
one carries a systematic rotation) and the mean 3-D shift error (image shift error rotated into the map frame), each as
|mean| (deg / px) with a z against random subsets of the same size. usage: leaf_pose_bias.py TRUTH REFINED[,REFINED2,...] STATES [FAMSUBSETS]"""
import sys, numpy as np
from scipy.spatial.transform import Rotation
def kv(line): return dict(t.split('=', 1) for t in line.split() if '=' in t)
def load(f, keys): return np.array([[float(kv(l).get(k, 0)) for k in keys] for l in open(f) if l.strip()])
def rot(e1, e2, e3):
    a, b, g = np.radians([e1, e2, e3]); ca, sa, cb, sb, cg, sg = np.cos(a), np.sin(a), np.cos(b), np.sin(b), np.cos(g), np.sin(g)
    return np.array([[cg * cb * ca - sg * sa, cg * cb * sa + sg * ca, -cg * sb], [-sg * cb * ca - cg * sa, -sg * cb * sa + cg * ca, sg * sb], [sb * ca, sb * sa, cb]])
K = ['e1', 'e2', 'e3', 'x', 'y', 'state']; T = load(sys.argv[1], K); refs = [load(f, K) for f in sys.argv[2].split(',')]
st = np.loadtxt(sys.argv[3]).astype(int); n = min(len(T), len(st), *[len(r) for r in refs]); T, st = T[:n], st[:n]; R = refs[0][:n].copy()
if len(sys.argv) > 4 and len(refs) > 1:   # per-family refined tables: row i takes the table of the family whose subset holds it
    for k, f in enumerate(sys.argv[4].split(',')):
        rows = np.loadtxt(f).astype(int) - 1; rows = rows[rows < n]; R[rows] = refs[k][:n][rows]
ok = (st > 0) & (T[:, 5] > 0); idx = np.where(ok)[0]
Rt = np.array([rot(*T[i, :3]) for i in idx]); Rr = np.array([rot(*R[i, :3]) for i in idx])
q1 = Rotation.from_matrix(np.einsum('nji,njk->nik', Rt, Rr)).as_rotvec(degrees=True); q2 = Rotation.from_matrix(np.einsum('nij,nkj->nik', Rr, Rt)).as_rotvec(degrees=True)
s_img = np.c_[R[idx, 3] - T[idx, 3], R[idx, 4] - T[idx, 4], np.zeros(len(idx))]
m1 = np.einsum('nji,nj->ni', Rt, s_img); m2 = np.einsum('nij,nj->ni', Rt, s_img)
lab = st[idx]; rng = np.random.default_rng(0); feats = {'rot Rt^T Rr': q1, 'rot Rr Rt^T': q2, 'shift R^T s': m1, 'shift R s': m2}
print('%d particles in %d states; whole: |mean| %s' % (len(idx), len(np.unique(lab)), ' '.join('%s %.3f' % (k, np.linalg.norm(v.mean(0))) for k, v in feats.items())))
print('%5s %6s' % ('state', 'n') + ''.join('%22s' % k for k in feats) + '   mean rotvec (Rt^T Rr, deg)')
for s in np.unique(lab):
    m = lab == s; line = '%5d %6d' % (s, m.sum())
    for k, v in feats.items():
        mm = np.linalg.norm(v[m].mean(0)); null = np.array([np.linalg.norm(v[rng.choice(len(v), m.sum(), replace=False)].mean(0)) for _ in range(100)])
        line += '%12.3f z%+7.1f' % (mm, (mm - null.mean()) / max(null.std(), 1e-12))
    print(line + '   (%+.2f %+.2f %+.2f)' % tuple(q1[m].mean(0)))
