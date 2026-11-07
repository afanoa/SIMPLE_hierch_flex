#!/usr/bin/env python3
# pose drift of a child node: incoming (n.simple = parent's stage-B poses) -> stage A -> stage B, against the parent's
# last-iteration spread (its stage-B dist / dist_inpl / shincarg over the child's rows). Usage:
#   pa_drift.py LABEL PARENT_B NODE_N NODE_A NODE_B SMPD ROWS_FILE (the node's rows in the parent's project)
import sys, numpy as np
KEYS = ('e1', 'e2', 'e3', 'x', 'y', 'state', 'dist', 'dist_inpl', 'shincarg', 'corr')

def read(fn):
    cols = {k: [] for k in KEYS}
    with open(fn) as f:
        for line in f:
            d = dict(t.split('=', 1) for t in line.split() if '=' in t)
            if 'e1' not in d:
                continue
            for k in KEYS:
                cols[k].append(float(d.get(k, 'nan')))
    return {k: np.array(v) for k, v in cols.items()}

def rotm(e1, e2, e3):
    a, b, c = np.deg2rad(e1), np.deg2rad(e2), np.deg2rad(e3)
    def rz(t):
        R = np.zeros((t.size, 3, 3)); R[:, 0, 0] = np.cos(t); R[:, 0, 1] = -np.sin(t); R[:, 1, 0] = np.sin(t); R[:, 1, 1] = np.cos(t); R[:, 2, 2] = 1
        return R
    def ry(t):
        R = np.zeros((t.size, 3, 3)); R[:, 0, 0] = np.cos(t); R[:, 0, 2] = np.sin(t); R[:, 2, 0] = -np.sin(t); R[:, 2, 2] = np.cos(t); R[:, 1, 1] = 1
        return R
    return rz(c) @ ry(b) @ rz(a)

def change(o1, o2, m, smpd):
    R1 = rotm(o1['e1'][m], o1['e2'][m], o1['e3'][m]); R2 = rotm(o2['e1'][m], o2['e2'][m], o2['e3'][m])
    D = R2 @ np.transpose(R1, (0, 2, 1))
    tr = np.trace(D, axis1=1, axis2=2)
    geo = np.rad2deg(np.arccos(np.clip((tr - 1) / 2, -1, 1)))
    ddir = np.rad2deg(np.arccos(np.clip(D[:, 2, 2], -1, 1)))
    tw = np.rad2deg(2 * np.arctan2(D[:, 1, 0] - D[:, 0, 1], 1 + tr))
    tw = (tw + 180) % 360 - 180
    dsh = smpd * np.hypot(o2['x'][m] - o1['x'][m], o2['y'][m] - o1['y'][m])
    return geo, ddir, np.abs(tw), dsh

def q(v):
    v = v[np.isfinite(v)]
    return 'med %6.2f p75 %6.2f p90 %6.2f p99 %7.2f' % tuple(np.percentile(v, [50, 75, 90, 99]))

lab, fp, fn, fa, fb, smpd = sys.argv[1:7]
smpd = float(smpd)
P, N, A, B = read(fp), read(fn), read(fa), read(fb)
m = N['state'] > 0
rows = np.loadtxt(sys.argv[7], dtype=int).ravel()
if P['state'].size != N['state'].size:
    P = {k: v[rows - 1] for k, v in P.items()}
print('## %s  n %d  (rows %d)  smpd %.4f' % (lab, m.sum(), m.size, smpd))
print('parent stage-B last-iteration spread over the child rows (dist = EMA of the per-iteration direction change, deg; dist_inpl deg; shincarg px):')
for k in ('dist', 'dist_inpl', 'shincarg'):
    print('  %-9s %s' % (k, q(P[k][m])))
print('  shincarg in A: %s' % q(smpd * P['shincarg'][m]))
print('incoming == parent stage B? max |de| %.3g' % np.nanmax(np.abs(np.r_[N['e1'][m] - P['e1'][m], N['e2'][m] - P['e2'][m], N['e3'][m] - P['e3'][m]])))
for name, o1, o2 in (('in->A', N, A), ('A->B', A, B), ('in->B', N, B)):
    geo, ddir, tw, dsh = change(o1, o2, m, smpd)
    print('%-6s geodesic %s' % (name, q(geo)))
    print('%-6s direction %s' % ('', q(ddir)))
    print('%-6s in-plane  %s   frac>10 %.3f  >30 %.3f  >90 %.3f' % ('', q(tw), (tw > 10).mean(), (tw > 30).mean(), (tw > 90).mean()))
    print('%-6s shift(A)  %s   frac>5A %.3f  >10A %.3f' % ('', q(dsh), (dsh > 5).mean(), (dsh > 10).mean()))
print('child stage-B last-iteration: dist %s | dist_inpl %s | shincarg(A) %s' % (q(B['dist'][m]), q(B['dist_inpl'][m]), q(smpd * B['shincarg'][m])))
