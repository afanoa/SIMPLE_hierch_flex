#!/usr/bin/env python3
"""Pose error of a ptcl3D table against the toy truth (all.txt): median/p90 angle over conformer particles, shift error.
usage: pose_error.py TABLE TRUTH [--states states.txt (restrict to state>0)] [--tag name]"""
import sys, argparse, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('table'); ap.add_argument('truth'); ap.add_argument('--states', default=''); ap.add_argument('--tag', default='poses'); a = ap.parse_args()
def rd(fn, keys): return np.array([[float(dict(t.split('=', 1) for t in l.split() if '=' in t).get(k, 0)) for k in keys] for l in open(fn) if l.strip() and not l.startswith('#')])
def rotmat(ang, axis):
    r = np.radians(ang); c, s = np.cos(r), np.sin(r); R = np.eye(3)
    if axis == 3: R[0, 0], R[0, 1], R[1, 0], R[1, 1] = c, -s, s, c
    else: R[0, 0], R[0, 2], R[2, 0], R[2, 2] = c, -s, s, c
    return R
def e2m(e): return rotmat(e[2], 3) @ rotmat(e[1], 2) @ rotmat(e[0], 3)
r = rd(a.table, ['e1', 'e2', 'e3', 'x', 'y', 'state']); t = rd(a.truth, ['e1', 'e2', 'e3', 'x', 'y', 'state']); n = min(len(r), len(t))
good = t[:n, 5] > 0
if a.states: st = np.loadtxt(a.states).astype(int)[:n]; good &= st > 0
ang = np.array([np.degrees(np.arccos(np.clip((np.trace(e2m(r[i, :3]) @ e2m(t[i, :3]).T) - 1) / 2, -1, 1))) for i in range(n)])
print('%s: angle median %.2f deg (p90 %.2f) over %d conformer particles; shift error median %.2f px' % (a.tag, np.median(ang[good]), np.percentile(ang[good], 90), good.sum(), np.median(np.hypot(r[:n, 3] - t[:n, 3], r[:n, 4] - t[:n, 4])[good])))
