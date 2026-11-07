#!/usr/bin/env python3
"""Do groups of particles differ in their consensus pose ERROR (truth known)? Per particle the map-frame rotation error between the
consensus pose and the true pose (rotation vector, deg) and the shift error (px); per group the mean error vector (a coherent mean
rotates / shifts the group's map, a reproducible nuisance difference), its spread, and the truth composition.
usage: pose_err_groups.py TRUE_ORIS CONS_ORIS LABELS[,LABELS2] [--gt gt_particles.txt]"""
import sys, argparse, numpy as np
from scipy.spatial.transform import Rotation
ap = argparse.ArgumentParser(); ap.add_argument('true'); ap.add_argument('cons'); ap.add_argument('labels'); ap.add_argument('--gt', default=''); a = ap.parse_args()
def kv(l): return dict(t.split('=', 1) for t in l.split() if '=' in t)
def rd(f): return np.array([[float(kv(l).get(k, 0)) for k in ('e1', 'e2', 'e3', 'x', 'y')] for l in open(f) if l.strip()])
T, C = rd(a.true), rd(a.cons); lab = sum(np.loadtxt(f, dtype=int) for f in a.labels.split(','))
def R(e): return Rotation.from_euler('ZYZ', e[:, :3], degrees=True)
Rt, Rc = R(T), R(C)
for nm, E in (('left (Rc Rt^-1)', Rc * Rt.inv()), ('right (Rt^-1 Rc)', Rt.inv() * Rc)):
    v = E.as_rotvec(degrees=True); ang = np.linalg.norm(v, axis=1)
    print('%s: rotation error median %.2f deg (p90 %.2f) over %d particles' % (nm, np.median(ang), np.quantile(ang, 0.9), len(ang)))
    for g in np.unique(lab[lab > 0]):
        m = lab == g; mu = v[m].mean(0); se = v[m].std(0).mean() / np.sqrt(m.sum())
        print('   group %d (n %d): mean error vector %.2f %.2f %.2f deg (|mean| %.2f, se %.3f), mean |error| %.2f' % (g, m.sum(), mu[0], mu[1], mu[2], np.linalg.norm(mu), se, ang[m].mean()))
ds = C[:, 3:5] - T[:, 3:5]; d3 = np.c_[ds, np.zeros(len(ds))]
for nm, V in (('R d', Rt.apply(d3)), ('R^T d', Rt.inv().apply(d3))):   # the image-frame shift error as a 3-D translation of the map (both conventions)
    print('3-D shift error (%s), px:' % nm)
    for g in np.unique(lab[lab > 0]): m = lab == g; mu = V[m].mean(0); print('   group %d (n %d): mean %.3f %.3f %.3f (|mean| %.3f, se %.3f)' % (g, m.sum(), mu[0], mu[1], mu[2], np.linalg.norm(mu), V[m].std(0).mean() / np.sqrt(m.sum())))
print('shift error (cons - true): median |d| %.2f px' % np.median(np.linalg.norm(ds, axis=1)))
for g in np.unique(lab[lab > 0]):
    m = lab == g; print('   group %d: mean shift error %.3f %.3f px, mean |d| %.2f' % (g, ds[m, 0].mean(), ds[m, 1].mean(), np.linalg.norm(ds[m], axis=1).mean()))
