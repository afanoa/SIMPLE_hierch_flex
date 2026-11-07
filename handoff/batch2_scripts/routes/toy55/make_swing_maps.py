#!/usr/bin/env python3
"""Synthetic swing of the closed state: the mobile domain of the deposited closed map (variants/domain_mask.mrc) rotated
rigidly about its hinge (centroid of the domain/rest interface) toward its position in the open conformers, the domain tip
travelling TRAVEL A in NSTEPS equal steps (step 0 = the deposited closed map). Writes OUT/conf<ID0+k>_swing_m160.mrc and
OUT/truth_swing.txt (step, angle, tip travel, domain-centroid travel; A). No particles: the maps only.
usage: make_swing_maps.py CLOSED DOMAIN_MASK OUT --open MAP[,MAP...] [--travel 30 --nsteps 31 --id0 200 --smpd 1.464]"""
import sys, os, argparse, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from model_class import read_mrc, write_mrc
from scipy.ndimage import gaussian_filter, binary_dilation, affine_transform
ap = argparse.ArgumentParser(); ap.add_argument('closed'); ap.add_argument('domain'); ap.add_argument('out'); ap.add_argument('--open', required=True)
ap.add_argument('--travel', type=float, default=30.); ap.add_argument('--nsteps', type=int, default=31); ap.add_argument('--id0', type=int, default=200)
ap.add_argument('--smpd', type=float, default=1.464); a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
C = read_mrc(a.closed).astype(np.float32); Dm = read_mrc(a.domain) > 0.5; A = np.clip(gaussian_filter(Dm.astype(np.float32), 1.0), 0, 1)
Cl = gaussian_filter(C, 1.5); dens = Cl > Cl.mean() + 1.5 * Cl.std(); dom = dens & Dm; rest = dens & ~Dm
idx = np.indices(C.shape).reshape(3, -1).T.astype(float)
def cen(m, w=None): w = (np.ones(m.sum()) if w is None else w[m.ravel()]); p = idx[m.ravel()]; return (p * w[:, None]).sum(0) / w.sum()
M = cen(rest, np.maximum(C.ravel(), 0)); dd = np.linalg.norm(idx[dom.ravel()] - M, axis=1)   # hinge = the domain's base: its voxels nearest the rest's centre
base_ = np.zeros(dom.sum(), bool); base_[dd <= np.percentile(dd, 10)] = True; h = idx[dom.ravel()][base_].mean(0); g = cen(dom, np.maximum(C.ravel(), 0))
O = np.mean([read_mrc(f).astype(np.float32) for f in a.open.split(',')], 0); near = binary_dilation(dom, iterations=8) & ~rest
d = gaussian_filter(O - C, 1.5); pos = near & (d > d[near].mean() + 2 * d[near].std()); o = cen(pos, np.maximum(d.ravel(), 0)) if pos.sum() > 20 else g + np.array([0., 0., 10.])
r = g - h; rh = r / np.linalg.norm(r); t = o - g; t = t - t.dot(rh) * rh; th = t / max(np.linalg.norm(t), 1e-9); u = np.cross(rh, th); u /= np.linalg.norm(u)
dist = np.linalg.norm(idx[dom.ravel()] - h, axis=1); rtip = np.percentile(dist, 95) * a.smpd; tmax = a.travel / rtip
def rot(u, t):
    K = np.array([[0, -u[2], u[1]], [u[2], 0, -u[0]], [-u[1], u[0], 0]]); return np.eye(3) + np.sin(t) * K + (1 - np.cos(t)) * K @ K
P = C * A; base = C * (1 - A); lines = []
for k in range(a.nsteps):
    tk = tmax * k / (a.nsteps - 1); R = rot(u, tk)
    if (R @ (g - h)).dot(th) < 0 and k: R = rot(u, -tk)   # rotate toward the open position
    Rk = affine_transform(P, R.T, offset=h - R.T @ h, order=3, mode='constant') if k else P
    write_mrc('%s/conf%d_swing_m160.mrc' % (a.out, a.id0 + k), (base + Rk).astype(np.float32), a.smpd)
    gk = R @ (g - h) + h; lines.append('%d %d %.3f %.2f %.2f' % (k, a.id0 + k, np.degrees(tk), rtip * tk, np.linalg.norm(gk - g) * a.smpd))
with open(a.out + '/truth_swing.txt', 'w') as f:
    f.write('# step conformer angle_deg tip_travel_A centroid_travel_A (hinge %s, axis %s, tip radius %.1f A, domain %.0f nm^3)\n' % (np.round(h, 1), np.round(u, 3), rtip, dom.sum() * a.smpd ** 3 / 1000))
    f.write('\n'.join(lines) + '\n')
print('swing: domain %.0f nm^3, hinge %s, tip radius %.1f A -> %.1f deg for %.0f A of tip travel in %d steps; centroid travel %.1f A; open-direction offset %.1f A' % (
    dom.sum() * a.smpd ** 3 / 1000, np.round(h, 1), rtip, np.degrees(tmax), a.travel, a.nsteps, float(lines[-1].split()[4]), np.linalg.norm(o - g) * a.smpd))
