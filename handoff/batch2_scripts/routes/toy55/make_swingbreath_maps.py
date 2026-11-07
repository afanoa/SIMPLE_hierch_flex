#!/usr/bin/env python3
"""Two motions on one particle: the closed state's mobile domain swung rigidly about its hinge (as make_swing_maps.py, tip travel
TRAVEL A in NK steps) and then the whole periphery breathing toward the density centroid (as make_breath_maps.py, alpha up to
AMAX in NJ steps). One map per (k, j): OUT/conf<ID0 + 10 k + j>_sb_m160.mrc; OUT/truth_grid.txt (conformer k j tip_travel_A alpha).
usage: make_swingbreath_maps.py CLOSED DOMAIN_MASK OUT --open MAP[,MAP...] [--travel 30 --nk 16 --amax 0.10 --nj 6 --id0 400]"""
import sys, os, argparse, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from model_class import read_mrc, write_mrc
from scipy.ndimage import gaussian_filter, binary_dilation, affine_transform, map_coordinates
ap = argparse.ArgumentParser(); ap.add_argument('closed'); ap.add_argument('domain'); ap.add_argument('out'); ap.add_argument('--open', required=True)
ap.add_argument('--travel', type=float, default=30.); ap.add_argument('--nk', type=int, default=16); ap.add_argument('--amax', type=float, default=0.10); ap.add_argument('--nj', type=int, default=6)
ap.add_argument('--r0', type=float, default=25.); ap.add_argument('--r1', type=float, default=55.); ap.add_argument('--id0', type=int, default=400); ap.add_argument('--smpd', type=float, default=1.464)
a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
C = read_mrc(a.closed).astype(np.float32); Dm = read_mrc(a.domain) > 0.5; A = np.clip(gaussian_filter(Dm.astype(np.float32), 1.0), 0, 1)
Cl = gaussian_filter(C, 1.5); dens = Cl > Cl.mean() + 1.5 * Cl.std(); dom = dens & Dm; rest = dens & ~Dm
idx = np.indices(C.shape).reshape(3, -1).T.astype(float)
def cen(m, w=None): w = (np.ones(m.sum()) if w is None else w[m.ravel()]); p = idx[m.ravel()]; return (p * w[:, None]).sum(0) / w.sum()
Mc = cen(rest, np.maximum(C.ravel(), 0)); dd = np.linalg.norm(idx[dom.ravel()] - Mc, axis=1); h = idx[dom.ravel()][dd <= np.percentile(dd, 10)].mean(0); g = cen(dom, np.maximum(C.ravel(), 0))
O = np.mean([read_mrc(f).astype(np.float32) for f in a.open.split(',')], 0); near = binary_dilation(dom, iterations=8) & ~rest
d = gaussian_filter(O - C, 1.5); pos = near & (d > d[near].mean() + 2 * d[near].std()); o = cen(pos, np.maximum(d.ravel(), 0)) if pos.sum() > 20 else g + np.array([0., 0., 10.])
r = g - h; rh = r / np.linalg.norm(r); t = o - g; t = t - t.dot(rh) * rh; th = t / max(np.linalg.norm(t), 1e-9); u = np.cross(rh, th); u /= np.linalg.norm(u)
rtip = np.percentile(np.linalg.norm(idx[dom.ravel()] - h, axis=1), 95) * a.smpd; tmax = a.travel / rtip
def rot(u, t):
    K = np.array([[0, -u[2], u[1]], [u[2], 0, -u[0]], [-u[1], u[0], 0]]); return np.eye(3) + np.sin(t) * K + (1 - np.cos(t)) * K @ K
P = C * A; base = C * (1 - A)
cw = np.maximum(C, 0) * dens; ii = np.indices(C.shape).astype(np.float32); c0 = np.array([(ii[i] * cw).sum() / cw.sum() for i in range(3)])
X = ii - c0[:, None, None, None]; rr = np.sqrt((X ** 2).sum(0)) * a.smpd; s = np.clip((rr - a.r0) / (a.r1 - a.r0), 0, 1); gb = s * s * (3 - 2 * s); lines = []
for k in range(a.nk):
    tk = tmax * k / (a.nk - 1); R = rot(u, tk)
    if k and (R @ (g - h)).dot(th) < 0: R = rot(u, -tk)
    S = base + (affine_transform(P, R.T, offset=h - R.T @ h, order=3, mode='constant') if k else P)
    for j in range(a.nj):
        al = a.amax * j / (a.nj - 1)
        V = map_coordinates(S, (c0[:, None, None, None] + X / np.maximum(1 - al * gb, 1e-3)).reshape(3, -1), order=3, mode='constant').reshape(S.shape) if j else S
        cid = a.id0 + 10 * k + j; write_mrc('%s/conf%d_sb_m160.mrc' % (a.out, cid), V.astype(np.float32), a.smpd); lines.append('%d %d %d %.2f %.4f' % (cid, k, j, rtip * tk, al))
with open(a.out + '/truth_grid.txt', 'w') as f:
    f.write('# conformer k j tip_travel_A alpha (swing: hinge %s axis %s tip radius %.1f A; breathing: centre %s R0 %.0f R1 %.0f A)\n' % (np.round(h, 1), np.round(u, 3), rtip, np.round(c0, 1), a.r0, a.r1))
    f.write('\n'.join(lines) + '\n')
print('%d maps: swing %d steps (%.1f A per step at the tip), breathing %d steps (alpha step %.3f)' % (a.nk * a.nj, a.nk, a.travel / (a.nk - 1), a.nj, a.amax / (a.nj - 1)))
