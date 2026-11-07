#!/usr/bin/env python3
"""Do a node's flex states sort particles by their pose error (problem 3: a latent fitted at imperfect poses absorbs rotation modes
(w x r).grad V, which differ between states over the whole molecule)? Per state of a labelling: pose error vs truth (geodesic, SIMPLE
euler2m), and the mean signed rotation error (mean of the error rotation vectors R_ref R_true^T as axis * angle, deg) - a state whose
particles share a systematic rotation reconstructs a rotated map. Random split of the same sizes = the sampling floor.
usage: p3_poserr.py DUMP LABELS TRUE_ORIS [--gt GT --name N]"""
import re, argparse, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('dump'); ap.add_argument('labels'); ap.add_argument('true'); ap.add_argument('--gt', default=''); ap.add_argument('--name', default='set'); a = ap.parse_args()
def oris(fn, keys):
    pat = {k: re.compile(r'(?:^| )%s=([-0-9.eE+]+)' % k) for k in keys}
    return np.array([[float(pat[k].search(l).group(1)) for k in keys] for l in open(fn) if l.strip()])
X, Tr = oris(a.dump, ('e1', 'e2', 'e3', 'x', 'y')), oris(a.true, ('e1', 'e2', 'e3')); lab = np.loadtxt(a.labels, dtype=int)
def rz(t): c, s = np.cos(t), np.sin(t); z = np.zeros_like(t); o = np.ones_like(t); return np.stack([np.stack([c, s, z], -1), np.stack([-s, c, z], -1), np.stack([z, z, o], -1)], -2)
def ry(t): c, s = np.cos(t), np.sin(t); z = np.zeros_like(t); o = np.ones_like(t); return np.stack([np.stack([c, z, -s], -1), np.stack([z, o, z], -1), np.stack([s, z, c], -1)], -2)
def rot(e): e = np.deg2rad(e); return rz(e[:, 2]) @ ry(e[:, 1]) @ rz(e[:, 0])
Q = rot(X[:, :3]).transpose(0, 2, 1) @ rot(Tr[:, :3])   # error rotation in the volume frame: R_ref^T R_true
ang = np.arccos(np.clip((np.trace(Q, axis1=1, axis2=2) - 1) / 2, -1, 1)); s_ = np.sin(ang); ax = np.stack([Q[:, 2, 1] - Q[:, 1, 2], Q[:, 0, 2] - Q[:, 2, 0], Q[:, 1, 0] - Q[:, 0, 1]], 1) / np.maximum(2 * s_, 1e-9)[:, None]
rv = np.rad2deg(ax * ang[:, None]); ad = np.rad2deg(ang); ok = lab > 0
if a.gt:
    G = np.loadtxt(a.gt, comments='#'); conf = np.full(len(lab), -1); gi = G[:, 0].astype(int) - 1; v = (gi >= 0) & (gi < len(lab)); conf[gi[v]] = G[v, 1].astype(int); ok &= conf > 0
rng = np.random.default_rng(1); rl = lab.copy(); rl[ok] = rng.permutation(lab[ok])
def rows(l):
    return [(s, (l[ok] == s).sum(), np.median(ad[ok][l[ok] == s]), np.linalg.norm(rv[ok][l[ok] == s].mean(0)), np.hypot(*X[ok][l[ok] == s][:, 3:5].mean(0))) for s in np.unique(l[ok])]
for tag, l in (('states', lab), ('random', rl)):
    r = rows(l); print('%s %s: median pose error %s deg | |mean rotation error| %s deg | |mean shift| %s px' % (a.name, tag, ' '.join('%d:%.2f' % (x[0], x[2]) for x in r),
          ' '.join('%.2f' % x[3] for x in r), ' '.join('%.2f' % x[4] for x in r)))
