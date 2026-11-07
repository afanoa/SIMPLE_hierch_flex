#!/usr/bin/env python3
"""Validity check of the sub-class split-half reliability on a HOMOGENEOUS set (truth needed): keep only the particles of one
conformer (drop junk and everything else by truth), keep their cls_expansion sub-classes, pool the (raw / bin-centred) 3D latent per
sub-class over even / odd members and compare the top split-half reliabilities with the within-bin member-shuffle null. On one
conformer every reproducible sub-class difference is an artefact (shared image noise between the 2D assignment and the 3D latent,
view or CTF nuisance). Also the same for the 2D latent itself if a reg dir (subcls_reg.npz) is given.
usage: subcls_nullcheck.py EMBED_DIR --conf 49 [--reg REG_DIR] [--minh 5] [--nnull 20]"""
import argparse, os, numpy as np
from scipy.linalg import eigh
ap = argparse.ArgumentParser(); ap.add_argument('d'); ap.add_argument('--conf', type=int, required=True); ap.add_argument('--reg', default='')
ap.add_argument('--minh', type=int, default=5); ap.add_argument('--nnull', type=int, default=20); a = ap.parse_args(); rng = np.random.default_rng(0)
E = np.load(os.path.join(a.d, 'subcls_embed.npz')); keep = E['conf'] == a.conf
sid0, eo, Z, binn = E['sid'][keep], E['eo'][keep], E['Z'][keep], E['sbin'][E['sid']][keep]
_, sid = np.unique(sid0, return_inverse=True); S = sid.max() + 1
print('%d particles of conformer %d (of %d node members), %d sub-classes' % (keep.sum(), a.conf, len(keep), S))
def gm(X, g, G, w=None):
    w = np.ones(len(g)) if w is None else w; c = np.bincount(g, w, minlength=G)
    return np.stack([np.bincount(g, w * X[:, j], minlength=G) for j in range(X.shape[1])], 1) / np.maximum(c, 1e-12)[:, None], c
def reliab(X, g):
    Me, ne = gm(X, g, S, (eo == 0).astype(float)); Mo, no = gm(X, g, S, (eo == 1).astype(float)); ok = (ne >= a.minh) & (no >= a.minh); w = (ne * no / np.maximum(ne + no, 1))[ok]
    W = w / w.sum(); A = Me[ok] - (W[:, None] * Me[ok]).sum(0); B = Mo[ok] - (W[:, None] * Mo[ok]).sum(0)
    Cx = (A * W[:, None]).T @ B; Cx = (Cx + Cx.T) / 2; Cw = ((A * W[:, None]).T @ A + (B * W[:, None]).T @ B) / 2
    lam = eigh(Cx, Cw + 1e-9 * np.trace(Cw) / len(Cw) * np.eye(len(Cw)))[0]; return np.sort(lam)[::-1], ok.sum()
def shuffle(g):
    out = g.copy()
    for b in np.unique(binn):
        m = np.where(binn == b)[0]; out[m] = g[m][rng.permutation(len(m))]
    return out
ub, bi = np.unique(binn, return_inverse=True)
feats = {'3D latent raw': Z, '3D latent bin-centred': Z - gm(Z, bi, len(ub))[0][bi]}
if a.reg:
    R = np.load(os.path.join(a.reg, 'subcls_reg.npz')); pos = {p: i for i, p in enumerate(R['p'])}; idx = np.array([pos.get(p, -1) for p in E['p'][keep]])
    if (idx >= 0).all(): feats['3D latent predicted from the 2D latent (cross-fitted)'] = R['Zh'][idx] - gm(R['Zh'][idx], bi, len(ub))[0][bi]
for nm, X in feats.items():
    lam, n = reliab(X, sid); nul = np.array([reliab(X, shuffle(sid))[0][:3] for _ in range(a.nnull)])
    print('%-55s %d sub-classes: top reliabilities %s | null (within-bin shuffle) top %.3f (95th pct %.3f) -> %d above'
          % (nm, n, ' '.join('%.3f' % l for l in lam[:5]), nul[:, 0].max(), np.percentile(nul[:, 0], 95), (lam > np.percentile(nul[:, 0], 95)).sum()))
