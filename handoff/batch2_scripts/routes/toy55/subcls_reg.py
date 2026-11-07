#!/usr/bin/env python3
"""2D -> 3D by per-view regression onto the node's 3D latent (prototype, truth readouts optional). Inside every direction bin the
members' cls_expansion 2D latent (SIMPLE_FLEXCLS_DUMP: posterior means z2, view-specific, arbitrary gauge per bin) is mapped into the
node's 3D latent frame (flex_pca coordinates z3, common to all views) by a ridge regression z3 ~ [1, z2] fitted on the OTHER half of
the bin's members (cross-fitting: a particle's prediction never sees its own z3 noise). Each bin is registered to the one common 3D
reference; no bin is synchronised with another. Outputs per particle: zhat = the 2D-predicted 3D coordinate. Readouts: per variable
(truth) held-out ridge R2 from z3, zhat, [z3 zhat], each also pooled over the cls_expansion sub-classes; truth-free: split-half
reliability over sub-classes (as subcls_embed2.py) of z3 vs zhat pooled.
usage: subcls_reg.py --dumps DIR(cls_expansion_latent_class*.txt) --labels cls_expansion_labels.txt --flex flex_pca_coordinates.txt --rows ROWS
       --out DIR [--gt gt (--codes swing:200:230 | --ttable T)] [--q 16] [--alpha 0.05]"""
import argparse, os, glob, numpy as np
from scipy.linalg import eigh
ap = argparse.ArgumentParser(); ap.add_argument('--dumps', required=True); ap.add_argument('--labels', required=True); ap.add_argument('--flex', required=True)
ap.add_argument('--rows', required=True); ap.add_argument('--out', required=True); ap.add_argument('--gt', default=''); ap.add_argument('--codes', default='')
ap.add_argument('--ttable', default=''); ap.add_argument('--q', type=int, default=16); ap.add_argument('--alpha', type=float, default=0.05); ap.add_argument('--minh', type=int, default=5)
ap.add_argument('--nnull', type=int, default=20); ap.add_argument('--kwin', type=int, default=8)
a = ap.parse_args(); os.makedirs(a.out, exist_ok=True); rng = np.random.default_rng(0)
rows = np.loadtxt(a.rows).astype(int); NR = len(rows)
L = np.loadtxt(a.labels, comments='#', ndmin=2).astype(int); sub_of = np.zeros(NR + 1, int); bin_of = np.zeros(NR + 1, int); sub_of[L[:, 0]] = L[:, 2]; bin_of[L[:, 0]] = L[:, 1]
F = np.loadtxt(a.flex, comments='#', ndmin=2); fi = F[:, 0].astype(int); z3_of = np.full((NR + 1, F.shape[1] - 5), np.nan); z3_of[fi] = F[:, 5:]; eo_of = np.full(NR + 1, -1); eo_of[fi] = F[:, 1].astype(int)
fl = sorted(glob.glob(os.path.join(a.dumps, 'cls_expansion_latent_class*.txt'))); assert fl, 'no dumps in ' + a.dumps
hdr = open(fl[0]).readline(); nc = int(hdr.split('z(1:')[1].split(')')[0]); q = min(a.q, nc)
P, Z2, B = [], [], []
for f in fl:
    X = np.loadtxt(f, comments='#', ndmin=2)
    if len(X) == 0: continue
    P.append(X[:, 0].astype(int)); Z2.append(X[:, 2:2 + q]); B.append(np.full(len(X), int(f.split('class')[-1].split('.')[0])))
P = np.concatenate(P); Z2 = np.concatenate(Z2); B = np.concatenate(B)
keep = (rows[P - 1] > 0) & np.isfinite(z3_of[P]).all(1) & (eo_of[P] >= 0) & (sub_of[P] > 0)
P, Z2, B = P[keep], Z2[keep], B[keep]; Z3 = z3_of[P]; eo = eo_of[P]; N, d3 = Z3.shape
us, sid = np.unique(sub_of[P], return_inverse=True); S = len(us)
print('%d members with a 2D latent (%d of %d components used) and a %d-dim 3D latent, %d bins, %d sub-classes' % (N, q, nc, d3, len(np.unique(B)), S))
# cross-fitted per-bin ridge z3 ~ [1, z2]
Zh = np.full_like(Z3, np.nan); r2in = []
for b in np.unique(B):
    m = np.where(B == b)[0]
    for h in (0, 1):
        tr = m[eo[m] != h]; te = m[eo[m] == h]
        if len(tr) < 3 * q or len(te) == 0: Zh[te] = Z3[tr].mean(0) if len(tr) else np.nan; continue
        mx = Z2[tr].mean(0); my = Z3[tr].mean(0); Xc = Z2[tr] - mx; lam = a.alpha * np.trace(Xc.T @ Xc) / q
        W = np.linalg.solve(Xc.T @ Xc + lam * np.eye(q), Xc.T @ (Z3[tr] - my)); Zh[te] = (Z2[te] - mx) @ W + my
        pr = Xc @ W + my; r2in.append(1 - ((Z3[tr] - pr) ** 2).sum() / ((Z3[tr] - my) ** 2).sum())
ok_ = np.isfinite(Zh).all(1); print('per-bin regression: in-sample R2 of z3 on z2 median %.3f (p10 %.3f, p90 %.3f); %d particles predicted' % (np.median(r2in), *np.percentile(r2in, [10, 90]), ok_.sum()))
Zh[~ok_] = Z3[~ok_]
oos = 1 - ((Z3 - Zh) ** 2).sum() / ((Z3 - Z3.mean(0)) ** 2).sum(); print('out-of-half R2 of z3 from the other half\'s regression: %.3f (the share of the 3D latent the 2D latent predicts)' % oos)
np.savez(os.path.join(a.out, 'subcls_reg.npz'), p=P, eo=eo, bin=B, sid=sid, Z3=Z3, Zh=Zh)

def gmeans(X, g, G, w=None):
    w = np.ones(len(g)) if w is None else w; c = np.bincount(g, w, minlength=G)
    return np.stack([np.bincount(g, w * X[:, j], minlength=G) for j in range(X.shape[1])], 1) / np.maximum(c, 1e-12)[:, None], c
def cv_ridge(X, y, groups, nfold=5, alpha=1e-3):
    ug = np.unique(groups); f = rng.permutation(len(ug)) % nfold; fold = f[np.searchsorted(ug, groups)]; pred = np.zeros(len(y))
    for k in range(nfold):
        tr, te = fold != k, fold == k; mu = X[tr].mean(0); A = X[tr] - mu; ym = y[tr].mean(); lamI = alpha * np.trace(A.T @ A) / A.shape[1] * np.eye(A.shape[1])
        pred[te] = (X[te] - mu) @ np.linalg.solve(A.T @ A + lamI, A.T @ (y[tr] - ym)) + ym
    return 1 - ((y - pred) ** 2).sum() / ((y - y.mean()) ** 2).sum()
def reliab(Me, Mo, w):
    W = w / w.sum(); mu = (W[:, None] * (Me + Mo) / 2).sum(0); A_ = Me - mu; B_ = Mo - mu
    Cx = (A_ * W[:, None]).T @ B_; Cx = (Cx + Cx.T) / 2; Cw = ((A_ * W[:, None]).T @ A_ + (B_ * W[:, None]).T @ B_) / 2
    lam, V = eigh(Cx, Cw + 1e-9 * np.trace(Cw) / len(Cw) * np.eye(len(Cw))); o = np.argsort(lam)[::-1]; return lam[o], V[:, o]
def binshuffle(g):
    out = g.copy()
    for b in np.unique(B):
        m = np.where(B == b)[0]; out[m] = g[m][rng.permutation(len(m))]
    return out
# bin-centring (view offsets) for both
def centre(X): bm, _ = gmeans(X, np.searchsorted(np.unique(B), B), len(np.unique(B))); return X - bm[np.searchsorted(np.unique(B), B)]
feats = {'z3': centre(Z3), 'zhat': centre(Zh), 'z3+zhat': np.c_[centre(Z3), centre(Zh)]}
for nm, X in feats.items():
    Me, ne = gmeans(X, sid, S, (eo == 0).astype(float)); Mo, no = gmeans(X, sid, S, (eo == 1).astype(float)); ok = (ne >= a.minh) & (no >= a.minh); w = ne * no / np.maximum(ne + no, 1)
    lam, V = reliab(Me[ok], Mo[ok], w[ok]); nul = []
    for t in range(a.nnull):
        g2 = binshuffle(sid); Me2, ne2 = gmeans(X, g2, S, (eo == 0).astype(float)); Mo2, no2 = gmeans(X, g2, S, (eo == 1).astype(float)); ok2 = (ne2 >= a.minh) & (no2 >= a.minh)
        nul.append(reliab(Me2[ok2], Mo2[ok2], (ne2 * no2 / np.maximum(ne2 + no2, 1))[ok2])[0][0])
    print('[%s] sub-class split-half reliability %s | null 95th pct %.3f -> %d reproducible dims' % (nm, ' '.join('%.3f' % l for l in lam[:6]), np.percentile(nul, 95), (lam > np.percentile(nul, 95)).sum()))
if a.gt:
    G = np.loadtxt(a.gt, comments='#', ndmin=2); conf_of = np.zeros(NR + 1, int); conf_of[G[:, 0].astype(int)] = G[:, 1].astype(int); conf = conf_of[P]; TV = {}
    for c in [c for c in a.codes.split(',') if c]:
        nm, lo, hi = c.split(':'); TV[nm] = np.where((conf >= int(lo)) & (conf <= int(hi)), conf - int(lo), np.nan)
    if a.ttable:
        TT = np.loadtxt(a.ttable, comments='#', ndmin=2); nv = TT.shape[1] - 1; tab = np.full((max(conf.max(), int(TT[:, 0].max())) + 1, nv), np.nan); tab[TT[:, 0].astype(int)] = TT[:, 1:]
        for j in range(nv): v = tab[conf, j].copy(); v[conf == 0] = np.nan; TV['a%d_%s' % (j // 2 + 1, ('th', 'ph')[j % 2])] = v
    from scipy.cluster.vq import kmeans2
    # k-means windows inside every bin on the 2D latent (standardised per bin): alternative sub-classes to the divisive leaves
    win = np.zeros(N, int); nxt = 0
    for b in np.unique(B):
        m = np.where(B == b)[0]; X = Z2[m] - Z2[m].mean(0); X /= X.std(0) + 1e-12
        _, lab = kmeans2(X, min(a.kwin, len(m)), minit='++', seed=int(b)); _, lab = np.unique(lab, return_inverse=True); win[m] = lab + nxt; nxt += lab.max() + 1
    _, wid = np.unique(win, return_inverse=True); W_ = wid.max() + 1
    for tn, v in TV.items():
        dv = np.isfinite(v); r2_2, r2_3 = [], []
        for b in np.unique(B):
            m = np.where((B == b) & dv)[0]
            if len(m) < 100: continue
            r2_2.append(cv_ridge(Z2[m], v[m], np.arange(len(m)))); r2_3.append(cv_ridge(Z3[m], v[m], np.arange(len(m))))
        print('TRUTH %-6s within-bin held-out ridge R2 (median over %d bins): from the 2D latent (%d comps) %.3f | from the 3D latent %.3f' % (tn, len(r2_2), q, np.median(r2_2), np.median(r2_3)))
        def e2(lab):
            u, inv = np.unique(lab[dv], return_inverse=True); mm = np.bincount(inv, v[dv]) / np.bincount(inv); return 1 - ((v[dv] - mm[inv]) ** 2).sum() / ((v[dv] - v[dv].mean()) ** 2).sum()
        Mw, _ = gmeans(feats['z3'], wid, W_); Ms, _ = gmeans(feats['z3'], sid, S)
        print('TRUTH %-6s partition eta2: divisive sub-classes %.3f | k-means windows on the 2D latent (k %d per bin, %d windows) %.3f; pooled 3D latent R2 over windows %.3f (over sub-classes %.3f)'
              % (tn, e2(sid), a.kwin, W_, e2(wid), cv_ridge(Mw[wid[dv]], v[dv], wid[dv]), cv_ridge(Ms[sid[dv]], v[dv], sid[dv])))
    np.savetxt(os.path.join(a.out, 'kwin_labels.txt'), np.c_[P, B, wid + 1], fmt='%d', header='particle bin window')
    for tn, v in TV.items():
        dv = np.isfinite(v); line = 'TRUTH %-6s held-out ridge R2 (per particle | pooled over sub-classes):' % tn
        for nm, X in feats.items():
            Ma, _ = gmeans(X, sid, S); line += '  %s %.3f | %.3f' % (nm, cv_ridge(X[dv], v[dv], sid[dv]), cv_ridge(Ma[sid[dv]], v[dv], sid[dv]))
        print(line)
