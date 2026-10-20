#!/usr/bin/env python3
"""Does the view-centred flex_pca latent of the direction-bin sub-class averages carry the modes? Per sub-class truth =
member mean of the arm angles (arms) or site occupancies (tethered); readouts: (1) cross-validated ridge R2 of each truth
column from the sub-class coordinates z1..zq (full fit), (2) even/odd reproducibility of the coordinates (per component
correlation between the even-average and odd-average fits, sign-aligned), (3) k-means placement (K states on the
whitened leading components) -> particle labels for label_truth_readout.
usage: avgflex_readout.py AVGFLEX_DIR CLSEXP_DIR GT_PARTICLES TRUTH MODE --k 6 --q 6"""
import sys, os, argparse, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('d'); ap.add_argument('clsexp'); ap.add_argument('gt'); ap.add_argument('truth'); ap.add_argument('mode'); ap.add_argument('--k', type=int, default=6); ap.add_argument('--q', type=int, default=6); a = ap.parse_args()
conf = {int(w[0]): int(w[1]) for w in (l.split() for l in open(a.gt)) if w and not w[0].startswith('#')}
T = {}
for w in (l.split() for l in open(a.truth)):
    if not w or w[0].startswith('#'): continue
    v = np.array([float(x) for x in w[1:]]); T[int(w[0])] = v[0::4] if a.mode.startswith('tethered') else v
cols = ['site%d' % (i + 1) for i in range(len(next(iter(T.values()))))] if a.mode.startswith('tethered') else sum([['arm%d_theta' % (i + 1), 'arm%d_phi' % (i + 1)] for i in range(len(next(iter(T.values()))) // 2)], [])
L = np.loadtxt(a.clsexp + '/cls_expansion_labels.txt', comments='#').astype(int); sub = L[:, 2]; pidx = L[:, 0]; N = sub.max()
truth_p = np.array([T[conf[p]] if conf.get(p) in T else np.full(len(cols), np.nan) for p in pidx])
ts = np.full((N + 1, len(cols)), np.nan); cnt = np.zeros(N + 1, int)
for g in range(1, N + 1):
    m = sub == g; cnt[g] = m.sum()
    if cnt[g]: ts[g] = np.nanmean(truth_p[m], 0)
def load(h):
    f = a.d + '/run_%s/flex_pca_coordinates.txt' % h
    if not os.path.isfile(f): return None
    X = np.loadtxt(f, comments='#'); return X[:, 0].astype(int), X[:, 5:5 + a.q]
Z = {h: load(h) for h in ('full', 'even', 'odd')}
def ridge_cv(X, y, lam=1e-2, nf=5):
    ok = np.isfinite(y); X, y = X[ok], y[ok]; n = len(y); idx = np.arange(n) % nf; pred = np.zeros(n)
    Xs = (X - X.mean(0)) / (X.std(0) + 1e-12); Xa = np.c_[Xs, np.ones(n)]
    for f in range(nf):
        tr, te = idx != f, idx == f; w = np.linalg.solve(Xa[tr].T @ Xa[tr] + lam * np.eye(Xa.shape[1]), Xa[tr].T @ y[tr]); pred[te] = Xa[te] @ w
    return 1 - np.sum((y - pred) ** 2) / np.sum((y - y.mean()) ** 2)
if Z['full'] is not None:
    gid, X = Z['full']; y = ts[gid]; print('sub-class latent (full fit, q=%d): cross-validated R2 per truth column: %s' % (a.q, ' '.join('%s %.2f' % (c, ridge_cv(X, y[:, j])) for j, c in enumerate(cols))))
if Z['even'] is not None and Z['odd'] is not None:
    (ge, Xe), (go, Xo) = Z['even'], Z['odd']; assert (ge == go).all()
    r = [abs(np.corrcoef(Xe[:, j], Xo[:, j])[0, 1]) for j in range(Xe.shape[1])]; print('even/odd coordinate reproducibility per component: %s' % ' '.join('%.2f' % v for v in r))
    # the odd truth predicted from a model fitted on the even coordinates (cross-half R2)
    def xhalf(j):
        ok = np.isfinite(ts[ge][:, j]) & np.isfinite(ts[go][:, j]); mu, sd = Xe[ok].mean(0), Xe[ok].std(0) + 1e-12
        A = np.c_[(Xe[ok] - mu) / sd, np.ones(ok.sum())]; Bm = np.c_[(Xo[ok] - mu) / sd, np.ones(ok.sum())]
        w = np.linalg.solve(A.T @ A + 1e-2 * np.eye(A.shape[1]), A.T @ ts[ge][ok, j]); yo = ts[go][ok, j]; return 1 - np.sum((yo - Bm @ w) ** 2) / np.sum((yo - yo.mean()) ** 2)
    print('cross-half: ridge fitted on the even coordinates, R2 on the odd ones: %s' % ' '.join('%s %.2f' % (c, xhalf(j)) for j, c in enumerate(cols)))
# k-means placement on the whitened leading components of the full fit -> particle labels
if Z['full'] is not None:
    gid, X = Z['full']; Xs = (X - X.mean(0)) / (X.std(0) + 1e-12); rng = np.random.default_rng(0); C = Xs[rng.choice(len(Xs), a.k, replace=False)]
    for it in range(50):
        lab = np.argmin(((Xs[:, None, :] - C[None]) ** 2).sum(2), 1)
        for k in range(a.k):
            if (lab == k).any(): C[k] = Xs[lab == k].mean(0)
    sl = np.zeros(N + 1, int); sl[gid] = lab + 1; pl = np.zeros(pidx.max(), int); pl[pidx - 1] = sl[sub]; np.savetxt(a.d + '/particle_states_kmeans.txt', pl, fmt='%d')
    print('k-means (K=%d) on the whitened latent: sub-classes per state %s -> %s/particle_states_kmeans.txt' % (a.k, np.bincount(lab, minlength=a.k).tolist(), a.d))
