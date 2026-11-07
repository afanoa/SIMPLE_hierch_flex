#!/usr/bin/env python3
"""Is a rare reference population separable in a node's flex latent? Node latent (flex_pca_coordinates z1..) centred per direction bin
(spherical k-means of the projection directions, as cheapkids.py); reference = rows of chosen states of an earlier run (e.g. f692r's
closed states). Reports: cross-validated LDA AUC (reference vs rest) in the full latent; per PC and per ICA component (FastICA on the
whitened latent, truth-free): skewness, excess kurtosis, AUC and the reference's share of the component's 2 % / 5 % tails; the truth-free
pick = the IC with the largest |skew|. usage: closed_latent.py --flex F --oris DUMP --rows ROWS --ref UNITS --ids 16,15 [--perbin 1300]"""
import argparse, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('--flex', required=True); ap.add_argument('--oris', required=True); ap.add_argument('--rows', required=True)
ap.add_argument('--ref', required=True); ap.add_argument('--ids', required=True); ap.add_argument('--perbin', type=int, default=1300); a = ap.parse_args()
rng = np.random.default_rng(0); rows = np.loadtxt(a.rows).astype(int); NR = len(rows); e1 = np.zeros(NR); e2 = np.zeros(NR); i = 0
for ln in open(a.oris):
    if 'e1=' not in ln: continue
    kv = dict(t.split('=', 1) for t in ln.split() if '=' in t); e1[i] = float(kv['e1']); e2[i] = float(kv['e2']); i += 1
F = np.loadtxt(a.flex, comments='#', ndmin=2); fi = F[:, 0].astype(int); keep = rows[fi - 1] > 0; p = fi[keep]; Z = F[keep, 5:]; N = len(p)
E1, E2 = np.radians(e1[p - 1]), np.radians(e2[p - 1]); v = np.c_[np.sin(E2) * np.cos(E1), np.sin(E2) * np.sin(E1), np.cos(E2)]
nb = max(20, int(round(N / a.perbin))); C = v[rng.choice(N, nb, replace=False)].copy()
for it in range(40):
    lab = np.argmax(v @ C.T, 1)
    for k in range(nb):
        m = lab == k; C[k] = v[m].mean(0) / np.linalg.norm(v[m].mean(0)) if m.any() else v[rng.integers(N)]
lab = np.argmax(v @ C.T, 1); _, bi = np.unique(lab, return_inverse=True); nbb = bi.max() + 1; cnt = np.bincount(bi)
bm = np.stack([np.bincount(bi, Z[:, j], minlength=nbb) for j in range(Z.shape[1])], 1) / cnt[:, None]; Zc = Z - bm[bi]; Zc -= Zc.mean(0)
ref = np.loadtxt(a.ref, dtype=int); ids = [int(x) for x in a.ids.split(',')]; y = np.isin(ref[p - 1], ids)
print('node %d particles, %d-dim latent, %d direction bins; reference states %s: %d particles in the node (%.2f %%)' % (N, Z.shape[1], nbb, a.ids, y.sum(), 100 * y.mean()))
def auc(s, y):
    r = np.argsort(np.argsort(s)) + 1.; n1 = y.sum(); n0 = len(y) - n1; return (r[y].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)
# cross-validated LDA (2 folds)
f = rng.random(N) < 0.5; aucs = []
for tr, te in ((f, ~f), (~f, f)):
    mu1, mu0 = Zc[tr & y].mean(0), Zc[tr & ~y].mean(0); S = np.cov(Zc[tr].T) + 1e-6 * np.eye(Zc.shape[1]); w = np.linalg.solve(S, mu1 - mu0); aucs.append(auc(Zc[te] @ w, y[te]))
print('LDA (cross-validated) AUC reference vs rest: %.3f / %.3f' % tuple(aucs))
def report(name, s):
    s = (s - s.mean()) / s.std(); sk = (s ** 3).mean(); ku = (s ** 4).mean() - 3; A = auc(s, y); sgn = 1 if A >= 0.5 else -1; t = sgn * s
    tails = ' '.join('top %d %%: %.0f %% ref (%.1fx)' % (q, 100 * y[t >= np.percentile(t, 100 - q)].mean(), y[t >= np.percentile(t, 100 - q)].mean() / max(y.mean(), 1e-9)) for q in (2, 5))
    print('  %-5s skew %+5.2f kurt %+6.2f | AUC %.3f | %s' % (name, sk, ku, max(A, 1 - A), tails)); return sk, ku, A
w_, U = np.linalg.eigh(np.cov(Zc.T)); U = U[:, ::-1]; w_ = w_[::-1]; P = Zc @ U
for j in range(min(4, P.shape[1])): report('PC%d' % (j + 1), P[:, j])
# FastICA (tanh) on the whitened latent
K = min(Z.shape[1], 12); Wh = U[:, :K] / np.sqrt(w_[:K]); X = Zc @ Wh; W = np.linalg.qr(rng.normal(size=(K, K)))[0]
for it in range(200):
    G = np.tanh(X @ W.T); Wn = (G.T @ X) / len(X) - np.diag((1 - G ** 2).mean(0)) @ W; u_, s_, vt_ = np.linalg.svd(Wn); Wn = u_ @ vt_
    if np.max(np.abs(np.abs(np.sum(Wn * W, 1)) - 1)) < 1e-6: W = Wn; break
    W = Wn
S_ = X @ W.T; res = [(report('IC%d' % (j + 1), S_[:, j]), j) for j in range(K)]
pick = max(res, key=lambda r: abs(r[0][0])); print('truth-free pick: IC%d (largest |skew| %.2f) AUC %.3f' % (pick[1] + 1, pick[0][0], max(pick[0][2], 1 - pick[0][2])))
