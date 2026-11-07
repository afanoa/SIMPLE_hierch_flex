#!/usr/bin/env python3
"""IgM ICA variants (follow-up of igm_ica.py): the arm angles are uniform (sub-Gaussian, kurtosis -1.2) but per-particle noise dilutes
them to ~-0.2, while nuisance/outlier directions are super-Gaussian (+0.5..+13), which plain logcosh FastICA picks. Variants on the
bin-centred latent whitened to its top NW PCs: FastICA logcosh / cube / exp, and SUB-GAUSSIAN PURSUIT (deflationary projected
gradient descent minimising the kurtosis of w'x, 6 random starts per component, keep the most negative). Scored per arm angle by the
Hungarian-matched |corr| on non-junk particles; the supervised ridge ceiling for reference.
usage: igm_ica2.py --flex F --bins bins.npz --embed subcls_embed.npz --gt gt --ttable truth"""
import argparse, numpy as np
from scipy.optimize import linear_sum_assignment
from sklearn.decomposition import FastICA
ap = argparse.ArgumentParser(); ap.add_argument('--flex', required=True); ap.add_argument('--bins', required=True); ap.add_argument('--embed', required=True)
ap.add_argument('--gt', required=True); ap.add_argument('--ttable', required=True); a = ap.parse_args(); rng = np.random.default_rng(0)
F = np.loadtxt(a.flex, comments='#', ndmin=2); p = F[:, 0].astype(int); Z = F[:, 5:]; N, d = Z.shape
lab = np.load(a.bins)['label']; bi = np.unique(lab[p - 1], return_inverse=True)[1]
G = np.loadtxt(a.gt, comments='#', ndmin=2); conf_of = np.zeros(G[:, 0].astype(int).max() + 1, int); conf_of[G[:, 0].astype(int)] = G[:, 1].astype(int); conf = conf_of[p]; junk = conf == 0
TT = np.loadtxt(a.ttable, comments='#', ndmin=2); tab = np.full((int(max(conf.max(), TT[:, 0].max())) + 1, TT.shape[1] - 1), np.nan); tab[TT[:, 0].astype(int)] = TT[:, 1:]
A = tab[conf]; A[junk] = np.nan; names = ['a%d_%s' % (j // 2 + 1, ('th', 'ph')[j % 2]) for j in range(A.shape[1])]; na = len(names)
E = np.load(a.embed); g_of = np.full(N, -2); pos = {q: i for i, q in enumerate(p)}; gp = E['bincentred_sgrp'][E['sid']]; g_of[np.array([pos[q] for q in E['p']])] = gp
ug, cg = np.unique(gp[gp >= 0], return_counts=True); main = ug[np.argmax(cg)]
def centre(m):
    b = bi[m]; nb = bi.max() + 1; c = np.bincount(b, minlength=nb); M = np.stack([np.bincount(b, Z[m][:, j], minlength=nb) for j in range(d)], 1) / np.maximum(c, 1)[:, None]
    Y = Z[m] - M[b]; return Y - Y.mean(0)
def matched(S, Am):
    ok = np.isfinite(Am).all(1); Sc = (S[ok] - S[ok].mean(0)) / S[ok].std(0); Ac = (Am[ok] - Am[ok].mean(0)) / Am[ok].std(0); Cm = np.abs(Sc.T @ Ac) / ok.sum()
    r, c = linear_sum_assignment(-Cm); out = np.zeros(na)
    for i, j in zip(r, c): out[j] = Cm[i, j]
    return out
def kurt(y): y = (y - y.mean()) / y.std(); return (y ** 4).mean() - 3
def subgauss(X, k, starts=6, iters=300, eta=0.05):   # X whitened (n x q); deflationary minimisation of kurtosis on the sphere
    W = []
    for c in range(k):
        best = None
        for s in range(starts):
            w = rng.standard_normal(X.shape[1])
            for _ in range(iters):
                if W: w -= np.array(W).T @ (np.array(W) @ w)
                w /= np.linalg.norm(w); y = X @ w; g = 4 * (X * (y ** 3)[:, None]).mean(0) - 12 * w; w = w - eta * g
            if W: w -= np.array(W).T @ (np.array(W) @ w)
            w /= np.linalg.norm(w); kk = kurt(X @ w)
            if best is None or kk < best[0]: best = (kk, w)
        W.append(best[1])
    return X @ np.array(W).T, [kurt(X @ w) for w in W]
for sn, m in (('main (truth-free junk removal)', g_of == main), ('oracle (truth junk removed)', ~junk)):
    Y = centre(m); Am = A[m]; w, U = np.linalg.eigh(np.cov(Y.T)); U = U[:, ::-1]; w = w[::-1]
    print('\n== %s: %d particles; matched |corr| per angle (%s) | mean' % (sn, m.sum(), ' '.join(names)))
    for nw in (10, 12, 16):
        X = (Y @ U[:, :nw]) / np.sqrt(w[:nw]); print('  whitened to the top %d PCs:' % nw)
        print('    %-26s %s | %.3f' % ('PCA', ' '.join('%.2f' % x for x in matched(X, Am)), matched(X, Am).mean()))
        for fun in ('logcosh', 'cube', 'exp'):
            S = FastICA(n_components=nw, whiten='unit-variance', fun=fun, random_state=0, max_iter=2000, tol=1e-5).fit_transform(X); mm = matched(S, Am)
            print('    %-26s %s | %.3f' % ('FastICA ' + fun, ' '.join('%.2f' % x for x in mm), mm.mean()))
        S, ks = subgauss(X, min(nw, na)); mm = matched(S, Am)
        print('    %-26s %s | %.3f   (kurtosis of the pursued directions %s)' % ('sub-Gaussian pursuit', ' '.join('%.2f' % x for x in mm), mm.mean(), ' '.join('%+.2f' % x for x in ks)))
