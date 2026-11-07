#!/usr/bin/env python3
"""Truth-free per-arm coordinates on the IgM-like root: the node's flex latent, centred per direction bin, whitened, FastICA. The arm
angles are independent and non-Gaussian, so independent components of the latent should each follow one arm angle. Particle sets:
'all' (no junk removal), 'main' (truth-free: members of the sub-class embedding's largest group, i.e. the junk groups dropped) and
'oracle' (truth junk removed). Scored on the non-junk particles of the set: per angle the |corr| of its Hungarian-matched component
(one component per angle), per component its best angle and the margin over its second-best; the same for PCA components; seed
stability (|corr| of matched components between FastICA seeds); the supervised ceiling = held-out ridge R2 from the whole latent.
usage: igm_ica.py --flex F --bins bins.npz --embed subcls_embed.npz --gt gt --ttable truth --out DIR [--seeds 3]"""
import argparse, os, numpy as np
from scipy.optimize import linear_sum_assignment
from sklearn.decomposition import FastICA
ap = argparse.ArgumentParser(); ap.add_argument('--flex', required=True); ap.add_argument('--bins', required=True); ap.add_argument('--embed', required=True)
ap.add_argument('--gt', required=True); ap.add_argument('--ttable', required=True); ap.add_argument('--out', required=True); ap.add_argument('--seeds', type=int, default=3)
a = ap.parse_args(); os.makedirs(a.out, exist_ok=True); rng = np.random.default_rng(0)
F = np.loadtxt(a.flex, comments='#', ndmin=2); p = F[:, 0].astype(int); Z = F[:, 5:]; N, d = Z.shape
lab = np.load(a.bins)['label']; bi = np.unique(lab[p - 1], return_inverse=True)[1]
G = np.loadtxt(a.gt, comments='#', ndmin=2); conf_of = np.zeros(G[:, 0].astype(int).max() + 1, int); conf_of[G[:, 0].astype(int)] = G[:, 1].astype(int); conf = conf_of[p]; junk = conf == 0
TT = np.loadtxt(a.ttable, comments='#', ndmin=2); tab = np.full((int(max(conf.max(), TT[:, 0].max())) + 1, TT.shape[1] - 1), np.nan); tab[TT[:, 0].astype(int)] = TT[:, 1:]
A = tab[conf]; A[junk] = np.nan; names = ['a%d_%s' % (j // 2 + 1, ('th', 'ph')[j % 2]) for j in range(A.shape[1])]; na = len(names)
E = np.load(a.embed); ep, esid, sg = E['p'], E['sid'], E['bincentred_sgrp']; g_of = np.full(N, -2); pos = {q: i for i, q in enumerate(p)}
gp = sg[esid]; idx = np.array([pos[q] for q in ep]); g_of[idx] = gp
ug, cg = np.unique(gp[gp >= 0], return_counts=True); main = ug[np.argmax(cg)]
sets = {'all': np.ones(N, bool), 'main (truth-free junk removal)': g_of == main, 'oracle (truth junk removed)': ~junk}
def centre(X, m):
    b = bi[m]; nb = bi.max() + 1; c = np.bincount(b, minlength=nb); M = np.stack([np.bincount(b, X[m][:, j], minlength=nb) for j in range(X.shape[1])], 1) / np.maximum(c, 1)[:, None]
    Y = X[m] - M[b]; return Y - Y.mean(0)
def cmat(S, Am):  # |corr| components x angles on rows with all angles defined
    ok = np.isfinite(Am).all(1); Sc = (S[ok] - S[ok].mean(0)) / S[ok].std(0); Ac = (Am[ok] - Am[ok].mean(0)) / Am[ok].std(0); return np.abs(Sc.T @ Ac) / ok.sum()
def match(Cm):
    r, c = linear_sum_assignment(-Cm); out = np.zeros(na); comp = np.zeros(na, int)
    for i, j in zip(r, c): out[j] = Cm[i, j]; comp[j] = i
    return out, comp
def ridge_r2(X, y, nfold=5):
    ok = np.isfinite(y); X = X[ok]; y = y[ok]; f = rng.permutation(len(y)) % nfold; pred = np.zeros(len(y))
    for k in range(nfold):
        tr, te = f != k, f == k; mu = X[tr].mean(0); Xa = X[tr] - mu; W = np.linalg.solve(Xa.T @ Xa + 1e-3 * np.trace(Xa.T @ Xa) / X.shape[1] * np.eye(X.shape[1]), Xa.T @ (y[tr] - y[tr].mean()))
        pred[te] = (X[te] - mu) @ W + y[tr].mean()
    return 1 - ((y - pred) ** 2).sum() / ((y - y.mean()) ** 2).sum()
print('IgM root: %d particles, %d-dim latent, %d bins; junk %d; sub-class groups %d, main group %d (%d particles, junk share %.2f)'
      % (N, d, bi.max() + 1, junk.sum(), len(ug), main, (g_of == main).sum(), junk[g_of == main].mean()))
for sn, m in sets.items():
    Y = centre(Z, m); Am = A[m]; print('\n== %s: %d particles, junk %d (%.1f %%)' % (sn, m.sum(), junk[m].sum(), 100 * junk[m].mean()))
    print('   supervised ceiling, held-out ridge R2 per angle (|corr| = sqrt): ' + ' '.join('%s %.2f' % (nm, ridge_r2(Y, Am[:, j])) for j, nm in enumerate(names)))
    w, U = np.linalg.eigh(np.cov(Y.T)); U = U[:, ::-1]; P = Y @ U; mp, _ = match(cmat(P, Am))
    print('   PCA components, matched |corr| per angle:        ' + ' '.join('%s %.2f' % (nm, x) for nm, x in zip(names, mp)) + ' | mean %.3f' % mp.mean())
    S0 = None; res = []
    for s in range(a.seeds):
        ica = FastICA(n_components=d, whiten='unit-variance', random_state=s, max_iter=2000, tol=1e-5); S = ica.fit_transform(Y)
        Cm = cmat(S, Am); mi, comp = match(Cm); res.append(mi)
        if s == 0:
            S0, C0, comp0 = S, Cm, comp; k = np.abs(((S - S.mean(0)) / S.std(0))).max(0)
            np.save(os.path.join(a.out, 'ica_%s_seed0.npy' % sn.split()[0]), S.astype(np.float32))
        else:
            X0 = (S0 - S0.mean(0)) / S0.std(0); X1 = (S - S.mean(0)) / S.std(0); cc = np.abs(X0.T @ X1) / len(S); r_, c_ = linear_sum_assignment(-cc)
            print('   seed %d vs seed 0: matched component |corr| median %.3f (min %.3f)' % (s, np.median(cc[r_, c_]), cc[r_, c_].min()))
    res = np.array(res)
    print('   ICA components, matched |corr| per angle (seed mean +- sd): ' + ' '.join('%s %.2f+-%.2f' % (nm, res[:, j].mean(), res[:, j].std()) for j, nm in enumerate(names)) + ' | mean %.3f' % res.mean())
    kurt = lambda x: ((x - x.mean()) ** 4).mean() / x.var() ** 2 - 3
    lines = []
    for j, nm in enumerate(names):
        i = comp0[j]; row = C0[i].copy(); best = row[j]; row[j] = 0; jj = int(np.argmax(row))
        lines.append('%s<-IC%d %.2f (next %s %.2f, kurtosis %+.2f)' % (nm, i + 1, best, names[jj], row[jj], kurt(S0[:, i])))
    print('   seed 0 assignment: ' + '; '.join(lines))
    un = [i for i in range(d) if i not in comp0]; print('   unassigned components (kurtosis, junk |corr|): ' + ' '.join('IC%d %+.2f %.2f' % (i + 1, kurt(S0[:, i]), abs(np.corrcoef(S0[:, i], junk[m].astype(float))[0, 1]) if junk[m].any() else 0) for i in un))
