#!/usr/bin/env python3
"""Is a node's flex latent stable across runs even when its K-state assignment is not? For pairs of flex_pca_coordinates.txt files of the
same node (rows matched by particle index): canonical correlations of the 10-D latents; NMI across runs of (a) the runs' own K-state labels,
(b) 8 quantile bins of each run's first canonical variate, (c) k-means (k 8, fixed init) on each run's first two canonical variates; and the
old f692r closed share per quantile bin of the first canonical variate (is it the arm-closing direction?).  usage: latent_cca.py NAME=FILE ..."""
import sys, numpy as np
K_ = '/mnt/beegfs/elmlund/afan/rnasep_klem'
old = np.loadtxt(K_ + '/f692r/coordtree_fam1/units.txt', dtype=int); oc = (old == 15) | (old == 16)
def nmi(a, b):
    ua, ia = np.unique(a, return_inverse=True); ub, ib = np.unique(b, return_inverse=True)
    C = np.zeros((len(ua), len(ub))); np.add.at(C, (ia, ib), 1); P = C / C.sum(); pa, pb = P.sum(1), P.sum(0); nz = P > 0
    I = (P[nz] * np.log(P[nz] / np.outer(pa, pb)[nz])).sum(); H = lambda p: -(p[p > 0] * np.log(p[p > 0])).sum(); return I / np.sqrt(H(pa) * H(pb))
def load(f):
    F = np.loadtxt(f, comments='#', ndmin=2); return F[:, 0].astype(int), F[:, 2].astype(int), F[:, 5:]
def qbin(x, k=8): return np.searchsorted(np.quantile(x, np.linspace(0, 1, k + 1)[1:-1]), x)
def kmeans(X, k=8, it=50):
    C = X[np.argsort(X[:, 0])[np.linspace(0, len(X) - 1, k).astype(int)]].copy()   # deterministic init: spread along the first variate
    for _ in range(it):
        lab = np.argmin(((X[:, None, :] - C[None]) ** 2).sum(2), 1); C = np.array([X[lab == j].mean(0) if (lab == j).any() else C[j] for j in range(k)])
    return lab
runs = dict(a.split('=', 1) for a in sys.argv[1:]); D = {n: load(f) for n, f in runs.items()}
names = list(runs)
for i in range(len(names)):
    for j in range(i + 1, len(names)):
        a, b = names[i], names[j]; pa, la, Za = D[a]; pb, lb, Zb = D[b]
        common, ia, ib = np.intersect1d(pa, pb, return_indices=True); X = Za[ia] - Za[ia].mean(0); Y = Zb[ib] - Zb[ib].mean(0)
        Qx, _ = np.linalg.qr(X); Qy, _ = np.linalg.qr(Y); U, s, Vt = np.linalg.svd(Qx.T @ Qy)
        u = Qx @ U; v = Qy @ Vt.T; sg = np.sign((u[:, 0] * v[:, 0]).sum()); u1, v1 = u[:, 0], sg * v[:, 0]
        r = common - 1; base = oc[r].mean(); qa = qbin(u1)
        cs = ' '.join('%.2f' % (oc[r][qa == q].mean() / base) for q in range(8))
        print('%s vs %s (%d particles): canonical correlations %s' % (a, b, len(common), ' '.join('%.2f' % x for x in s[:6])))
        print('   NMI across runs: own K-state labels %.2f | 8 quantile bins of the first canonical variate %.2f | k-means(8) on the first two %.2f'
              % (nmi(la[ia], lb[ib]), nmi(qa, qbin(v1)), nmi(kmeans(np.c_[u1, u[:, 1]] / u1.std()), kmeans(np.c_[v1, sg * v[:, 1]] / v1.std()))))
        print('   old-closed enrichment per quantile bin of the first canonical variate (run %s): %s' % (a, cs))
