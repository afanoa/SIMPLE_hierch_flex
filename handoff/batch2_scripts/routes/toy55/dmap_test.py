#!/usr/bin/env python3
"""usage: dmap_test.py NAME:KIND:TRUTH:FLEX_EMB:FFLEX_EMB ..   (KIND igm | cahra | cdisc)
Phase A, discrete vs continuous at a node by graph theory on its latent; root nodes with truth.
- Graph: per latent (flex, fflex) and k (15, 30), a symmetric kNN graph on the per-component standardised MAP z of all particles, with
  self-tuning Gaussian weights exp(-d^2 / (s_i s_j)), s_i = distance to the k-th neighbour. Diffusion map (random-walk normalised):
  the 10 leading nontrivial eigenvectors psi_j.
- Coordinates: each psi_j's best abs(Spearman) with a truth coordinate. Angles use genuine particles; binary truths use all.
- Discreteness: phi_j = the minimum sweep-cut conductance of psi_j over thresholds whose smaller side holds >= 2 % of the volume;
  the null phi_null_j is the same on a Gaussian surrogate with the latent's covariance and n; the ratio phi_j / phi_null_j.
  A second ratio, phi_j / min over the null's 10 eigenvectors, is the conservative one.
- Profile: particle occupancy along psi_j's path. The centroids of 20 equal-count bins in psi_j order trace the path; its first
  principal axis u is the motion direction; each particle's projection x.u gives the occupancy profile (30-bin histogram over the 1-99 %
  range). Modes are peaks >= 25 % of the maximum, after 3-bin smoothing; the dip is the minimum between the two highest over the lower
  one. Against a continuous truth t: the KS distance between the projection's and t's distributions, both scaled to [0, 1].
- IgM also gets 3 subsamples of genuine particles with an imposed occupancy along a4_th (two-peaked +-1 sd, skewed, uniform control),
  each rerun on fflex at k = 15."""
import sys, numpy as np, scipy.sparse as sp
from scipy.spatial import cKDTree
from scipy.sparse.linalg import eigsh
from scipy.stats import spearmanr
SPECS = [a.split(':') for a in sys.argv[1:]]; RNG = np.random.default_rng(11)
def read_emb(f):
    b = open(f, 'rb').read(); ver, n, nc, box = (int(v) for v in np.frombuffer(b[8:24], '<i4'))
    assert b[:8] == b'SIMPLFXC' and ver == 5, 'flex_pca embedding version 5 expected'
    o = 28; pinds = np.frombuffer(b[o:o + 4 * n], '<i4').astype(int); o += 4 * n
    return pinds - 1, np.frombuffer(b[o:o + 8 * n * nc], '<f8').reshape(nc, n).T.copy()
def truth(kind, tf, NM):
    T = {}
    if kind == 'cahra':
        sys.path.insert(0, '/mnt/beegfs/elmlund/afan/cahra_13991'); import cahra_score as cs
        T['open'] = ('bin', (np.array(cs.read_star(tf)['conformation'], int)[:NM] == 1).astype(float)); junk = np.zeros(NM, bool)
    else:
        hd = [l.split()[1:] for l in open(tf) if l.startswith('#')][0]; M = np.genfromtxt(tf, comments='#'); tr = M[:, 0].astype(int) - 1
        def col(c): v = np.zeros(NM); v[tr] = M[:, hd.index(c)]; return v
        junk = col('junk') > 0; T['junk'] = ('bin', junk.astype(float))
        if kind == 'igm':
            for c in hd[3:]: T[c] = ('ang', col(c))
        else:
            jt = col('jtype'); T['j1'] = ('bin', (junk & (jt == 1)).astype(float)); T['j2'] = ('bin', (junk & (jt == 2)).astype(float))
            for c, nm in (('A_domain', 'A'), ('C_closed', 'C'), ('B_blob', 'B')): T[nm] = ('bin', ((col(c) > 0) & ~junk).astype(float))
    return T, junk
def graph(X, k):
    d, nb = cKDTree(X).query(X, k=k + 1, workers=4); d, nb = d[:, 1:], nb[:, 1:]; s = np.maximum(d[:, -1], 1e-12); n = len(X)
    w = np.exp(-d ** 2 / (s[:, None] * s[nb])); A = sp.csr_matrix((w.ravel(), (np.repeat(np.arange(n), k), nb.ravel())), shape=(n, n))
    W = A.maximum(A.T).tocsr(); deg = np.asarray(W.sum(1)).ravel(); return W, deg
def dmap(W, deg, m=10):
    Dm = sp.diags(1 / np.sqrt(deg)); S = (Dm @ W @ Dm).tocsr(); vals, vecs = eigsh(S, k=m + 1, which='LA', tol=1e-6, ncv=max(40, 3 * m), maxiter=20000)
    o = np.argsort(-vals); vals, vecs = vals[o], vecs[:, o]; psi = vecs / np.sqrt(deg)[:, None]; return vals[1:], psi[:, 1:]
def sweep(W, deg, f, minfrac=0.02):
    n = len(f); o = np.argsort(f, kind='stable'); pos = np.empty(n, int); pos[o] = np.arange(n); U = sp.triu(W, 1).tocoo()
    lo, hi = np.minimum(pos[U.row], pos[U.col]), np.maximum(pos[U.row], pos[U.col]); dif = np.zeros(n + 1); np.add.at(dif, lo, U.data); np.add.at(dif, hi, -U.data)
    cut = np.cumsum(dif)[:n - 1]; vp = np.cumsum(deg[o])[:n - 1]; tot = deg.sum(); sm = np.minimum(vp, tot - vp); ok = sm >= minfrac * tot
    phi = np.where(ok, cut / np.maximum(sm, 1e-12), np.inf); t = int(np.argmin(phi)); return float(phi[t]), float(sm[t] / tot)
def axis_profile(f, X, nb=20):
    o = np.argsort(f, kind='stable'); cen = np.array([X[i].mean(0) for i in np.array_split(o, nb)]); c0 = cen - cen.mean(0)
    u = np.linalg.svd(c0, full_matrices=False)[2][0]; u = u if np.corrcoef(X @ u, f)[0, 1] >= 0 else -u; return X @ u
def modes(x):
    lo, hi = np.percentile(x, [1, 99]); h = np.histogram(x, bins=30, range=(lo, hi))[0].astype(float); sm = np.convolve(np.r_[h[0], h, h[-1]], np.ones(3) / 3, 'valid')
    pk = [i for i in range(len(sm)) if (i == 0 or sm[i] >= sm[i - 1]) and (i == len(sm) - 1 or sm[i] > sm[i + 1]) and sm[i] >= 0.25 * sm.max()]
    if len(pk) < 2: return len(pk), 1.0
    a, b = sorted(sorted(pk, key=lambda i: -sm[i])[:2]); return len(pk), float(sm[a:b + 1].min() / min(sm[a], sm[b]))
def ks_shape(x, t):
    def sc(v): lo, hi = np.percentile(v, [1, 99]); return np.clip((v - lo) / max(hi - lo, 1e-12), 0, 1)
    g = np.linspace(0, 1, 201); a, b = np.sort(sc(x)), np.sort(sc(t)); return float(np.abs(np.searchsorted(a, g) / len(a) - np.searchsorted(b, g) / len(b)).max())
def run(tag, X, T, G, k, out, focus=None):
    n = len(X); W, deg = graph(X, k); vals, psi = dmap(W, deg)
    C = np.cov(X.T); Xn = RNG.multivariate_normal(np.zeros(X.shape[1]), C, size=n); Wn, dn = graph(Xn, k); _, psin = dmap(Wn, dn)
    phn = [sweep(Wn, dn, psin[:, j])[0] for j in range(psin.shape[1])]; pmin = min(phn); rows = []
    for j in range(psi.shape[1]):
        f = psi[:, j]; best = ('-', 0.)
        for nm, (kd, v) in T.items():
            if focus and nm != focus: continue
            m = G if kd == 'ang' else np.ones(n, bool); r = abs(spearmanr(f[m], v[m])[0]) if v[m].std() > 0 else 0.
            if r > best[1]: best = (nm, r)
        phi, bal = sweep(W, deg, f); phin = phn[j]; xp = axis_profile(f[G], X[G]); nm_, dip = modes(xp)
        ks = ks_shape(xp, T[best[0]][1][G]) if best[0] != '-' and T[best[0]][0] == 'ang' else np.nan
        rows.append((j + 1, vals[j], best[0], best[1], phi, phin, phi / phin, phi / pmin, bal, nm_, dip, ks))
    print('== %s: n %d, k %d, eigenvalues %s' % (tag, n, k, ' '.join('%.4f' % v for v in vals)))
    print('   psi  lambda  best-truth abs-rho | phi    phi_null ratio ratio/min side | modes dip  | KS(profile, truth)')
    for j, lam, bn, br, phi, phin, ra, rm, bal, nm_, dip, ks in rows:
        print('   %2d  %.4f  %-8s %.3f | %.4f %.4f  %.2f  %.2f  %.2f | %d %.2f | %s' % (j, lam, bn, br, phi, phin, ra, rm, bal, nm_, dip, '%.3f' % ks if ks == ks else '-'))
        out.write('%s\t%d\t%d\t%.5f\t%s\t%.4f\t%.5f\t%.5f\t%.3f\t%.3f\t%.3f\t%d\t%.3f\t%s\n' % (tag, k, j, lam, bn, br, phi, phin, ra, rm, bal, nm_, dip, '%.4f' % ks if ks == ks else 'nan'))
    if any(T[nm][0] == 'ang' for nm in T):
        print('   coordinate recovery (best psi per truth angle): %s' % ' '.join('%s->psi%d(%.2f)' % (nm, *max(((j + 1, abs(spearmanr(psi[G, j], T[nm][1][G])[0])) for j in range(psi.shape[1])), key=lambda x: x[1]))
                                                                             for nm in T if T[nm][0] == 'ang'))
    sys.stdout.flush(); return rows
out = open('dmap_test.tsv', 'w'); out.write('tag\tk\tpsi\tlambda\tbest_truth\tabs_rho\tphi\tphi_null\tratio\tratio_null_min\tsmaller_side\tmodes\tdip\tks_profile_vs_truth\n')
for nm, kind, tf, fe, ffe in SPECS:
    for lat, ef in (('flex', fe), ('fflex', ffe)):
        rows, z = read_emb(ef); NM = rows.max() + 1; T, junk = truth(kind, tf, NM); X = (z - z.mean(0)) / z.std(0)
        Tr = {k_: (kd, v[rows]) for k_, (kd, v) in T.items()}; G = ~junk[rows]
        for k in (15, 30): run('%s %s' % (nm, lat), X, Tr, G, k, out)
        if kind == 'igm' and lat == 'fflex':
            t = Tr['a4_th'][1]; tz = (t - t[G].mean()) / t[G].std(); gi = np.where(G)[0]
            two = 0.5 * np.exp(-(tz[gi] + 1) ** 2 / (2 * 0.35 ** 2)) + 0.5 * np.exp(-(tz[gi] - 1) ** 2 / (2 * 0.35 ** 2)); two /= two.max()
            skw = np.exp(1.2 * tz[gi]); skw /= skw.max(); u = RNG.random(len(gi))
            for sname, keep in (('two-peaked', gi[u < two]), ('skewed', gi[u < skw]), ('uniform', gi[u < two.mean()])):
                Xs = X[keep]; Xs = (Xs - Xs.mean(0)) / Xs.std(0); Ts = {k_: (kd, v[keep]) for k_, (kd, v) in Tr.items() if kd == 'ang'}
                ts = Ts['a4_th'][1]; tzs = (ts - t[G].mean()) / t[G].std(); h = np.histogram(tzs, bins=8, range=(-1.75, 1.75))[0]
                print('-- subsample %s: %d genuine particles; imposed a4_th occupancy (8 bins over +-1.75 sd): %s' % (sname, len(keep), ' '.join('%d' % x for x in h)))
                run('%s %s subsample %s' % (nm, lat, sname), Xs, Ts, np.ones(len(keep), bool), 15, out, focus='a4_th')
out.close()
