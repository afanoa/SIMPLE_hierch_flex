#!/usr/bin/env python3
"""Python toy: how should several INDEPENDENT motions be represented? Particles are noisy 2D 'maps' (P = G x G pixels) of a base plus
localized domains; each domain has a binary state (blob at c - delta / c + delta) or a continuous one (blob position along the swing).
Representations compared, all scored against the truth factors:
  TREE   greedy top-down splits (K-way to depth KWAY, binary below; held-out gain z > ZCRIT; minpop; maxdepth) + an equivalence merge-up,
         the hetero_tree logic in miniature;
  FLAT   one K-way partition (K by the same held-out test);
  TPF    TREE + post-hoc factorization: the cross-half reproducible subspace of the unit means, rotated to spatially independent components
         (spatial ICA), then a factorial re-assignment of every particle (per-factor coordinates; 2-means per factor) with pooled M-steps;
  DF     the same factorial model started from the particles' own cross-half PCA (no tree).
usage: mm_toy.py SCENARIO [--dprime 1.5 --n 60000 --seed 1]   SCENARIO: mm3 | mm2c | one | cont1 | coupled | partial"""
import sys, argparse, numpy as np
from sklearn.decomposition import FastICA
ap = argparse.ArgumentParser(); ap.add_argument('scen'); ap.add_argument('--dprime', type=float, default=1.5); ap.add_argument('--n', type=int, default=60000)
ap.add_argument('--seed', type=int, default=1); ap.add_argument('--grid', type=int, default=32); ap.add_argument('--minpop', type=int, default=900)
ap.add_argument('--maxdepth', type=int, default=4); ap.add_argument('--kway', type=int, default=2); ap.add_argument('--zcrit', type=float, default=3.0)
ap.add_argument('--quiet', action='store_true'); a = ap.parse_args(); rng = np.random.default_rng(a.seed); G = a.grid; P = G * G
yy, xx = np.mgrid[:G, :G].astype(float)
def blob(cx, cy, w=1.8): return np.exp(-0.5 * ((xx - cx) ** 2 + (yy - cy) ** 2) / w ** 2).ravel()
# domain sites on a ring; a binary domain swings its blob between two positions, a continuous one along the same arc
SITES = [(8, 8), (24, 9), (16, 24), (25, 23)]
SC = {'mm3': [('b', .5), ('b', .65), ('b', .8)], 'mm2c': [('b', .5), ('b', .7), ('c', 0)], 'one': [('b', .65)], 'cont1': [('C', 0)],
      'coupled': [('b', .65)] * 3, 'partial': [('b', .5), ('b', .65), ('b', .8)]}[a.scen]
def dom_map(k, s, kind):   # s in [0, 1]: 0 = state 1, 1 = state 2; 'C' = a long curved swing (a single motion spanning > 1 direction)
    cx, cy = SITES[k]; amp = 4.0 if kind != 'C' else 7.0; th = (2 * s - 1) * (0.9 if kind == 'C' else 0.6)
    return blob(cx + amp * np.sin(th), cy + amp * (1 - np.cos(th)) * (1 if kind == 'C' else 0) - amp * 0.0)
NF = len(SC); N = a.n
S = np.zeros((N, NF))
for k, (kind, p) in enumerate(SC):
    S[:, k] = (rng.random(N) < p).astype(float) if kind == 'b' else rng.random(N)
if a.scen == 'coupled': S[:, 1] = S[:, 0]; S[:, 2] = S[:, 0]   # three domains move together: ONE motion
if a.scen == 'partial': keep = (S[:, 0] + S[:, 1] + S[:, 2]) % 2 == 0; S = S[keep]; N = len(S)   # only 4 of 8 combinations exist
base = np.zeros(P)
for cx, cy in [(16, 16), (12, 14), (20, 18)]: base += 1.5 * blob(cx, cy, 4.)
X = np.tile(base, (N, 1))
for k, (kind, p) in enumerate(SC):
    dk = 'C' if kind == 'C' else kind
    if kind == 'b': X += np.where(S[:, k:k + 1] > 0.5, dom_map(k, 1, dk), dom_map(k, 0, dk))
    else:
        lv = np.linspace(0, 1, 41); M = np.array([dom_map(k, t, dk) for t in lv]); X += M[np.rint(S[:, k] * 40).astype(int)]
# noise: per-particle d' of a binary domain's two states = DPRIME
d01 = np.linalg.norm(dom_map(0, 1, 'b') - dom_map(0, 0, 'b')); sig = d01 / a.dprime
Y = (X + rng.normal(0, sig, X.shape)).astype(np.float32); del X
# low-pass (the tree's lp): keep the orthonormal 2D DCT coefficients with kx^2 + ky^2 <= (G/LP)^2; white noise stays sigma per coefficient
from scipy.fft import dctn, idctn
ky, kx = np.mgrid[:G, :G]; KEEP = (kx ** 2 + ky ** 2 <= (G / 3.) ** 2).ravel()
def to_coef(Z): return dctn(Z.reshape(-1, G, G), axes=(1, 2), norm='ortho').reshape(len(Z), -1)[:, KEEP]
def to_pix(Cf): Z = np.zeros((len(Cf), P)); Z[:, KEEP] = Cf; return idctn(Z.reshape(-1, G, G), axes=(1, 2), norm='ortho').reshape(len(Cf), -1)
Y = to_coef(Y).astype(np.float32); P_eff = Y.shape[1]
half = rng.random(N) < 0.5
truth_combo = np.array([int(''.join('2' if S[i, k] > 0.5 else '1' for k in range(NF))) for i in range(N)])
def log(*x):
    if not a.quiet: print(*x, flush=True)
# ---------------- mixture fitting (spherical, known sigma: the likelihood is -|y - m|^2 / 2 sig^2; |y|^2 cancels between models) ----------------
def sqd(Yn, M): return (M ** 2).sum(1)[None] - 2. * (Yn @ M.T)
def pcs(Yn, q=8):
    sub = Yn[rng.choice(len(Yn), min(len(Yn), 20000), replace=False)]; sub = sub - sub.mean(0)
    return np.linalg.svd(sub, full_matrices=False)[2][:q]
def em(Yn, K, it=25, restarts=3, init=None):
    best = None
    for r in range(restarts if init is None else 1):
        if init is None:   # seeds from k-means++ in the node's top principal components (the flex_pca latent's role)
            Z = (Yn - Yn.mean(0)) @ pcs(Yn).T; c = [Z[rng.integers(len(Z))]]
            for _ in range(K - 1):
                d = np.maximum(np.min([((Z - ci) ** 2).sum(1) for ci in c], 0), 1e-12); c.append(Z[rng.choice(len(Z), p=d / d.sum())])
            c = np.array(c)
            for _ in range(20): lz = ((Z[:, None] - c[None]) ** 2).sum(2).argmin(1); c = np.array([Z[lz == k].mean(0) if (lz == k).any() else c[k] for k in range(K)])
            M = np.array([Yn[lz == k].mean(0) if (lz == k).any() else Yn.mean(0) for k in range(K)], np.float64)
        else: M = init.copy()
        w = np.full(K, 1. / K)
        for _ in range(it):
            L = -0.5 * sqd(Yn, M) / sig ** 2 + np.log(w)[None]; L0 = L.copy()
            L -= L.max(1, keepdims=True); R = np.exp(L); R /= R.sum(1, keepdims=True)
            nk = R.sum(0) + 1e-9; M = (R.T @ Yn) / nk[:, None]; w = nk / nk.sum()
        mx = L0.max(1); ll = (mx + np.log(np.exp(L0 - mx[:, None]).sum(1))).sum()
        if best is None or ll > best[0]: best = (ll, M, w)
    return best[1], best[2]
def loglik(Yn, M, w):
    L = -0.5 * sqd(Yn, M) / sig ** 2 + np.log(w)[None]; mx = L.max(1); return mx + np.log(np.exp(L - mx[:, None]).sum(1))
def heldout_gain(idx, K, M1=None):
    """fit K on the A half, per-particle log-lik gain over K-1 (and 1) on the B half; z of the gain over K-1"""
    A = idx[half[idx]]; Bh = idx[~half[idx]]
    MK, wK = em(Y[A], K)
    Mk1, wk1 = (Y[A].mean(0, keepdims=True), np.ones(1)) if K == 2 else em(Y[A], K - 1)
    g = loglik(Y[Bh], MK, wK) - loglik(Y[Bh], Mk1, wk1); z = g.mean() / (g.std() / np.sqrt(len(g)) + 1e-12)
    return z, MK, wK
def assign(idx, M, w):
    L = -0.5 * sqd(Y[idx], M) / sig ** 2 + np.log(w)[None]; return L.argmax(1)
def split_node(idx, kmax):
    """K* = the largest K (<= kmax) whose step from K-1 passes; children re-fitted on all node rows from the A-half means"""
    Ks, last, zs = 1, None, []
    for K in range(2, kmax + 1):
        if len(idx) < K * a.minpop: break
        z, MK, wK = heldout_gain(idx, K); zs.append('%d:%.1f' % (K, z))
        if z < a.zcrit: break
        Ks, last = K, (MK, wK)
    ZLOG.append((len(idx), ' '.join(zs)))
    if Ks == 1: return None
    M, w = em(Y[idx], Ks, it=15, init=last[0]); lab = assign(idx, M, w)
    ch = [idx[lab == k] for k in range(Ks)]
    if min(len(c) for c in ch) < a.minpop: return None
    return ch
ZLOG = []
def run_tree():
    leaves, splits = [], []
    def rec(idx, name, depth):
        ch = None if depth >= a.maxdepth else split_node(idx, 5 if depth < a.kway else 2)
        if ch is None: leaves.append((name, idx)); return
        splits.append((name, [Y[c].mean(0) for c in ch], [len(c) for c in ch], ch))
        for k, c in enumerate(ch): rec(c, name + str(k + 1), depth + 1)
    rec(np.arange(N), 'r', 0); return leaves, splits
def mergeup(cells):
    """fuse cell pairs whose held-out 2-class separation is NOT significant (|z| < 1), most similar first"""
    cells = [np.asarray(c) for c in cells]; changed = True
    while changed and len(cells) > 1:
        changed = False; mu = np.array([Y[c].mean(0) for c in cells]); best = None
        for i in range(len(cells)):
            for j in range(i + 1, len(cells)):
                A = np.r_[cells[i][half[cells[i]]], cells[j][half[cells[j]]]]; Bh = np.r_[cells[i][~half[cells[i]]], cells[j][~half[cells[j]]]]
                yb = np.r_[np.zeros((~half[cells[i]]).sum()), np.ones((~half[cells[j]]).sum())]
                mA = [Y[cells[i][half[cells[i]]]].mean(0), Y[cells[j][half[cells[j]]]].mean(0)]; dv = mA[1] - mA[0]
                s = (Y[Bh] - 0.5 * (mA[0] + mA[1])) @ dv / sig ** 2; g = np.where(yb > 0, s, -s)   # held-out log-odds gain
                z = g.mean() / (g.std() / np.sqrt(len(g)) + 1e-12)
                if z < 1 and (best is None or z < best[0]): best = (z, i, j)
        if best: _, i, j = best; cells[i] = np.r_[cells[i], cells[j]]; cells.pop(j); changed = True
    return cells
# ---------------- factorization ----------------
def xhalf_rank(MA, MB, wts):
    """cross-half reproducible rank of centred cell means: eigenvalues of the symmetrised cross-half Gram above 3x the |most negative|"""
    cA = MA - (wts[:, None] * MA).sum(0) / wts.sum(); cB = MB - (wts[:, None] * MB).sum(0) / wts.sum()
    Gm = np.sqrt(wts)[:, None] * (cA @ cB.T) * np.sqrt(wts)[None]; ev = np.linalg.eigvalsh(0.5 * (Gm + Gm.T))
    thr = 3 * max(abs(ev.min()), 1e-12); return int((ev > thr).sum()), ev
def spatial_ica(C, r):
    """r spatially independent difference patterns from the rows of C (centred means): SVD subspace, then FastICA over pixels"""
    U, s, Vt = np.linalg.svd(C, full_matrices=False); B = Vt[:r]
    if r == 1: return B
    ica = FastICA(n_components=r, whiten='unit-variance', random_state=0, max_iter=2000); Sx = ica.fit_transform(to_pix(B).T)   # pixels x r
    D = to_coef(Sx.T); D /= np.linalg.norm(D, axis=1, keepdims=True); return D
def total_corr(pi, r):
    Sc = np.array([[(c >> k) & 1 for k in range(r)] for c in range(2 ** r)]); m = np.array([pi[Sc[:, k] == 1].sum() for k in range(r)])
    prod = np.prod(np.where(Sc == 1, m[None], 1 - m[None]), 1); return float((pi * np.log((pi + 1e-12) / (prod + 1e-12))).sum())
def group_components(Q, T):
    """motions = connected groups of components whose labels or coordinates are dependent (label NMI > 0.05 or |coordinate r| > 0.2);
    independent motions give ~0 on both (toy: NMI <= 0.004, |r| <= 0.09), one curved motion spread over components gives 0.15+ / 0.3+"""
    r = Q.shape[1]; adj = np.eye(r, dtype=bool); cr = np.corrcoef(T.T) if r > 1 else np.ones((1, 1))
    for i in range(r):
        for j in range(i + 1, r):
            if nmi(Q[:, i] > .5, Q[:, j] > .5) > 0.05 or abs(cr[i, j]) > 0.2: adj[i, j] = adj[j, i] = True
    grp = -np.ones(r, int); g = 0
    for i in range(r):
        if grp[i] >= 0: continue
        stack = [i]; grp[i] = g
        while stack:
            u = stack.pop()
            for v in np.where(adj[u])[0]:
                if grp[v] < 0: grp[v] = g; stack.append(v)
        g += 1
    return grp
def group_scores(Q, T, grp, R):
    """per truth factor: the best motion group's joint-state labels (combinations of its components, MAP from R marginalised)"""
    r = Q.shape[1]; Sc = np.array([[(c >> k) & 1 for k in range(r)] for c in range(2 ** r)]); out = []
    glab = {}
    for g in np.unique(grp):
        cols = np.where(grp == g)[0]; code = (Sc[:, cols] * (2 ** np.arange(len(cols)))[None]).sum(1)
        P = np.zeros((N, 2 ** len(cols)))
        for c in range(2 ** r): P[:, code[c]] += R[:, c]
        glab[g] = P.argmax(1)
    for k, (kind, p) in enumerate(SC):
        if a.scen == 'coupled' and k > 0: break
        if kind == 'b':
            st = S[:, k] > 0.5; best = max((sum(max(st[glab[g] == c].sum(), (~st[glab[g] == c]).sum()) for c in np.unique(glab[g])) / N, g) for g in glab)
            out.append('f%d acc %.3f (motion %d)' % (k + 1, *best))
        else:
            t = S[:, k]; res = []
            for g, lab in glab.items():
                u = np.unique(lab); mu = np.array([t[lab == c].mean() for c in u]); cnt = np.array([(lab == c).sum() for c in u])
                eta = (cnt * (mu - t.mean()) ** 2).sum() / ((t - t.mean()) ** 2).sum(); lo, hi = np.percentile(t, [10, 90])
                tf = [max(2. * ((t <= lo) & (lab == c)).sum() / ((lab == c).sum() + (t <= lo).sum()) for c in u), max(2. * ((t >= hi) & (lab == c)).sum() / ((lab == c).sum() + (t >= hi).sum()) for c in u)]
                res.append((eta, g, tf))
            eta, g, tf = max(res); out.append('f%d eta2 %.3f tailF1 %.2f/%.2f (motion %d, %d states)' % (k + 1, eta, *tf, g, len(np.unique(glab[g]))))
    return out
def fact_ll(Yn, mu, D, pi):
    r = len(D); Sc = np.array([[(c >> k) & 1 for k in range(r)] for c in range(2 ** r)], float); Mc = mu[None] + Sc @ D
    Lg = -0.5 * sqd(Yn, Mc) / sig ** 2 + np.log(pi + 1e-300)[None]; mx = Lg.max(1); return mx + np.log(np.exp(Lg - mx[:, None]).sum(1))
def prune(D0, label):
    """drop factors that do not raise the held-out likelihood: fit the factorial model on half A with and without each factor,
    per-particle log-lik gain on half B, keep a factor only if its z >= ZCRIT; the weakest goes first, then refit"""
    D = D0.copy(); A = np.where(half)[0]; Bh = np.where(~half)[0]
    while len(D) > 1:
        _, _, Dk, pik, _, muk = factorial(D, rows=A, ret_mu=True); base = fact_ll(Y[Bh], muk, Dk, pik); zs = []
        for j in range(len(D)):
            keep = [q for q in range(len(D)) if q != j]; _, _, Dj, pij, _, muj = factorial(D[keep], rows=A, ret_mu=True)
            g = base - fact_ll(Y[Bh], muj, Dj, pij); zs.append(g.mean() / (g.std() / np.sqrt(len(g)) + 1e-12))
        j = int(np.argmin(zs)); log('%s held-out factor z: %s' % (label, ' '.join('%.1f' % z for z in zs)))
        if zs[j] >= a.zcrit: break
        D = np.delete(D, j, 0)
    return D
def factorial(D0, iters=15, rows=None, ret_mu=False):
    """additive factorial mixture: y = mu + sum_k s_k d_k + noise, s_k in {0, 1}, a FREE joint prior over the 2^r combinations (coupled
    motions allowed), exact E-step over the combinations, LS M-step for mu and D from E[s], E[s s^T]; started from the patterns D0 with
    per-factor 2-means on the LS coordinates. Returns the LS coordinates T, the per-factor marginal posteriors Q (N x r) and D."""
    Yf = Y if rows is None else Y[rows]; Nf = len(Yf)
    r = len(D0); mu = Yf.mean(0); T = np.linalg.lstsq(D0.T, (Yf - mu).T, rcond=None)[0].T; L = np.zeros_like(T)
    for k in range(r):
        t = T[:, k]; c = np.percentile(t, [25, 75])
        for _ in range(20): lab = np.abs(t - c[1]) < np.abs(t - c[0]); c = np.array([t[~lab].mean(), t[lab].mean()])
        L[:, k] = lab
    coef = np.linalg.lstsq(np.c_[np.ones(Nf), L], Yf, rcond=None)[0]; mu, D = coef[0], coef[1:]
    Sc = np.array([[(c >> k) & 1 for k in range(r)] for c in range(2 ** r)], float); pi = np.full(2 ** r, 1. / 2 ** r)
    for _ in range(iters):
        Mc = mu[None] + Sc @ D; Lg = -0.5 * sqd(Yf, Mc) / sig ** 2 + np.log(pi + 1e-300)[None]
        Lg -= Lg.max(1, keepdims=True); R = np.exp(Lg); R /= R.sum(1, keepdims=True)
        pi = R.mean(0); Es = R @ Sc; Ess = Sc.T @ (R.sum(0)[:, None] * Sc)   # E[s] per particle, sum E[s s^T]
        A = np.zeros((r + 1, r + 1)); A[0, 0] = Nf; A[0, 1:] = A[1:, 0] = Es.sum(0); A[1:, 1:] = Ess
        rhs = np.r_[Yf.sum(0)[None], Es.T @ Yf]; coef = np.linalg.solve(A, rhs); mu, D = coef[0], coef[1:]
    Q = R @ Sc; T = np.linalg.lstsq(D.T, (Yf - mu).T, rcond=None)[0].T
    return (T, Q, D, pi, R, mu) if ret_mu else (T, Q, D, pi, R)
# ---------------- scoring ----------------
def nmi(x, y):
    x = np.unique(x, return_inverse=True)[1]; y = np.unique(y, return_inverse=True)[1]; p = np.zeros((x.max() + 1, y.max() + 1)); np.add.at(p, (x, y), 1); p /= p.sum()
    px, py = p.sum(1), p.sum(0); nz = p > 0; mi = (p[nz] * np.log(p[nz] / (px[:, None] * py[None])[nz])).sum()
    hx = -(px[px > 0] * np.log(px[px > 0])).sum(); hy = -(py[py > 0] * np.log(py[py > 0])).sum(); return mi / max(np.sqrt(hx * hy), 1e-12)
def labels_from_cells(cells):
    lab = np.zeros(N, int)
    for k, c in enumerate(cells): lab[c] = k + 1
    return lab
def factor_scores(lab):
    """per factor: leaf-majority accuracy (binary) or eta^2 (continuous); combination NMI; best single-cell F1 per binary factor state"""
    out = []
    for k, (kind, p) in enumerate(SC):
        if a.scen == 'coupled' and k > 0: break
        if kind == 'b':
            st = S[:, k] > 0.5; acc = sum(max(st[lab == c].sum(), (~st[lab == c]).sum()) for c in np.unique(lab)) / N
            f1s = []
            for s_ in (False, True):
                tg = st == s_; f1s.append(max(2. * (tg & (lab == c)).sum() / ((lab == c).sum() + tg.sum()) for c in np.unique(lab)))
            out.append('f%d acc %.3f cellF1 %.2f/%.2f' % (k + 1, acc, *f1s))
        else:
            t = S[:, k]; mu = np.array([t[lab == c].mean() for c in np.unique(lab)]); cnt = np.array([(lab == c).sum() for c in np.unique(lab)])
            eta = (cnt * (mu - t.mean()) ** 2).sum() / ((t - t.mean()) ** 2).sum(); lo, hi = np.percentile(t, [10, 90])
            tf = [max(2. * ((t <= lo) & (lab == c)).sum() / ((lab == c).sum() + (t <= lo).sum()) for c in np.unique(lab)),
                  max(2. * ((t >= hi) & (lab == c)).sum() / ((lab == c).sum() + (t >= hi).sum()) for c in np.unique(lab))]
            out.append('f%d eta2 %.3f tailF1 %.2f/%.2f' % (k + 1, eta, *tf))
    return out
def contrast(groups_by_factor):
    """map-level contrast of factor k between two delivered maps A, B: |frac(state 2 in A) - frac(state 2 in B)| x sqrt(nA nB / (nA + nB))
    (the domain difference's SNR in the difference of the two maps, up to d'); groups_by_factor[k] = list of index arrays, best pair"""
    out = []
    for k, groups in enumerate(groups_by_factor):
        if SC[k][0] != 'b' or (a.scen == 'coupled' and k > 0): out.append(np.nan); continue
        st = S[:, k] > 0.5; fr = [(st[g].mean(), len(g)) for g in groups if len(g) > 0]; best = 0.
        for i in range(len(fr)):
            for j in range(i + 1, len(fr)):
                best = max(best, abs(fr[i][0] - fr[j][0]) * np.sqrt(fr[i][1] * fr[j][1] / (fr[i][1] + fr[j][1])))
        out.append(best)
    return out
def cells_contrast(lab): return contrast([[np.where(lab == c)[0] for c in np.unique(lab)]] * NF)
def factorial_scores(T, Q, D):
    """match each truth factor to its best component: accuracy of the marginal posterior > 0.5 / R^2 and tail purity of the coordinate;
    plus the factor-state maps' contrast (each factor's two pooled groups)"""
    out = []; groups = []
    for k, (kind, p) in enumerate(SC):
        if a.scen == 'coupled' and k > 0: groups.append([]); break
        if kind == 'b':
            st = S[:, k] > 0.5; accs = [max(((Q[:, j] > .5) == st).mean(), ((Q[:, j] > .5) != st).mean()) for j in range(Q.shape[1])]; j = int(np.argmax(accs))
            groups.append([np.where(Q[:, j] > .5)[0], np.where(Q[:, j] <= .5)[0]]); out.append('f%d acc %.3f (comp %d)' % (k + 1, accs[j], j + 1))
        else:
            t = S[:, k]; Xd = np.c_[np.ones(N), T]; r2s = [np.corrcoef(T[:, j], t)[0, 1] ** 2 for j in range(T.shape[1])]
            res = t - Xd @ np.linalg.lstsq(Xd, t, rcond=None)[0]; r2all = 1 - res.var() / t.var()
            lo, hi = np.percentile(t, [10, 90]); j = int(np.argmax(r2s)); tj = T[:, j] * np.sign(np.corrcoef(T[:, j], t)[0, 1])
            n10 = int(0.1 * N); tl = np.argsort(tj)[:n10]; th = np.argsort(tj)[-n10:]; groups.append([])
            out.append('f%d R2 best comp %.3f all comps %.3f tail purity %.2f/%.2f' % (k + 1, r2s[j], r2all, (t[tl] <= lo).mean(), (t[th] >= hi).mean()))
    cs = contrast(groups + [[]] * (NF - len(groups)))
    out.append('map contrast ' + ' '.join('%.0f' % c for c in cs if np.isfinite(c)))
    return out
log('== scenario %s (%d low-pass coefficients): %d particles, %d factors %s, d\' %.2f (sigma %.3f), grid %d, minpop %d, maxdepth %d, K-way to depth %d, zcrit %.1f' %
    (a.scen, P_eff, N, NF, SC, a.dprime, sig, G, a.minpop, a.maxdepth, a.kway, a.zcrit))
leaves, splits = run_tree(); log('node tests (rows: K:z): ' + ' | '.join('%d %s' % zl for zl in ZLOG)); lab_t = labels_from_cells([c for _, c in leaves])
units = mergeup([c for _, c in leaves]); lab_u = labels_from_cells(units)
rep = {}
for name, means, cnts, ch in splits:
    labs = np.zeros(N, int)
    for k, c in enumerate(ch): labs[c] = k + 1
    m = labs > 0; Vs = []
    for k in range(NF):
        st = S[m, k] > 0.5 if SC[k][0] == 'b' else S[m, k] > np.median(S[:, k])
        Vs.append(nmi(labs[m], st))
    fol = [k + 1 for k in range(NF) if Vs[k] > 0.1]
    log('split %-6s K %d sizes %s  NMI per factor %s -> follows %s' % (name, len(ch), cnts, ' '.join('%.2f' % v for v in Vs), fol))
    for f in fol: rep.setdefault(f, []).append(name)
deep = {f: [n for n in ns if not any(o != n and o.startswith(n) for o in ns)] for f, ns in rep.items()}
log('TREE  leaves %d units %d; leaf sizes %s; combo NMI leaves %.3f units %.3f; repeats %s' % (len(leaves), len(units), sorted([len(c) for _, c in leaves])[::-1],
    nmi(lab_t, truth_combo), nmi(lab_u, truth_combo), ' '.join('f%d:%d' % (f, len(v)) for f, v in sorted(deep.items()))))
log('TREE  leaves: ' + '; '.join(factor_scores(lab_t)) + '; map contrast ' + ' '.join('%.0f' % c for c in cells_contrast(lab_t) if np.isfinite(c)))
log('TREE  units : ' + '; '.join(factor_scores(lab_u)) + '; map contrast ' + ' '.join('%.0f' % c for c in cells_contrast(lab_u) if np.isfinite(c)))
# FLAT
idx = np.arange(N); Ks, last = 1, None
for K in range(2, 13):
    z, MK, wK = heldout_gain(idx, K)
    if z < a.zcrit: break
    Ks, last = K, (MK, wK)
lab_f = assign(idx, *em(Y, Ks, it=15, init=last[0])) + 1 if Ks > 1 else np.ones(N, int)
log('FLAT  K %d; combo NMI %.3f; %s; map contrast %s' % (Ks, nmi(lab_f, truth_combo), '; '.join(factor_scores(lab_f)), ' '.join('%.0f' % c for c in cells_contrast(lab_f) if np.isfinite(c))))
# TPF: factorize the unit means
wts = np.array([len(c) for c in units], float)
MA = np.array([Y[c[half[c]]].mean(0) for c in units]); MB = np.array([Y[c[~half[c]]].mean(0) for c in units]); Mall = np.array([Y[c].mean(0) for c in units])
r, ev = xhalf_rank(MA, MB, wts); log('TPF   %d units -> cross-half rank %d (top eigenvalues %s)' % (len(units), r, ' '.join('%.3g' % e for e in ev[::-1][:6])))
if r >= 1:
    C = np.sqrt(wts)[:, None] * (Mall - (wts[:, None] * Mall).sum(0) / wts.sum()); D0 = prune(spatial_ica(C, r), 'TPF  '); r = len(D0); T, Q, D, pi, R = factorial(D0)
    lab_c = R.argmax(1)
    log('TPF   factors %d; combo NMI %.3f; %s' % (r, nmi(lab_c, truth_combo), '; '.join(factorial_scores(T, Q, D))))
    grp = group_components(Q, T); log('TPF   motions (dependent components grouped): %d %s; %s' % (grp.max() + 1, grp.tolist(), '; '.join(group_scores(Q, T, grp, R))))
    if r >= 2: log('TPF   joint-prior total correlation %.4f nats (0 = independent factors)' % total_corr(pi, r)); cc = np.corrcoef(T.T); log('TPF   coordinate correlations between factors: max |r| %.2f; label NMI between factors max %.3f; joint prior %s' % (np.abs(cc - np.eye(r)).max(), max(nmi(Q[:, i] > .5, Q[:, j] > .5) for i in range(r) for j in range(i + 1, r)), ' '.join('%.3f' % x for x in pi)))
# DF: factorial from the particles' own cross-half PCA (components from half A, kept when their variance replicates in half B)
mu = Y.mean(0); Ya = Y[half] - mu; Yb = Y[~half] - mu
V2 = np.linalg.svd(Ya[rng.choice(len(Ya), min(len(Ya), 30000), replace=False)], full_matrices=False)[2][:20]
vb = ((Yb @ V2.T) ** 2).mean(0); zb = (vb - sig ** 2) / (sig ** 2 * np.sqrt(2. / len(Yb))); rdf = int(np.argmax(zb < a.zcrit)) if (zb < a.zcrit).any() else 20
if rdf >= 1:
    B = V2[:rdf]; D0 = prune(spatial_ica(B, rdf) if rdf > 1 else B, 'DF   '); rdf = len(D0); T2, Q2, D2, pi2, R2 = factorial(D0)
    lab_d = R2.argmax(1)
    log('DF    rank %d (half-B excess variance z %s; total corr %.4f); combo NMI %.3f; %s' % (rdf, ' '.join('%.1f' % x for x in zb[:rdf + 2]), total_corr(pi2, rdf), nmi(lab_d, truth_combo), '; '.join(factorial_scores(T2, Q2, D2))))
    grp2 = group_components(Q2, T2); log('DF    motions: %d %s; %s' % (grp2.max() + 1, grp2.tolist(), '; '.join(group_scores(Q2, T2, grp2, R2))))
else: log('DF    rank 0 (z %s)' % ' '.join('%.1f' % x for x in zb[:4]))
