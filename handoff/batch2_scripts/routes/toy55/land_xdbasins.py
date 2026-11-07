#!/usr/bin/env python3
"""land_xdbasins.py: BASINS of the DECONVOLVED flex latent density (no path, no state count).
From a run's flex_pca_embedding.bin: raw MAP latents z_i, their stored precision Q_i, the eigenvalues (prior P = 1/eig) and the
calibrated noise scale a. SIMPLE's own deconvolution model (simple_flex_pca_deconv: R_i = Q_i^-1 D_i, N_i = a Q_i^-1 D_i Q_i^-1,
D_i = Q_i - P) is used in its equivalent ML form x_i = D_i^-1 Q_i z_i with noise a D_i^-1, and a K-Gaussian mixture is fitted through
it by extreme deconvolution (Bovy et al. 2011; K by particle-half held-out likelihood, or SIMPLE's K). The basins are the modes of the
fitted (deconvolved) density, found by mean shift from every component and from samples; two modes join at the best (maximin) path
between them on a kNN graph of density samples (edges scored at interior points), and the weakest join below THR kT merges, repeatedly.
Masses from density samples; a particle goes to the basin its posterior mean flows to. A split counts only if each half (eo), fitted
alone, gives the same number of basins with matching membership (ARI >= 0.8 on its particles).
Modes: --gt G --codes 200:230 (imposed shapes along a known coordinate), --gt G --discrete (contingency with truth groups),
or runs (--emb A,B,C --names a,b,c; --labout PREFIX --nm N writes one basin per master row for land_rec.sbatch)."""
import os, argparse, numpy as np
from multiprocessing import Pool
from scipy.spatial import cKDTree
from scipy.special import logsumexp
ap = argparse.ArgumentParser(); ap.add_argument('--emb', required=True); ap.add_argument('--names', default='')
ap.add_argument('--k', default='select', help="'select' (held-out halves), 'simple' (the run's flex_pca_deconv_prior.txt K) or an integer")
ap.add_argument('--kmax', type=int, default=10); ap.add_argument('--thr', default='0.5,1.0,1.5,2.0,3.0'); ap.add_argument('--gt', default='')
ap.add_argument('--codes', default='200:230'); ap.add_argument('--discrete', action='store_true'); ap.add_argument('--shapes', default='sim,flat,broad,dw2,dw4,shoulder')
ap.add_argument('--seed', type=int, default=0); ap.add_argument('--exclude', default=''); ap.add_argument('--labout', default='')
ap.add_argument('--labthr', default='1.0,2.0'); ap.add_argument('--nm', type=int, default=0); ap.add_argument('--nproc', type=int, default=8)
ap.add_argument('--maxit', type=int, default=150); ap.add_argument('--nsamp', type=int, default=20000); a = ap.parse_args()
def log(*s): print(*s, flush=True)
THR = [float(x) for x in a.thr.split(',')]; L2P = np.log(2 * np.pi)

def read_emb(path):
    """raw block (pinds, MAP z, eigenvalues, stored precision) and the trailing block's noise scale / deconvolved z"""
    b = open(path, 'rb').read(); n, D = np.frombuffer(b, np.int32, 2, 12); n, D = int(n), int(D); o = 28
    pind = np.frombuffer(b, np.int32, n, o).copy(); o += 4 * n
    z = np.frombuffer(b, np.float64, n * D, o).reshape(D, n).T.copy(); o += 8 * n * D
    eig = np.frombuffer(b, np.float64, D, o).copy(); o += 8 * D + 24 * n
    Q = np.frombuffer(b, np.float64, n * D * D, o).reshape(n, D, D).copy(); o += 8 * n * D * D + 8
    ns, zd = None, None
    if len(b) >= o + 20 and b[o:o + 8] == b'SIMPLFXD':
        ns = float(np.frombuffer(b, np.float64, 1, o + 8)[0])
        if int(np.frombuffer(b, np.int32, 1, o + 16)[0]) == 1: zd = np.frombuffer(b, np.float64, n * D, o + 28).reshape(D, n).T.copy()
    return pind, z, eig, 0.5 * (Q + np.transpose(Q, (0, 2, 1))), ns, zd
def ml_form(z, Q, eig, ns):
    """SIMPLE's noise model in ML form: x = D^-1 Q z, noise a D^-1 (D = Q - P, eigenvalues floored)"""
    P = 1. / np.maximum(eig, 1e-300); Dm = Q - np.diag(P)[None]; w, V = np.linalg.eigh(Dm); fl = 1e-3 * P.min(); bad = (w < fl).any(1)
    Di = (V / np.maximum(w, fl)[:, None, :]) @ np.transpose(V, (0, 2, 1)); x = (Di @ (Q @ z[..., None]))[..., 0]
    return x, ns * Di, bad.mean()

def kmeanspp(X, K, r):
    C = [X[r.integers(len(X))]]
    for _ in range(1, K):
        d2 = np.min(((X[:, None, :] - np.array(C)[None]) ** 2).sum(2), 1); C.append(X[r.choice(len(X), p=d2 / d2.sum())])
    C = np.array(C)
    for _ in range(20):
        lab = np.argmin(((X[:, None, :] - C[None]) ** 2).sum(2), 1); C = np.array([X[lab == k].mean(0) if (lab == k).any() else C[k] for k in range(K)])
    return C, lab
def xd_fit(X, S, K, Xinit, maxit=150, tol=1e-6, seed=0, pimin=1e-4):
    """extreme deconvolution, x_i = v_i + noise(S_i), v ~ sum_k pi_k N(mu_k, Sig_k); init by k-means++ on Xinit (posterior means)"""
    n, D = X.shape; r = np.random.default_rng(seed); sub = r.choice(n, min(n, 20000), replace=False); C, lab = kmeanspp(Xinit[sub], K, r)
    pi = np.array([max((lab == k).mean(), pimin) for k in range(K)]); pi /= pi.sum(); mu = C.copy(); vtot = np.cov(Xinit.T) + 1e-6 * np.eye(D)
    Sig = np.array([np.cov(Xinit[sub][lab == k].T) + 0.05 * vtot if (lab == k).sum() > D + 1 else vtot.copy() for k in range(K)]); llo = -np.inf
    for it in range(maxit):
        lp = np.empty((n, K)); Ti = []; Tid = []
        for k in range(K):
            T = S + Sig[k][None]; Lc = np.linalg.cholesky(T); Ci = np.linalg.inv(Lc); Tk = np.transpose(Ci, (0, 2, 1)) @ Ci
            d = X - mu[k]; td = (Tk @ d[..., None])[..., 0]; lp[:, k] = np.log(pi[k]) - 0.5 * ((d * td).sum(1) + 2 * np.log(np.diagonal(Lc, axis1=1, axis2=2)).sum(1) + D * L2P)
            Ti.append(Tk); Tid.append(td)
        ll = logsumexp(lp, 1); R = np.exp(lp - ll[:, None]); tot = ll.sum()
        for k in range(K):
            rk = R[:, k]; Nk = rk.sum()
            if Nk <= pimin * n: pi[k] = pimin; continue
            bk = mu[k] + Tid[k] @ Sig[k]; mnew = (rk[:, None] * bk).sum(0) / Nk; db = bk - mnew; TW = np.einsum('n,nij->ij', rk, Ti[k])
            Sig[k] = (np.einsum('n,ni,nj->ij', rk, db, db) + Nk * Sig[k] - Sig[k] @ TW @ Sig[k]) / Nk; Sig[k] = 0.5 * (Sig[k] + Sig[k].T) + 1e-6 * np.trace(vtot) / D * np.eye(D)
            mu[k] = mnew; pi[k] = Nk / n
        pi = np.maximum(pi, pimin); pi /= pi.sum()
        if it > 5 and abs(tot - llo) <= tol * abs(tot): break
        llo = tot
    post = sum(R[:, k][:, None] * (mu[k] + Tid[k] @ Sig[k]) for k in range(K))   # posterior means (last E-step)
    return dict(pi=pi, mu=mu, Sig=Sig, ll=tot / n, it=it + 1, post=post, resp=R)
def xd_score(X, S, f):
    lp = np.empty((len(X), len(f['pi'])))
    for k in range(len(f['pi'])):
        T = S + f['Sig'][k][None]; sgn, ld = np.linalg.slogdet(T); d = X - f['mu'][k]; sol = np.linalg.solve(T, d[..., None])[..., 0]
        lp[:, k] = np.log(f['pi'][k]) - 0.5 * ((d * sol).sum(1) + ld + X.shape[1] * L2P)
    return logsumexp(lp, 1).sum()
def _ladder(arg):
    K, XA, SA, IA, XB, SB, IB, maxit = arg; fa = xd_fit(XA, SA, K, IA, maxit=maxit); fb = xd_fit(XB, SB, K, IB, maxit=maxit)
    return K, xd_score(XB, SB, fa) + xd_score(XA, SA, fb)
def select_k(X, S, Xi, eo, kmax, maxit, seed=0):
    r = np.random.default_rng(seed); idx = [r.permutation(np.where(eo == h)[0])[:10000] for h in (0, 1)]
    args = [(K, X[idx[0]], S[idx[0]], Xi[idx[0]], X[idx[1]], S[idx[1]], Xi[idx[1]], maxit) for K in range(1, kmax + 1)]
    with Pool(min(a.nproc, kmax)) as p: sc = dict(p.map(_ladder, args))
    kb = max(sc, key=sc.get); return kb, sc

def gmm_lp(X, f):
    Si = np.linalg.inv(f['Sig']); ld = np.linalg.slogdet(f['Sig'])[1]; lp = np.empty((len(X), len(f['pi'])))
    for k in range(len(f['pi'])):
        d = X - f['mu'][k]; lp[:, k] = np.log(f['pi'][k]) - 0.5 * (np.einsum('ni,ij,nj->n', d, Si[k], d) + ld[k] + X.shape[1] * L2P)
    return lp
def mean_shift(X, f, iters=400, tol=1e-7):
    Si = np.linalg.inv(f['Sig']); Sim = np.einsum('kij,kj->ki', Si, f['mu']); X = X.copy(); sc = np.sqrt(np.trace(np.average(f['Sig'], axis=0, weights=f['pi'])) / f['mu'].shape[1] + 1e-12)
    act = np.ones(len(X), bool)
    for _ in range(iters):
        if not act.any(): break
        Xa = X[act]; lp = gmm_lp(Xa, f); w = np.exp(lp - logsumexp(lp, 1, keepdims=True)); A = np.einsum('nk,kij->nij', w, Si); bb = w @ Sim
        Xn = np.linalg.solve(A, bb[..., None])[..., 0]; ch = np.abs(Xn - Xa).max(1); X[act] = Xn; ii = np.where(act)[0]; act[ii[ch < tol * sc]] = False
    return X
def gmm_basins(f, seed=0):
    """modes of the mixture density and, per threshold, the surviving basins (mode -> basin) and sample masses; full diagram at inf"""
    r = np.random.default_rng(seed); K, D = f['mu'].shape; M = a.nsamp; comp = r.choice(K, M, p=f['pi'])
    Ls = np.linalg.cholesky(f['Sig']); Y = f['mu'][comp] + np.einsum('nij,nj->ni', Ls[comp], r.standard_normal((M, D)))
    st = np.vstack([f['mu'], Y[r.choice(M, min(M, 400), replace=False)]]); cv = mean_shift(st, f); lpc = logsumexp(gmm_lp(cv, f), 1)
    tot = np.cov(Y.T) + 1e-9 * np.eye(D); W = np.linalg.inv(np.linalg.cholesky(tot)); sc = np.sqrt(np.trace(tot) / D)
    modes = []
    for i in np.argsort(-lpc):
        if all(np.linalg.norm(cv[i] - m) > 1e-3 * sc * np.sqrt(D) for m in modes):
            pr = cv[i] + 0.02 * sc * r.standard_normal((24, D))
            if (logsumexp(gmm_lp(pr, f), 1) <= lpc[i] + 1e-9).all(): modes.append(cv[i])
    modes = np.array(modes) if modes else cv[[int(np.argmax(lpc))]]; nm = len(modes); Nn = np.vstack([Y, modes]); lpn = logsumexp(gmm_lp(Nn, f), 1); Wn = Nn @ W.T
    _, nb = cKDTree(Wn).query(Wn, k=13, workers=-1); u = np.repeat(np.arange(len(Nn)), 12); v = nb[:, 1:].ravel(); kk = np.unique(np.minimum(u, v) * len(Nn) + np.maximum(u, v))
    u, v = kk // len(Nn), kk % len(Nn); lev = np.minimum(lpn[u], lpn[v])
    for t in (0.25, 0.5, 0.75): lev = np.minimum(lev, logsumexp(gmm_lp((1 - t) * Nn[u] + t * Nn[v], f), 1))
    order = np.argsort(-lev); ul, vl, ll = u[order].tolist(), v[order].tolist(), lev[order].tolist(); ismode = np.zeros(len(Nn), bool); ismode[M:] = True
    def run(tau):
        par = list(range(len(Nn))); top = [i if ismode[i] else -1 for i in range(len(Nn))]; pers = []
        def fd(x):
            rt = x
            while par[rt] != rt: rt = par[rt]
            while par[x] != rt: par[x], x = rt, par[x]
            return rt
        for x, y, l_ in zip(ul, vl, ll):
            rx, ry = fd(x), fd(y)
            if rx == ry: continue
            tx, ty = top[rx], top[ry]
            if tx < 0 or ty < 0: par[rx if tx < 0 else ry] = ry if tx < 0 else rx; continue
            lo, hi = (rx, ry) if lpn[tx] < lpn[ty] else (ry, rx); prom = lpn[top[lo]] - l_
            if prom < tau: par[lo] = hi; pers.append(prom)
        tops = np.array([top[fd(i)] for i in range(len(Nn))]); return tops, sorted(pers, reverse=True)
    tinf, pinf = run(np.inf); res = {}
    for t in THR:
        tops, _ = run(t); surv = sorted(set(tops[M:].tolist()), key=lambda m: -(tops[:M] == m).sum())
        res[t] = dict(mode2b=np.array([surv.index(tops[M + j]) for j in range(nm)]), mass=np.array([(tops[:M] == m).mean() for m in surv]))
    return dict(modes=modes, lpm=lpn[M:], pinf=pinf, res=res, nm=nm)
def assign(post, f, B):
    """particle -> mode its posterior mean flows to -> surviving basin, per threshold"""
    ms = mean_shift(post, f); dm = ((ms[:, None, :] - B['modes'][None]) ** 2).sum(2); mi = np.argmin(dm, 1)
    return {t: B['res'][t]['mode2b'][mi] for t in THR}
def ari(x, y):
    _, xi = np.unique(x, return_inverse=True); _, yi = np.unique(y, return_inverse=True); C = np.zeros((xi.max() + 1, yi.max() + 1)); np.add.at(C, (xi, yi), 1)
    s = lambda v: (v * (v - 1) / 2.).sum(); n = len(x); t = s(C); ta = s(C.sum(1)); tb = s(C.sum(0)); e = ta * tb / (n * (n - 1) / 2.)
    return 1. if ta == tb == t else (t - e) / (0.5 * (ta + tb) - e)
def _fit(arg):
    X, S, Xi, K, seed, maxit = arg; f = xd_fit(X, S, K, Xi, maxit=maxit, seed=seed); B = gmm_basins(f, seed); return f, B, assign(f['post'], f, B)
def analyse(lab, X, S, Xi, eo, K, extra=None):
    """fit all particles and each half (parallel), basins at every threshold, half check; returns full-set labels per threshold"""
    with Pool(3) as p: out = p.map(_fit, [(X, S, Xi, K, a.seed, a.maxit)] + [(X[eo == h], S[eo == h], Xi[eo == h], K, a.seed + 1 + h, a.maxit) for h in (0, 1)])
    (f, B, Lf), H = out[0], out[1:]
    log('  %-9s n %6d K %2d: loglik/ptcl %.3f (%d it), weights %s; modes %d, prominences %s | halves modes %d / %d, prominences %s / %s' % (lab, len(X), K, f['ll'], f['it'],
        ' '.join('%.3f' % x for x in sorted(f['pi'], reverse=True)), B['nm'], ' '.join('%.2f' % x for x in B['pinf'][:5]) or '-', H[0][1]['nm'], H[1][1]['nm'],
        ' '.join('%.2f' % x for x in H[0][1]['pinf'][:4]) or '-', ' '.join('%.2f' % x for x in H[1][1]['pinf'][:4]) or '-'))
    for t in THR:
        lf = Lf[t]; nf = len(B['res'][t]['mass']); nh = [len(h[1]['res'][t]['mass']) for h in H]; ah = [ari(lf[eo == h], H[h][2][t]) for h in (0, 1)]
        rep = nh[0] == nf == nh[1] and min(ah) >= 0.8; ex = extra(t, lf) if extra else ''
        log('RESULT %-9s K %2d thr %.2f: %d basin(s) [halves %d / %d, ARI %.2f / %.2f, %s] masses %s (particles %s)%s' % (lab, K, t, nf, nh[0], nh[1], ah[0], ah[1],
            'REPRODUCED' if rep else 'not reproduced', ' '.join('%.3f' % x for x in B['res'][t]['mass'][:6]), ' '.join('%.3f' % x for x in (np.bincount(lf, minlength=nf) / len(lf))[:6]), ex))
    return Lf, f, B
def choose_k(X, S, Xi, eo, ksimple=None):
    if a.k == 'simple' and ksimple: return ksimple, 'SIMPLE'
    if a.k not in ('select', 'simple'): return int(a.k), 'fixed'
    kb, sc = select_k(X, S, Xi, eo, a.kmax, a.maxit); return kb, 'held-out (%s)' % ' '.join('%d:%.0f' % (k, sc[k] - max(sc.values())) for k in sorted(sc))
def bas1d(p, thr):
    K = len(p); G = -np.log(np.maximum(p, 1e-12)); mx = [i for i in range(K) if (i == 0 or G[i] <= G[i - 1]) and (i == K - 1 or G[i] <= G[i + 1])]
    B = [dict(lo=0, hi=K - 1, m=mx[0])]
    for j in range(1, len(mx)):
        sd = mx[j - 1] + int(np.argmax(G[mx[j - 1]:mx[j] + 1])); B[-1]['hi'] = sd; B.append(dict(lo=sd + 1, hi=K - 1, m=mx[j]))
    while len(B) > 1:
        bar = [G[B[k]['hi']] - max(G[B[k]['m']], G[B[k + 1]['m']]) for k in range(len(B) - 1)]; k = int(np.argmin(bar))
        if bar[k] >= thr: break
        m = B[k]['m'] if G[B[k]['m']] <= G[B[k + 1]['m']] else B[k + 1]['m']; B[k:k + 2] = [dict(lo=B[k]['lo'], hi=B[k + 1]['hi'], m=m)]
    return len(B)
def load_run(path):
    pind, z, eig, Q, ns, zd = read_emb(path); cf = os.path.join(os.path.dirname(path), 'flex_pca_coordinates.txt'); C = np.loadtxt(cf, comments='#')
    assert (C[:, 0].astype(int) == pind).all(), 'coordinates and embedding rows differ'
    pf = os.path.join(os.path.dirname(path), 'flex_pca_deconv_prior.txt'); ks, wts = None, ''
    if os.path.isfile(pf):
        h = open(pf).readline(); ks = int(h.split('K=')[1].split()[0]); ns = ns if ns else float(h.split('noise_scale=')[1].split()[0])
        wts = ' '.join('%.3f' % x for x in sorted(np.loadtxt(pf, comments='#', ndmin=2)[:, 1], reverse=True))
    x, S, fb = ml_form(z, Q, eig, ns if ns else 1.)
    log('%s: %d particles, D %d, noise scale %s, SIMPLE deconv K %s (weights %s); non-PD D_i %.3f; z_dec = coordinates file: %s' % (path, len(pind), z.shape[1], ns, ks, wts, fb,
        'yes' if zd is not None and np.allclose(zd, C[:, 5:], rtol=1e-5, atol=1e-6) else 'no'))
    return pind - 1, C[:, 1].astype(int), x, S, (zd if zd is not None else C[:, 5:]), ks

paths = a.emb.split(','); names = a.names.split(',') if a.names else ['run%d' % i for i in range(len(paths))]
if a.gt:
    rows, eo, X, S, Xi, ks = load_run(paths[0]); conf = np.loadtxt(a.gt, comments='#')[:, 1].astype(int)[rows]; keep = np.ones(len(rows), bool)
    if a.exclude: keep &= np.loadtxt(a.exclude).astype(int)[rows] == 0
    if a.discrete:
        ok = keep & (conf > 0); G = conf[ok]; log('discrete truth: %d non-junk; groups %s' % (ok.sum(), ' '.join('%d:%d' % (g, (G == g).sum()) for g in np.unique(G))))
        K, how = choose_k(X[ok], S[ok], Xi[ok], eo[ok], ks); log('K = %d by %s' % (K, how))
        ex = lambda t, lf: ' || ARI vs truth %.3f; basin: groups %s' % (ari(lf, G), ' ; '.join('b%d %s' % (b, ' '.join('%d=%d' % (g, ((lf == b) & (G == g)).sum()) for g in np.unique(G)
            if ((lf == b) & (G == g)).sum() > 0.02 * max(1, (lf == b).sum()))) for b in range(min(lf.max() + 1, 6))))
        analyse('nonjunk', X[ok], S[ok], Xi[ok], eo[ok], K, ex)
    else:
        i0, i1 = [int(x) for x in a.codes.split(':')]; KT = i1 - i0 + 1; tt = conf - i0; ok = keep & (tt >= 0) & (tt < KT); Nt = np.bincount(tt[ok], minlength=KT); tg = np.arange(KT, dtype=float)
        def target(name):
            if name == 'flat': g = np.zeros(KT)
            elif name.startswith('dw'): h = float(name[2:]); g = h * (((tg - 15.) / 10.) ** 2 - 1.) ** 2
            elif name == 'broad': g = 0.5 * ((tg - 15.) / 7.) ** 2
            elif name == 'shoulder':
                m_ = np.exp(-0.5 * ((tg - 20.) / 4.) ** 2); p = 0.95 * m_ / m_.sum() + 0.05 * (tg <= 3) / (tg <= 3).sum(); return p / p.sum()
            else: return None
            p = np.exp(-g); return p / p.sum()
        for s in a.shapes.split(','):
            cod = np.where(ok)[0]; p = target(s)
            if p is None: idx = cod
            else:
                sc = np.min(Nt / np.maximum(p, 1e-12)); nk = np.floor(sc * p).astype(int); r = np.random.default_rng(a.seed)
                idx = np.concatenate([r.choice(cod[tt[cod] == q], nk[q], replace=False) for q in range(KT)])
            ph = np.bincount(tt[idx], minlength=KT) / len(idx); K, how = choose_k(X[idx], S[idx], Xi[idx], eo[idx]); log('%s: K = %d by %s' % (s, K, how))
            ex = lambda t, lf, idx=idx, ph=ph: ' || truth (1-D) %d basin(s); basin mean t %s' % (bas1d(ph, t), ' '.join('%.1f' % tt[idx][lf == b].mean() for b in range(lf.max() + 1)))
            analyse(s, X[idx], S[idx], Xi[idx], eo[idx], K, ex)
else:
    R = {}
    for nm_, pth in zip(names, paths):
        rows, eo, X, S, Xi, ks = load_run(pth); keep = np.ones(len(rows), bool)
        if a.exclude: keep &= np.loadtxt(a.exclude).astype(int)[rows] == 0
        K, how = choose_k(X[keep], S[keep], Xi[keep], eo[keep], ks); log('%s: K = %d by %s' % (nm_, K, how))
        Lf, f, B = analyse(nm_, X[keep], S[keep], Xi[keep], eo[keep], K); R[nm_] = dict(rows=rows[keep], L=Lf)
        if a.labout and a.nm:
            for t in [float(x) for x in a.labthr.split(',')]:
                lab = np.zeros(a.nm, int); lab[rows[keep]] = Lf[t] + 1; np.savetxt('%s_%s_thr%.2f.txt' % (a.labout, nm_, t), lab, fmt='%d')
    common = np.array(sorted(set.intersection(*[set(r['rows'].tolist()) for r in R.values()]))); nn = list(R)
    for t in THR:
        M = {x: dict(zip(R[x]['rows'].tolist(), R[x]['L'][t].tolist())) for x in nn}; M = {x: np.array([M[x][c] for c in common]) for x in nn}
        mt = []
        for b in range(M[nn[0]].max() + 1):
            sel = M[nn[0]] == b; row_ = ['%s b%d %.3f' % (nn[0], b, sel.mean())]
            for y in nn[1:]:
                cnt = np.bincount(M[y][sel], minlength=M[y].max() + 1); bb = int(np.argmax(cnt)); row_.append('%s b%d %.3f (Jaccard %.2f)' % (y, bb, (M[y] == bb).mean(), cnt[bb] / float((sel | (M[y] == bb)).sum())))
            mt.append(', '.join(row_))
        log('CROSS thr %.2f: ARI %s | matched: %s' % (t, ' '.join('%s-%s %.3f' % (x, y, ari(M[x], M[y])) for i, x in enumerate(nn) for y in nn[i + 1:]), ' ;; '.join(mt[:6])))
log('DONE')
