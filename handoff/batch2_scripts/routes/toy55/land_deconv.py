#!/usr/bin/env python3
"""land_deconv.py (module): smooth deconvolution of an occupancy landscape on a grid from per-particle likelihoods L (N x K, rows scaled
freely). Maximises mean_i log(L_i . p) - lam ||D2 G||^2 over G = -ln p (p = softmax(-G); the second-difference penalty leaves a constant and
a linear tilt of G free), lam chosen by CROSS-HALF held-out likelihood (fit on one half, score the other, both ways), bootstrap bands.
npmle() = the unpenalised grid NPMLE (EM + SQUAREM)."""
import numpy as np
from scipy.optimize import minimize
LAMS = (0., 1e-6, 3e-6, 1e-5, 3e-5, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 0.1)
def _pk(G): e = np.exp(-(G - G.min())); return e / e.sum()
def fit_pen(L, lam, G0=None, maxiter=400):
    N, K = L.shape; D2 = np.diff(np.eye(K), 2, axis=0); P = D2.T @ D2
    def f(G):
        p = _pk(G); d = np.maximum(L @ p, 1e-300); r = (L * p[None]) / d[:, None]
        return -np.log(d).mean() + lam * G @ P @ G, (r.mean(0) - p) + 2. * lam * (P @ G)
    x0 = np.zeros(K) if G0 is None else G0 - G0.mean()
    r = minimize(f, x0, jac=True, method='L-BFGS-B', options=dict(maxiter=maxiter, gtol=1e-10, ftol=1e-14))
    G = r.x - r.x.min(); return G, _pk(G)
def heldout(L, p): return float(np.log(np.maximum(L @ p, 1e-300)).mean())
def cv_lambda(L, half, lams=LAMS):
    """half: 0/1 per particle. Returns (best lam, table of summed held-out log-lik per lam)."""
    tab = []
    for lam in lams:
        s = 0.
        for h in (0, 1):
            _, p = fit_pen(L[half == h], lam); s += heldout(L[half != h], p) * (half != h).sum()
        tab.append(s / len(half))
    return float(lams[int(np.argmax(tab))]), np.array(tab)
def fit_cv_boot(L, half, nboot=20, seed=0, lams=LAMS):
    lam, tab = cv_lambda(L, half, lams); G, p = fit_pen(L, lam); rng = np.random.default_rng(seed); B = []
    for b in range(nboot):
        ib = rng.integers(0, len(L), len(L)); B.append(fit_pen(L[ib], lam, G0=G)[1])
    return dict(lam=lam, cvtab=tab, p=p, G=G, boot=np.array(B))
def _F(L, p): d = L @ p; return p * (L.T @ (1. / d)) / L.shape[0]
def _obj(L, p): return np.log(L @ p).mean()
def npmle(L, p0=None, nit=3000, tol=1e-10):
    p = np.full(L.shape[1], 1. / L.shape[1]) if p0 is None else p0.copy()
    for _ in range(nit):
        p1 = _F(L, p); p2 = _F(L, p1); r = p1 - p; v = p2 - 2. * p1 + p; nv = np.linalg.norm(v)
        if nv < 1e-14: p = p2; break
        al = min(-np.linalg.norm(r) / nv, -1.); pq = np.maximum(p - 2. * al * r + al * al * v, 0.); pq /= pq.sum(); pq = _F(L, pq)
        if not np.isfinite(_obj(L, pq)) or _obj(L, pq) < _obj(L, p2): pq = p2
        done = np.abs(pq - p).max() < tol; p = pq
        if done: break
    return p
def sm3(p): return np.convolve(p, np.ones(3), 'same') / np.convolve(np.ones_like(p), np.ones(3), 'same')
def metrics(p, pt, wells=(5, 25), top=15, minor=4):
    G = -np.log(np.maximum(sm3(p), 1e-12)); return dict(tv=0.5 * float(np.abs(p - pt).sum()), low=float(p[:top].sum() + 0.5 * p[top]), minor=float(p[:minor].sum()),
                                                     barrier=float(G[top] - min(G[wells[0]], G[wells[1]])))
def band(boot, pt, q=(10, 50, 90)):
    return {k: [float(np.percentile([metrics(b, pt)[k] for b in boot], qq)) for qq in q] for k in ('tv', 'low', 'minor', 'barrier')}
