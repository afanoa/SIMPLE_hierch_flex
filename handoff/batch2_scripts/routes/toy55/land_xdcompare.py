#!/usr/bin/env python3
"""land_xdcompare.py: is the DECONVOLVED latent landscape the same across independent runs? Per run: SIMPLE's noise model in ML form and
an extreme-deconvolution mixture at the run's own K (as land_xdbasins.py), fitted on all particles and on each half (eo). The runs are put
in the first run's frame by the least-squares affine map between their ML latents over the common particles (x_ref ~ M x_run + c), which
maps each mixture exactly (mu -> M mu + c, Sig -> M Sig M^T). Reports the Bhattacharyya coefficient BC = E_p[sqrt(q/p)] (Monte Carlo)
between runs' densities, with each run's even/odd BC as the noise ceiling; and the deconvolved marginals along the reference density's
principal axes (1/5/50/95/99 % points, in reference sd units) for every run and half: the extremes. Optional --fig PNG of the marginals.
usage: land_xdcompare.py --emb A,B,C --names a,b,c [--fig out.png --nmc 200000]"""
import os, argparse, numpy as np
from multiprocessing import Pool
from scipy.special import logsumexp
ap = argparse.ArgumentParser(); ap.add_argument('--emb', required=True); ap.add_argument('--names', required=True); ap.add_argument('--fig', default='')
ap.add_argument('--nmc', type=int, default=200000); ap.add_argument('--maxit', type=int, default=150); ap.add_argument('--seed', type=int, default=0)
ap.add_argument('--naxes', type=int, default=4); ap.add_argument('--save', default=''); a = ap.parse_args()
def log(*s): print(*s, flush=True)
L2P = np.log(2 * np.pi)
def read_emb(path):
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
    P = 1. / np.maximum(eig, 1e-300); w, V = np.linalg.eigh(Q - np.diag(P)[None]); fl = 1e-3 * P.min()
    Di = (V / np.maximum(w, fl)[:, None, :]) @ np.transpose(V, (0, 2, 1)); return (Di @ (Q @ z[..., None]))[..., 0], ns * Di
def kmeanspp(X, K, r):
    C = [X[r.integers(len(X))]]
    for _ in range(1, K):
        d2 = np.min(((X[:, None, :] - np.array(C)[None]) ** 2).sum(2), 1); C.append(X[r.choice(len(X), p=d2 / d2.sum())])
    C = np.array(C)
    for _ in range(20):
        lab = np.argmin(((X[:, None, :] - C[None]) ** 2).sum(2), 1); C = np.array([X[lab == k].mean(0) if (lab == k).any() else C[k] for k in range(K)])
    return C, lab
def xd_fit(arg):
    X, S, K, Xinit, maxit, seed = arg; tol = 1e-6; pimin = 1e-4
    n, D = X.shape; Xs = X  # SIMPLE's start (simple_xd_gmm init): equal-mass quantile means along the widest axis, Sigma = cov - mean noise
    jb = int(np.argmax(X.var(0))); grp = np.array_split(np.argsort(X[:, jb]), K); mu = np.array([X[g].mean(0) for g in grp]); pi = np.full(K, 1. / K)
    co = np.cov(X.T); w_, V_ = np.linalg.eigh(co - S.mean(0)); rid = 1e-3 * np.trace(co) / D; S0 = (V_ * np.maximum(w_, rid)) @ V_.T
    Sig = np.array([S0.copy() for _ in range(K)]); vtot = S0 + rid * np.eye(D); llo = -np.inf
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
    return dict(pi=pi, mu=mu, Sig=Sig, ll=tot / n, it=it + 1)
def lp_gmm(X, f):
    Si = np.linalg.inv(f['Sig']); ld = np.linalg.slogdet(f['Sig'])[1]
    return logsumexp(np.stack([np.log(f['pi'][k]) - 0.5 * (np.einsum('ni,ij,nj->n', X - f['mu'][k], Si[k], X - f['mu'][k]) + ld[k] + X.shape[1] * L2P) for k in range(len(f['pi']))], 1), 1)
def sample(f, m, r):
    c = r.choice(len(f['pi']), m, p=f['pi']); Ls = np.linalg.cholesky(f['Sig']); return f['mu'][c] + np.einsum('nij,nj->ni', Ls[c], r.standard_normal((m, f['mu'].shape[1])))
def bc(f, g, r):
    X = sample(f, a.nmc, r); return float(np.mean(np.exp(0.5 * (lp_gmm(X, g) - lp_gmm(X, f)))))
def amap(f, M, c): return dict(pi=f['pi'], mu=f['mu'] @ M.T + c, Sig=M[None] @ f['Sig'] @ M.T[None], ll=f['ll'], it=f['it'])
def marg_q(f, u, qs):
    """quantiles of the 1-D marginal of a mixture along unit vector u"""
    m = f['mu'] @ u; s = np.sqrt(np.einsum('i,kij,j->k', u, f['Sig'], u)); lo, hi = (m - 6 * s).min(), (m + 6 * s).max(); g = np.linspace(lo, hi, 4001)
    from scipy.stats import norm
    cdf = (f['pi'][None] * norm.cdf((g[:, None] - m[None]) / s[None])).sum(1); return np.interp(qs, cdf, g), g, (f['pi'][None] * norm.pdf((g[:, None] - m[None]) / s[None]) / s[None]).sum(1)

paths = a.emb.split(','); names = a.names.split(','); RUN = {}; jobs = []
for nm_, pth in zip(names, paths):
    pind, z, eig, Q, ns, zd = read_emb(pth); C = np.loadtxt(os.path.join(os.path.dirname(pth), 'flex_pca_coordinates.txt'), comments='#'); eo = C[:, 1].astype(int)
    hf = open(os.path.join(os.path.dirname(pth), 'flex_pca_deconv_prior.txt')).readline(); K = int(hf.split('K=')[1].split()[0]); ns = ns or float(hf.split('noise_scale=')[1].split()[0])
    X, S = ml_form(z, Q, eig, ns); RUN[nm_] = dict(pind=pind, X=X, eo=eo, K=K); log('%s: %d particles, K %d, noise scale %.3f' % (nm_, len(pind), K, ns))
    for part, sel in (('all', np.ones(len(X), bool)), ('even', eo == 0), ('odd', eo == 1)): jobs.append(((nm_, part), (X[sel], S[sel], K, zd[sel], a.maxit, a.seed)))
with Pool(min(9, len(jobs))) as p: fits = dict(zip([j[0] for j in jobs], p.map(xd_fit, [j[1] for j in jobs])))
for k_, f in fits.items(): log('fit %-6s %-4s: loglik/ptcl %.3f (%d it), weights %s' % (k_[0], k_[1], f['ll'], f['it'], ' '.join('%.3f' % x for x in sorted(f['pi'], reverse=True))))
ref = names[0]; Mp = {ref: (np.eye(RUN[ref]['X'].shape[1]), np.zeros(RUN[ref]['X'].shape[1]))}
for nm_ in names[1:]:
    common, ia, ib = np.intersect1d(RUN[ref]['pind'], RUN[nm_]['pind'], return_indices=True); Xa = RUN[ref]['X'][ia]; Xb = RUN[nm_]['X'][ib]
    B = np.linalg.lstsq(np.hstack([Xb, np.ones((len(Xb), 1))]), Xa, rcond=None)[0]; M, c = B[:-1].T, B[-1]; res = Xa - (Xb @ M.T + c)
    r2 = 1 - res.var(0) / Xa.var(0); sv = np.linalg.svd(M, compute_uv=False); Mp[nm_] = (M, c)
    log('map %s -> %s over %d common particles: R^2 per axis %s; singular values of M %s' % (nm_, ref, len(common), ' '.join('%.3f' % x for x in r2), ' '.join('%.2f' % x for x in sv)))
F = {k_: amap(f, *Mp[k_[0]]) for k_, f in fits.items()}; r = np.random.default_rng(a.seed + 7)
for nm_ in names: log('BC %-6s even vs odd (noise ceiling): %.4f' % (nm_, bc(F[(nm_, 'even')], F[(nm_, 'odd')], r)))
for i, x in enumerate(names):
    for y in names[i + 1:]: log('BC %s vs %s (all particles, %s frame): %.4f | halves even-even %.4f odd-odd %.4f' % (x, y, ref, bc(F[(x, 'all')], F[(y, 'all')], r),
        bc(F[(x, 'even')], F[(y, 'even')], r), bc(F[(x, 'odd')], F[(y, 'odd')], r)))
fr = F[(ref, 'all')]; Ctot = np.einsum('k,kij->ij', fr['pi'], fr['Sig'] + np.einsum('ki,kj->kij', fr['mu'], fr['mu'])) - np.outer(fr['pi'] @ fr['mu'], fr['pi'] @ fr['mu'])
ev, U = np.linalg.eigh(Ctot); U = U[:, ::-1]; ev = ev[::-1]; m0 = fr['pi'] @ fr['mu']; qs = np.array([0.01, 0.05, 0.5, 0.95, 0.99]); curves = {}
log('reference axes: variance fractions %s' % ' '.join('%.3f' % x for x in ev[:a.naxes] / ev.sum()))
for j in range(a.naxes):
    u = U[:, j]; sd = np.sqrt(ev[j])
    for k_ in sorted(F):
        q, g, pdf = marg_q(F[k_], u, qs); q = (q - m0 @ u) / sd; curves[(j,) + k_] = ((g - m0 @ u) / sd, pdf * sd)
        log('MARG axis %d %-6s %-4s: 1/5/50/95/99 %% points %s (sd units of the reference)' % (j + 1, k_[0], k_[1], ' '.join('%+.2f' % x for x in q)))
if a.save: np.savez(a.save, **{'%s_%s_%s' % (k_[0], k_[1], v): F[k_][v] for k_ in F for v in ('pi', 'mu', 'Sig')}, axes=U, ev=ev, m0=m0)
if a.fig:
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, a.naxes, figsize=(3.6 * a.naxes, 3.0)); col = dict(zip(names, ['#2a6fdb', '#d9822b', '#3a9d5d', '#8e5bd0']))
    for j in range(a.naxes):
        for nm_ in names:
            g, pdf = curves[(j, nm_, 'all')]; ax[j].plot(g, pdf, color=col[nm_], lw=1.6, label=nm_)
            for h, ls in (('even', ':'), ('odd', '--')): g, pdf = curves[(j, nm_, h)]; ax[j].plot(g, pdf, color=col[nm_], lw=0.8, ls=ls)
        ax[j].set_xlabel('axis %d (reference sd)' % (j + 1)); ax[j].set_yticks([]); ax[j].set_xlim(-4, 4)
    ax[0].set_ylabel('deconvolved density'); ax[0].legend(frameon=False, fontsize=8); fig.suptitle('r1 deconvolved latent marginals: runs (solid), even (dotted) / odd (dashed) halves', fontsize=9)
    fig.tight_layout(); fig.savefig(a.fig, dpi=150); log('figure %s' % a.fig)
log('DONE')
