#!/usr/bin/env python3
"""States as modes of the DECONVOLVED latent density (RECOVAR's reading, instead of lattice cells or k-means). Every
particle's posterior over the board has a mean z_i and a covariance S_i; the density of the means is the conformational
density blurred by the S_i. Two deconvolutions, both read for modes:
 (1) extreme deconvolution (Bovy, Hogg, Roweis 2011): a Gaussian mixture p(z) = sum_k pi_k N(mu_k, V_k) fitted by EM to
     x_i ~ N(z, S_i) from K0 components; its modes found by mean-shift on the mixture density (components climbing to the
     same peak merge); particles assigned by the noise-aware posterior p(k | x_i) ~ pi_k N(x_i; mu_k, V_k + S_i).
 (2) Richardson-Lucy deconvolution of the 2-D histogram of the means by the mean posterior kernel, then modes persistent
     across bandwidths (scale space); particles assigned by hill-climbing from their mean.
Readouts vs truth: modes found, one-mode-per-conformer purity (Hungarian) and recall, against k-means K on the means,
the argmax cells and the same mixture WITHOUT deconvolution (S_i = 0).
usage: density_modes.py KLEM_DIR [--gt gt_particles.txt] [--k0 20] [--kmeans 10] [--nfit 100000]"""
import sys, os, argparse, numpy as np
from scipy.optimize import linear_sum_assignment
from scipy.ndimage import gaussian_filter, maximum_filter
ap = argparse.ArgumentParser(); ap.add_argument('klem'); ap.add_argument('--gt', default=''); ap.add_argument('--k0', type=int, default=20); ap.add_argument('--kmeans', type=int, default=10)
ap.add_argument('--nfit', type=int, default=100000); ap.add_argument('--iters', type=int, default=150); ap.add_argument('--prune', type=float, default=0.005); a = ap.parse_args()
rng = np.random.default_rng(0); Z = np.load(a.klem + '/klem_resp.npz'); r = Z['r'].astype(float); pinds = Z['pinds']; i1 = Z['i1'].astype(float); i2 = Z['i2'].astype(float); ngrid = len(i1)
rg = r[:, :ngrid]; mass = rg.sum(1); on = mass > 0.5; rn = rg[on] / mass[on][:, None]; C = np.c_[i1, i2]; x = rn @ C
M2 = rn @ np.einsum('mi,mj->mij', C, C).reshape(ngrid, 4); S = M2.reshape(-1, 2, 2) - np.einsum('ni,nj->nij', x, x); S += 1e-3 * np.eye(2)[None]
am = rg[on].argmax(1); pind = pinds[on]; n = len(x)
gt = {}
if a.gt:
    for l in open(a.gt):
        if l.strip() and not l.startswith('#'): t = l.split(); gt[int(t[0])] = int(t[1])
truth = np.array([gt.get(int(p), 0) for p in pind]); good = truth > 0
print('%s: %d particles on the grid (%d with a true conformer); posterior sd per particle median %.2f cells (5-95 pct %.2f-%.2f)' % (a.klem, n, good.sum(), np.median(np.sqrt(np.trace(S, axis1=1, axis2=2) / 2)), *np.percentile(np.sqrt(np.trace(S, axis1=1, axis2=2) / 2), [5, 95])))


def score(lab, tag):
    """one-state-per-conformer purity (Hungarian), recall at that state, states used"""
    L_ = lab[good]; T_ = truth[good]; us, ui = np.unique(L_, return_inverse=True); cs = np.unique(T_); Cm = np.zeros((len(us), len(cs)))
    for j, c in enumerate(cs): Cm[:, j] = np.bincount(ui[T_ == c], minlength=len(us))
    P = Cm / np.maximum(Cm.sum(1, keepdims=True), 1); R = Cm / np.maximum(Cm.sum(0, keepdims=True), 1); ri, ci = linear_sum_assignment(-P); pur = np.zeros(len(cs)); rec = np.zeros(len(cs)); pur[ci] = P[ri, ci]; rec[ci] = R[ri, ci]
    owned = len(set(P.argmax(1)))
    print('  %-44s %2d states | purity mean %.2f (min %.2f) recall mean %.2f | conformers owning a state %d/%d | %s' % (tag, len(us), pur.mean(), pur.min(), rec.mean(), owned, len(cs), ' '.join('%d:%.2f' % (c, pur[j]) for j, c in enumerate(cs))))
    return pur.mean()


def inv2(A):
    d = A[..., 0, 0] * A[..., 1, 1] - A[..., 0, 1] * A[..., 1, 0]; B = np.empty_like(A); B[..., 0, 0] = A[..., 1, 1]; B[..., 1, 1] = A[..., 0, 0]; B[..., 0, 1] = -A[..., 0, 1]; B[..., 1, 0] = -A[..., 1, 0]
    return B / np.maximum(d, 1e-30)[..., None, None], d


def kmeans(X, K, iters=100, seeds=5):
    best = None
    for s in range(seeds):
        Cc = X[rng.choice(len(X), K, replace=False)]
        for _ in range(iters):
            d = ((X[:, None, :] - Cc[None]) ** 2).sum(2); lab = d.argmin(1); Cn = np.array([X[lab == k].mean(0) if (lab == k).any() else Cc[k] for k in range(K)])
            if np.allclose(Cn, Cc): break
            Cc = Cn
        cost = d[np.arange(len(X)), lab].sum()
        if best is None or cost < best[0]: best = (cost, lab, Cc)
    return best[1], best[2]


def xd_fit(X, Sx, K, iters, deconv=True):
    """extreme deconvolution EM: mixture of Gaussians for the TRUE positions given noisy x_i with covariances S_i"""
    lab, mu = kmeans(X, K); pi = np.bincount(lab, minlength=K) / len(lab) + 1e-6; pi /= pi.sum(); V = np.array([np.cov(X[lab == k].T) + 0.05 * np.eye(2) if (lab == k).sum() > 2 else 0.5 * np.eye(2) for k in range(K)])
    Sz = Sx if deconv else np.zeros_like(Sx)
    for it in range(iters):
        Tm = V[None] + Sz[:, None]; Ti, det = inv2(Tm); dx = X[:, None, :] - mu[None]; q = np.einsum('nki,nkij,nkj->nk', dx, Ti, dx)
        logp = np.log(pi)[None] - 0.5 * q - 0.5 * np.log(np.maximum(det, 1e-30)); logp -= logp.max(1, keepdims=True); q = np.exp(logp); q /= q.sum(1, keepdims=True)
        Nk = q.sum(0) + 1e-12; VTi = np.einsum('kij,nkjl->nkil', V, Ti); b = mu[None] + np.einsum('nkij,nkj->nki', VTi, dx); Bm = V[None] - np.einsum('nkij,kjl->nkil', VTi, V)
        mu_n = np.einsum('nk,nki->ki', q, b) / Nk[:, None]; db = mu_n[None] - b; V = (np.einsum('nk,nki,nkj->kij', q, db, db) + np.einsum('nk,nkij->kij', q, Bm)) / Nk[:, None, None] + 1e-4 * np.eye(2)[None]; mu = mu_n; pi = Nk / Nk.sum()
        keep = pi > a.prune
        if keep.sum() < len(pi): pi, mu, V = pi[keep] / pi[keep].sum(), mu[keep], V[keep]
    return pi, mu, V


def mixture_modes(pi, mu, V):
    """modes of a Gaussian mixture by fixed-point mean-shift started at every component mean; returns mode index per component and the modes"""
    Vi, det = inv2(V); starts = mu.copy(); ends = []
    for z in starts:
        for _ in range(200):
            dx = z[None] - mu; q = np.einsum('ki,kij,kj->k', dx, Vi, dx); wk = pi * np.exp(-0.5 * q) / np.sqrt(det); A = np.einsum('k,kij->ij', wk, Vi); bvec = np.einsum('k,kij,kj->i', wk, Vi, mu); zn = np.linalg.solve(A + 1e-12 * np.eye(2), bvec)
            if np.linalg.norm(zn - z) < 1e-5: z = zn; break
            z = zn
        ends.append(z)
    ends = np.array(ends); modes = []; mode_of = np.zeros(len(mu), int)
    for k, e in enumerate(ends):
        for j, m in enumerate(modes):
            if np.linalg.norm(e - m) < 0.15: mode_of[k] = j; break
        else: modes.append(e); mode_of[k] = len(modes) - 1
    return mode_of, np.array(modes)


def xd_assign(X, Sx, pi, mu, V, mode_of, deconv=True):
    Sz = Sx if deconv else np.zeros_like(Sx); Tm = V[None] + Sz[:, None]; Ti, det = inv2(Tm); dx = X[:, None, :] - mu[None]; q = np.einsum('nki,nkij,nkj->nk', dx, Ti, dx)
    logp = np.log(pi)[None] - 0.5 * q - 0.5 * np.log(np.maximum(det, 1e-30)); logp -= logp.max(1, keepdims=True); q = np.exp(logp); q /= q.sum(1, keepdims=True)
    nm = mode_of.max() + 1; pm = np.zeros((len(X), nm))
    for k in range(len(pi)): pm[:, mode_of[k]] += q[:, k]
    return pm


# ---- baselines
score(am, 'argmax cells')
lab_km, _ = kmeans(x[rng.choice(n, min(n, a.nfit), replace=False)], a.kmeans); _, Ck = kmeans(x[rng.choice(n, min(n, 50000), replace=False)], a.kmeans)
score(((x[:, None, :] - Ck[None]) ** 2).sum(2).argmin(1), 'k-means %d on the posterior means' % a.kmeans)
# ---- (1) extreme deconvolution, with and without the per-particle covariances
fit = rng.choice(n, min(n, a.nfit), replace=False)
for deconv in (False, True):
    pi, mu, V = xd_fit(x[fit], S[fit], a.k0, a.iters, deconv=deconv); mode_of, modes = mixture_modes(pi, mu, V); pm = xd_assign(x, S, pi, mu, V, mode_of, deconv=deconv); lab = pm.argmax(1)
    pop = np.bincount(lab, minlength=len(modes)) / n
    tag = 'mixture K0=%d %s: %d comps -> %d modes' % (a.k0, 'DECONVOLVED (S_i)' if deconv else 'plain (S_i = 0)', len(pi), len(modes))
    score(lab, tag)
    if deconv:
        dom = []
        for j in range(len(modes)):
            m = (lab == j) & good; dom.append('%d(%s %.2f, n %.3f)' % (j + 1, np.bincount(truth[m]).argmax() if m.any() else '-', (np.bincount(truth[m]).max() / m.sum()) if m.any() else 0., pop[j]))
        print('     modes at (z1,z2): %s' % ' '.join('%d:(%.1f,%.1f)' % (j + 1, *modes[j]) for j in range(len(modes))))
        print('     mode: dominant conformer, its share, mode mass: %s' % ' '.join(dom))
        print('     mean max mode posterior %.2f (particles beyond 0.9: %.2f)' % (pm.max(1).mean(), np.mean(pm.max(1) > 0.9)))
# ---- (2) Richardson-Lucy deconvolution of the histogram of the means by the mean posterior kernel + scale-space modes
n1 = int(round(i1.max())) + 1; n2 = int(round(i2.max())) + 1; res = 10; e1 = np.linspace(-0.5, n1 - 0.5, res * n1 + 1); e2 = np.linspace(-0.5, n2 - 0.5, res * n2 + 1)
Hh, _, _ = np.histogram2d(x[:, 0], x[:, 1], bins=[e1, e2]); Sbar = S.mean(0); sd = np.sqrt(np.diag(Sbar)) * res   # mean posterior sd in histogram pixels per axis
blur = lambda v: gaussian_filter(v, sd, mode='constant'); est = np.full_like(Hh, Hh.mean()); obs = Hh + 1e-9
for _ in range(40): est *= blur(obs / np.maximum(blur(est), 1e-9))
found = {}
for h in (0.25, 0.35, 0.5, 0.7, 1.0):
    sm = gaussian_filter(est, h * res, mode='constant'); pk = (sm == maximum_filter(sm, size=int(2 * h * res) + 1)) & (sm > 0.02 * sm.max()); ij = np.argwhere(pk); found[h] = ij
persist = []   # modes at the finest bandwidth that have a partner within h at the next two bandwidths
hs = sorted(found)
for ij in found[hs[0]]:
    okp = all(np.min(np.linalg.norm(found[h2] - ij, axis=1), initial=1e9) <= 1.5 * h2 * res for h2 in hs[1:3]) if all(len(found[h2]) for h2 in hs[1:3]) else False
    if okp: persist.append(ij)
persist = np.array(persist) if persist else np.zeros((0, 2), int)
print('  RL deconvolution (mean posterior sd %.2f, %.2f cells): modes per bandwidth %s; persistent modes at the finest: %d' % (np.sqrt(Sbar[0, 0]), np.sqrt(Sbar[1, 1]), ' '.join('%.2f:%d' % (h, len(found[h])) for h in hs), len(persist)))
if len(persist):
    sm = gaussian_filter(est, hs[0] * res, mode='constant'); gi = np.clip(np.rint((x[:, 0] + 0.5) * res).astype(int), 0, sm.shape[0] - 1); gj = np.clip(np.rint((x[:, 1] + 0.5) * res).astype(int), 0, sm.shape[1] - 1)
    for _ in range(3 * res):   # discrete hill-climb on the smoothed deconvolved density
        best_i, best_j = gi.copy(), gj.copy(); bv = sm[gi, gj]
        for di in (-1, 0, 1):
            for dj in (-1, 0, 1):
                ii = np.clip(gi + di, 0, sm.shape[0] - 1); jj = np.clip(gj + dj, 0, sm.shape[1] - 1); v = sm[ii, jj]; mv = v > bv; best_i[mv], best_j[mv], bv[mv] = ii[mv], jj[mv], v[mv]
        if np.all((best_i == gi) & (best_j == gj)): break
        gi, gj = best_i, best_j
    end = np.c_[gi, gj]; lab = np.linalg.norm(end[:, None, :] - persist[None], axis=2).argmin(1)
    score(lab, 'RL-deconvolved density, %d persistent modes' % len(persist))
    print('     modes at (z1,z2): %s' % ' '.join('%d:(%.1f,%.1f)' % (j + 1, persist[j][0] / res - 0.5, persist[j][1] / res - 0.5) for j in range(len(persist))))
