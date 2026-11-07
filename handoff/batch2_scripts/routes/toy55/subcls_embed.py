#!/usr/bin/env python3
"""Sub-class embedding of one tree node (prototype, truth readouts optional). The node's 2D sub-classes (cls_expansion inside
its direction bins) are pooling bags; every bag is placed in the node's 3D latent frame = the flex_pca coordinates of its members,
each computed at the member's own pose and CTF, pooled over all / even / odd members (no image is averaged across poses and no
per-view embedding is synchronised: the 3D basis is the common frame).
Truth-free readouts: (1) split-half reliability of the pooled coordinates (maximum-reliability directions, generalised eigenproblem
of the cross-half vs within-half covariance over sub-classes) against a null that shuffles members between the sub-classes of the
same bin -> number of reproducible dimensions d (0 = no structure: stop); (2) extreme-deconvolution Gaussian mixture of the
sub-class means in those d dimensions (per sub-class noise = within-sub-class member covariance / members; K by held-out cross-half
log-likelihood) -> components chained when their Ashman separation D < 2 -> groups; a group with reproducible internal spread
(deconvolved sd > noise floor) is a continuum, without it a discrete state; (3) junk = sub-classes whose members' latent-model
residual is a robust outlier; (4) children = groups, or equal-width bins of a continuum's leading coordinate over sub-classes.
Truth (gt_particles.txt; conformer 200..230 = swing 0..30 A, 0 = junk): the sub-class partition's swing eta2 (ceiling of any
sub-class labelling), per-particle vs pooled latent (held-out ridge R2 of the swing), children eta2, junk recall, conformer NMI.
usage: subcls_embed.py --labels cls_expansion_labels.txt --flex flex_pca_coordinates.txt --rows rows_node.txt [--gt gt] --out DIR"""
import argparse, os, numpy as np
from scipy.linalg import eigh
from scipy.sparse.csgraph import connected_components
ap = argparse.ArgumentParser()
ap.add_argument('--labels', required=True); ap.add_argument('--flex', required=True); ap.add_argument('--rows', required=True)
ap.add_argument('--gt', default=''); ap.add_argument('--out', required=True); ap.add_argument('--minh', type=int, default=5)
ap.add_argument('--nchild', type=int, default=10); ap.add_argument('--nnull', type=int, default=20); ap.add_argument('--kmax', type=int, default=8)
ap.add_argument('--refs', default='', help='comma list name=per-row label file (particle-route labels scored on the same members)')
ap.add_argument('--seed', type=int, default=0)
a = ap.parse_args(); os.makedirs(a.out, exist_ok=True); rng = np.random.default_rng(a.seed)

# ---------------- inputs: members with a sub-class and a latent
L = np.loadtxt(a.labels, comments='#', ndmin=2).astype(int)                     # particle(1-based row) parent_bin global_subclass
F = np.loadtxt(a.flex, comments='#', ndmin=2)                                    # particle eo label residual mean_residual z1..zK
rows = np.loadtxt(a.rows).astype(int); NR = len(rows)
sub_of = np.zeros(NR + 1, int); bin_of = np.zeros(NR + 1, int); sub_of[L[:, 0]] = L[:, 2]; bin_of[L[:, 0]] = L[:, 1]
fi = F[:, 0].astype(int); keep = (rows[fi - 1] > 0) & (sub_of[fi] > 0)
p = fi[keep]; eo = F[keep, 1].astype(int); res = F[keep, 3]; mres = F[keep, 4]; Z = F[keep, 5:]
sub = sub_of[p]; binn = bin_of[p]; N, nz = Z.shape
us, sid = np.unique(sub, return_inverse=True); S = len(us); nS = np.bincount(sid, minlength=S)
sbin = np.zeros(S, int); sbin[sid] = binn
print('node: %d members with a sub-class and a %d-dim latent; %d sub-classes in %d bins; members per sub-class p10/50/90 %d/%d/%d (min %d, max %d)'
      % (N, nz, S, len(np.unique(binn)), *np.percentile(nS, [10, 50, 90]).astype(int), nS.min(), nS.max()))

def gmeans(X, g, G, w=None):
    w = np.ones(len(g)) if w is None else w
    c = np.bincount(g, w, minlength=G); M = np.stack([np.bincount(g, w * X[:, j], minlength=G) for j in range(X.shape[1])], 1)
    return M / np.maximum(c, 1e-12)[:, None], c

def pooled(Zx, g):  # all / even / odd member means and counts per sub-class
    Ma, na = gmeans(Zx, g, S); Me, ne = gmeans(Zx, g, S, (eo == 0).astype(float)); Mo, no = gmeans(Zx, g, S, (eo == 1).astype(float))
    return Ma, Me, Mo, na, ne, no

def reliab(Me, Mo, w):  # maximum split-half reliability directions: Cx v = lam Cw v over sub-classes (weights w)
    W = w / w.sum(); mu = (W[:, None] * (Me + Mo) / 2).sum(0); A = Me - mu; B = Mo - mu
    Cx = (A * W[:, None]).T @ B; Cx = (Cx + Cx.T) / 2; Cw = ((A * W[:, None]).T @ A + (B * W[:, None]).T @ B) / 2
    lam, V = eigh(Cx, Cw + 1e-9 * np.trace(Cw) / len(Cw) * np.eye(len(Cw))); o = np.argsort(lam)[::-1]
    return lam[o], V[:, o]

def binshuffle(g):  # members permuted between the sub-classes of their bin (sizes and eo kept)
    out = g.copy()
    for b in np.unique(binn):
        m = np.where(binn == b)[0]; out[m] = g[m][rng.permutation(len(m))]
    return out

def eta2(y, lab):  # correlation ratio squared (raw, size-adjusted epsilon^2)
    u, inv = np.unique(lab, return_inverse=True); K = len(u); n = len(y)
    m = np.bincount(inv, y) / np.bincount(inv); W = ((y - m[inv]) ** 2).sum(); T = ((y - y.mean()) ** 2).sum()
    return 1 - W / T, 1 - (W / max(n - K, 1)) / (T / (n - 1))

def cv_ridge(X, y, groups, nfold=5, alpha=1e-3):  # held-out R2 and corr of a ridge fit; folds over groups (no member leakage)
    ug = np.unique(groups); f = rng.permutation(len(ug)) % nfold; fold = f[np.searchsorted(ug, groups)]; pred = np.zeros(len(y))
    for k in range(nfold):
        tr, te = fold != k, fold == k; mu = X[tr].mean(0); A = X[tr] - mu; ym = y[tr].mean()
        lamI = alpha * np.trace(A.T @ A) / A.shape[1] * np.eye(A.shape[1])
        pred[te] = (X[te] - mu) @ np.linalg.solve(A.T @ A + lamI, A.T @ (y[tr] - ym)) + ym
    return 1 - ((y - pred) ** 2).sum() / ((y - y.mean()) ** 2).sum(), np.corrcoef(pred, y)[0, 1]

def xd_fit(Y, Sn, K, nit=200, init=None):  # extreme deconvolution (Bovy 2011): Y_a ~ sum_k pi_k N(mu_k, V_k + Sn_a)
    n, d = Y.shape
    if init is None:
        o = np.argsort(Y[:, 0]); mu = np.stack([Y[o[int((k + 0.5) * n / K)]] for k in range(K)])
    else: mu = init.copy()
    V = np.stack([np.cov(Y.T).reshape(d, d) / K ** 2 + 1e-6 * np.eye(d) for _ in range(K)]); pi = np.ones(K) / K; ll_old = -np.inf
    for it in range(nit):
        T = V[None] + Sn[:, None]; Ti = np.linalg.inv(T); r = Y[:, None] - mu[None]
        lg = -0.5 * (np.einsum('aki,akij,akj->ak', r, Ti, r) + np.linalg.slogdet(T)[1] + d * np.log(2 * np.pi)) + np.log(pi + 1e-300)
        mx = lg.max(1, keepdims=True); R = np.exp(lg - mx); sR = R.sum(1, keepdims=True); ll = (mx[:, 0] + np.log(sR[:, 0])).sum(); R /= sR
        b = mu[None] + np.einsum('kij,akjl,akl->aki', V, Ti, r); B = V[None] - np.einsum('kij,akjl,klm->akim', V, Ti, V)
        rk = R.sum(0) + 1e-12; pi = rk / n; mu = np.einsum('ak,aki->ki', R, b) / rk[:, None]
        dm = mu[None] - b; V = (np.einsum('ak,aki,akj->kij', R, dm, dm) + np.einsum('ak,akij->kij', R, B)) / rk[:, None, None] + 1e-9 * np.eye(d)
        if ll - ll_old < 1e-6 * abs(ll): break
        ll_old = ll
    return pi, mu, V, ll

def xd_ll(Y, Sn, pi, mu, V):
    d = Y.shape[1]; T = V[None] + Sn[:, None]; Ti = np.linalg.inv(T); r = Y[:, None] - mu[None]
    lg = -0.5 * (np.einsum('aki,akij,akj->ak', r, Ti, r) + np.linalg.slogdet(T)[1] + d * np.log(2 * np.pi)) + np.log(pi + 1e-300)
    mx = lg.max(1, keepdims=True); return (mx[:, 0] + np.log(np.exp(lg - mx).sum(1))).sum(), lg

# ---------------- truth
GT = None
if a.gt:
    G = np.loadtxt(a.gt, comments='#', ndmin=2).astype(float); conf_of = np.zeros(NR + 1, int); conf_of[G[:, 0].astype(int)] = G[:, 1].astype(int)
    conf = conf_of[p]; isw = (conf >= 200) & (conf <= 230); swing = np.where(isw, conf - 200, np.nan); junk = conf == 0
    GT = dict(conf=conf, isw=isw, swing=swing, junk=junk)
    uc, cc = np.unique(conf, return_counts=True); o = np.argsort(-cc)
    print('TRUTH node composition: %d swing (%.1f %%), %d junk (%.1f %%), other conformers %s' % (isw.sum(), 100 * isw.mean(), junk.sum(), 100 * junk.mean(),
          ' '.join('%d:%d' % (uc[i], cc[i]) for i in o if uc[i] > 0 and not (200 <= uc[i] <= 230))))
    e_raw, e_adj = eta2(swing[isw], sid[isw]); eb = eta2(swing[isw], binn[isw])
    sw_sub = np.bincount(sid, isw, minlength=S) / nS; jk_sub = np.bincount(sid, junk, minlength=S) / nS
    print('TRUTH sub-class partition: swing eta2 %.3f (size-adjusted %.3f; the bins alone %.3f / %.3f) = ceiling of any sub-class labelling (corr %.3f)'
          % (e_raw, e_adj, *eb, np.sqrt(max(e_adj, 0))))
    print('TRUTH sub-class swing share p10/50/90 %.2f/%.2f/%.2f; junk: %d particles, %.0f %% of them in sub-classes >= 50 %% junk, %d such sub-classes'
          % (*np.percentile(sw_sub, [10, 50, 90]), junk.sum(), 100 * (jk_sub[sid][junk] >= 0.5).mean() if junk.any() else 0, (jk_sub >= 0.5).sum()))

# ---------------- embeddings: raw latent and bin-centred latent
bm, _ = gmeans(Z, binn - binn.min(), binn.max() - binn.min() + 1); Zc = Z - bm[binn - binn.min()]
summary = {}
for tag, Zx in (('raw', Z), ('bincentred', Zc)):
    Ma, Me, Mo, na, ne, no = pooled(Zx, sid); ok = (ne >= a.minh) & (no >= a.minh); w = ne * no / np.maximum(ne + no, 1)
    lam, V = reliab(Me[ok], Mo[ok], w[ok]); null1 = []
    for t in range(a.nnull):
        g2 = binshuffle(sid); _, Me2, Mo2, _, ne2, no2 = pooled(Zx, g2); ok2 = (ne2 >= a.minh) & (no2 >= a.minh); w2 = ne2 * no2 / np.maximum(ne2 + no2, 1)
        null1.append(reliab(Me2[ok2], Mo2[ok2], w2[ok2])[0][0])
    thr = np.percentile(null1, 95); d = int((lam > thr).sum())
    print('\n[%s] split-half reliability over %d sub-classes (>= %d members per half): %s | null (within-bin member shuffle, %d draws) top %.3f (95th pct %.3f) -> %d reproducible dimension(s)'
          % (tag, ok.sum(), a.minh, ' '.join('%.3f' % l for l in lam[:6]), a.nnull, np.max(null1), thr, d))
    mu0 = Zx.mean(0); Ys = (Ma - mu0) @ V; Ye = (Me - mu0) @ V; Yo = (Mo - mu0) @ V; Yp = (Zx - mu0) @ V   # sub-class / half / particle coordinates
    if GT is not None:
        sw = GT['isw']
        for j in range(min(3, nz)):
            print('  TRUTH direction %d (reliability %.3f): corr with swing: sub-class means %.3f (members inherit) | per particle %.3f'
                  % (j + 1, lam[j], np.corrcoef(Ys[sid[sw], j], GT['swing'][sw])[0, 1], np.corrcoef(Yp[sw, j], GT['swing'][sw])[0, 1]))
        r2p, cp = cv_ridge(Zx[sw], GT['swing'][sw], sid[sw]); r2s, cs = cv_ridge(Ma[sid[sw]], GT['swing'][sw], sid[sw])
        print('  TRUTH held-out ridge of the swing on the full %d-dim latent: per particle R2 %.3f (corr %.3f) | pooled over sub-classes R2 %.3f (corr %.3f)' % (nz, r2p, cp, r2s, cs))
    if d == 0:
        print('  STOP: no reproducible sub-class structure in the node latent'); summary[tag] = dict(d=0); continue
    # noise model in the reproducible subspace: within-sub-class member covariance from the half differences
    dd = max(d, 1); Dh = (Ye[:, :dd] - Yo[:, :dd])[ok]; f = (1. / ne[ok] + 1. / no[ok]); Sig = (Dh.T @ Dh) / f.sum()   # E[(ye-yo)(ye-yo)^T] = Sig (1/ne + 1/no)
    Y = Ys[ok, :dd]; Sn = Sig[None] / na[ok][:, None, None]; Sne = Sig[None] / ne[ok][:, None, None]; Sno = Sig[None] / no[ok][:, None, None]
    best = None; rows_k = []
    for K in range(1, a.kmax + 1):
        pe = xd_fit(Ye[ok, :dd], Sne, K); po = xd_fit(Yo[ok, :dd], Sno, K)
        ho = 0.5 * (xd_ll(Yo[ok, :dd], Sno, *pe[:3])[0] + xd_ll(Ye[ok, :dd], Sne, *po[:3])[0]); rows_k.append(ho)
        if best is None or ho > best[1] + 2.0: best = (K, ho)                   # a component must buy > 2 nats of held-out likelihood
    K = best[0]; pi, mu, Vk, ll = xd_fit(Y, Sn, K)
    print('  held-out cross-half log-likelihood K=1..%d: %s -> K = %d' % (a.kmax, ' '.join('%.1f' % (h - rows_k[0]) for h in rows_k), K))
    _, lg = xd_ll(Y, Sn, pi, mu, Vk); resp = np.exp(lg - lg.max(1, keepdims=True)); resp /= resp.sum(1, keepdims=True); comp = resp.argmax(1)
    # chain components with Ashman D < 2 along the line between their means -> groups
    Dm = np.zeros((K, K))
    for i in range(K):
        for j in range(i + 1, K):
            u = mu[j] - mu[i]; nu = np.linalg.norm(u); u = u / max(nu, 1e-12); Dm[i, j] = Dm[j, i] = nu / np.sqrt(max(0.5 * (u @ Vk[i] @ u + u @ Vk[j] @ u), 1e-12))
    grp = connected_components((Dm < 2) & ~np.eye(K, dtype=bool), directed=False)[1]
    ug = np.unique(grp); noise_sd = np.sqrt(np.median(Sig.diagonal() / np.median(na[ok])))
    print('  components (share, mean on dir 1, deconvolved sd on dir 1; noise sd of a median sub-class mean %.3f):' % noise_sd)
    for k in np.argsort(mu[:, 0]): print('    comp %d: %.3f  %+.3f  %.3f  group %d' % (k + 1, pi[k], mu[k, 0], np.sqrt(Vk[k][0, 0]), np.where(ug == grp[k])[0][0] + 1))
    print('  -> %d group(s); adjacent separations D: %s' % (len(ug), ' '.join('%.1f' % Dm[i, j] for i, j in zip(np.argsort(mu[:, 0])[:-1], np.argsort(mu[:, 0])[1:]))))
    # group nature: reproducible internal spread (reliability of the group's leading direction among its sub-classes vs the same null)
    okidx = np.where(ok)[0]; sgrp = -np.ones(S, int); sgrp[okidx] = np.searchsorted(ug, grp[comp])
    for gi in range(len(ug)):
        m = okidx[sgrp[okidx] == gi]; share = na[m].sum() / na[ok].sum()
        if len(m) < 10: print('  group %d: %d sub-classes (%.1f %% of members): too few for an internal test' % (gi + 1, len(m), 100 * share)); continue
        lg_, _ = reliab(Me[m], Mo[m], w[m]); kind = 'CONTINUUM (reproducible internal spread)' if lg_[0] > thr else 'DISCRETE (no reproducible internal spread)'
        extra = ''
        if GT is not None:
            mm = np.isin(sid, m); swm = mm & GT['isw']
            extra = ' | TRUTH swing share %.2f, junk share %.2f, mean swing %.1f' % (GT['isw'][mm].mean(), GT['junk'][mm].mean(), np.nanmean(GT['swing'][swm]) if swm.any() else np.nan)
        print('  group %d: %d sub-classes (%.1f %% of members), internal reliability %.3f vs null %.3f -> %s%s' % (gi + 1, len(m), 100 * share, lg_[0], thr, kind, extra))
    summary[tag] = dict(d=d, K=K, ngroups=len(ug), lam=lam, thr=thr, Ys=Ys, Yp=Yp, ok=ok, sgrp=sgrp, mu=mu, pi=pi, V=Vk)
    # children along the leading reproducible coordinate: equal-width bins over sub-classes (1-99 % weighted range), members inherit
    y1 = Ys[:, 0]; o = np.argsort(y1); cw = np.cumsum(na[o]) / na.sum(); lo, hi = y1[o][np.searchsorted(cw, 0.01)], y1[o][np.searchsorted(cw, 0.99)]
    edges = np.linspace(lo, hi, a.nchild + 1); ch_s = np.clip(np.digitize(y1, edges[1:-1]), 0, a.nchild - 1); child = ch_s[sid] + 1
    yp1 = Yp[:, 0]; lo2, hi2 = np.percentile(yp1, [1, 99]); child_p = np.clip(np.digitize(yp1, np.linspace(lo2, hi2, a.nchild + 1)[1:-1]), 0, a.nchild - 1) + 1
    out = np.zeros(NR, int); out[p - 1] = child; np.savetxt(os.path.join(a.out, 'children_%s.txt' % tag), out, fmt='%d')
    rs = np.corrcoef(Ye[ok, 0], Yo[ok, 0])[0, 1]
    print('  children: %d equal-width bins of direction 1 over sub-classes, populations %s; even/odd sub-class coordinate corr %.3f'
          % (a.nchild, ' '.join(str((child == c).sum()) for c in range(1, a.nchild + 1)), rs))
    if GT is not None:
        sw = GT['isw']; e1 = eta2(GT['swing'][sw], child[sw]); e2 = eta2(GT['swing'][sw], child_p[sw])
        mus = [np.nanmean(GT['swing'][sw & (child == c)]) if (sw & (child == c)).any() else np.nan for c in range(1, a.nchild + 1)]
        shr = [GT['isw'][child == c].mean() if (child == c).any() else np.nan for c in range(1, a.nchild + 1)]
        print('  TRUTH children swing eta2 %.3f (adjusted %.3f) vs the same direction per particle %.3f; per child mean swing (swing share): %s'
              % (*e1, e2[0], ' '.join('%.1f(%.0f%%)' % (m_, 100 * s_) for m_, s_ in zip(mus, shr))))
np.savez(os.path.join(a.out, 'subcls_embed.npz'), p=p, eo=eo, sid=sid, us=us, nS=nS, sbin=sbin, Z=Z, res=res, mres=mres,
         **({'conf': GT['conf']} if GT is not None else {}))
# junk readout: sub-class mean residual ratio (latent-model residual / mean-only residual) as a robust z-score
rr = res / np.maximum(mres, 1e-12); rs_, _ = gmeans(np.c_[rr, mres], sid, S); med = np.median(rs_, 0); mad = 1.4826 * np.median(np.abs(rs_ - med), 0) + 1e-12; zj = (rs_ - med) / mad
for j, nm in enumerate(('residual ratio', 'mean residual')):
    fl = zj[:, j] > 3.5
    line = '\njunk flag on the sub-class %s (robust z > 3.5): %d sub-classes, %d members' % (nm, fl.sum(), nS[fl].sum())
    if GT is not None and GT['junk'].any():
        inf = fl[sid]; line += ' | TRUTH junk recall %.2f, flagged-member junk share %.2f' % (inf[GT['junk']].mean(), GT['junk'][inf].mean() if inf.any() else 0)
    print(line)
if a.refs and GT is not None:
    sw = GT['isw']
    for item in a.refs.split(','):
        nm, fn = item.split('='); lab = np.loadtxt(fn).astype(int)[p - 1]; m = sw & (lab > 0)
        print('REFERENCE %s on the same members: swing eta2 %.3f (%d labels)' % (nm, eta2(GT['swing'][m], lab[m])[0], len(np.unique(lab[lab > 0]))))
