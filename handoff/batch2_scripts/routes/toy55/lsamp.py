#!/usr/bin/env python3
"""lsamp.py (landscape sampling, step 2): samples of the reproducible variability of a tree's score space, tail to tail at equal steps, the
number of samples set by resolvability.
(1) AXES: the centred score vectors (lsamp_score.py; per particle minus its mean over the used leaves) against the noise covariance G of the
    scores (what a particle with no heterogeneity shows): generalised eigenvectors C v = lam G v. lam - 1 = the map-space signal variance
    along the axis in units of the per-particle noise. Computed per eo half (each half scored against the other half's leaf maps) and pooled;
    an axis is reproducible when lam - 1 > 0 in both halves and its loadings agree (G-metric cosine).
(2) SAMPLES: coordinates t = centred scores . v on NAX axes; per axis n samples at equal steps from the 1 % to the 99 % quantile (--ends centre:
    sample k centred at q01 + k (q99 - q01) / (n - 1), the end samples open-ended; tail: the end samples = the particles beyond q01 / q99,
    n - 2 equal-width interior samples), grids for 2-3 axes. A cell whose half-maps do not reproduce to the score band limit (FSC at EM_LIM
    < MINFSC) merges into its nearest kept neighbour (only such cells merge).
(3) RESOLVABILITY: n = the largest candidate for which every pair of neighbouring samples is distinguishable by splitbench16's test: maps of
    each pair from the equal N of its smaller member in each half (crop 48, direct Fourier insertion, one fixed regulariser), docked by a 3-D
    translation to the other half's consensus on the tested band, statistic = best 3-shell window of the cross-half FSC of the pair
    difference over shells beyond the score band (kb + RB_MARGIN .. Nyquist of the crop: the coordinates never saw them), threshold = the
    maximum over NSHUF view x defocus stratified label shuffles.
(4) CHECKS at the chosen sampling: every sample reproducible (FSC at EM_LIM); neighbours distinguishable (margins); end samples vs a direct
    reconstruction of the tail particles (beyond q01 / q99) by their cross-half position along the end-to-end direction on the tested band
    (0 = low end sample, 1 = high end sample); neighbour map distances d2 = cross-half difference energy / sqrt(signal_p signal_q) (equal N,
    docked) on the full band and on the tested band, and their CV. With --gt: per-sample truth (coded coordinate, closed share, conformers).
usage: lsamp.py --scores DIR --sim .. --mr .. --cache CP --poses POSES --out OUT [--use r22_,r21_ | --exclude r_s3] [--naxes 1]
       [--ends centre|tail] [--ns 4,6,8,10,12,16,20] [--gt GT --codes swing:200:230] [--focus MASK] [--ctfabs 0]"""
import os, argparse, numpy as np
from scipy.linalg import eigh
from scipy.optimize import linear_sum_assignment
ap = argparse.ArgumentParser()
ap.add_argument('--scores', required=True); ap.add_argument('--sim', required=True); ap.add_argument('--mr', required=True); ap.add_argument('--cache', required=True)
ap.add_argument('--poses', required=True); ap.add_argument('--out', required=True); ap.add_argument('--use', default=''); ap.add_argument('--exclude', default='')
ap.add_argument('--naxes', type=int, default=0, help='0: report the axes only (no sampling)'); ap.add_argument('--axes', default='', help='comma list of axis numbers (1-based) to sample instead of the leading NAXES')
ap.add_argument('--ends', default='centre'); ap.add_argument('--ns', default='4,6,8,10,12,14,16,20'); ap.add_argument('--grid_ns', default='2,3,4,5,6,8')
ap.add_argument('--qlo', type=float, default=1.); ap.add_argument('--em_lim', type=float, default=20.); ap.add_argument('--rb_margin', type=int, default=2); ap.add_argument('--bandlo', type=float, default=40.)
ap.add_argument('--nshuf', type=int, default=10); ap.add_argument('--minfsc', type=float, default=0.143); ap.add_argument('--minhalf', type=int, default=150); ap.add_argument('--focus', default='')
ap.add_argument('--gt', default=''); ap.add_argument('--codes', default=''); ap.add_argument('--crop', type=int, default=48); ap.add_argument('--mskdiam', type=float, default=230.)
ap.add_argument('--smpd', type=float, default=1.464); ap.add_argument('--kv', type=float, default=300.); ap.add_argument('--cs', type=float, default=2.7); ap.add_argument('--fraca', type=float, default=0.1)
ap.add_argument('--ctfabs', type=int, default=0); ap.add_argument('--eps', type=float, default=0.05); ap.add_argument('--nview', type=int, default=40); ap.add_argument('--seed', type=int, default=0)
ap.add_argument('--tag', default='lsamp'); ap.add_argument('--select', default='axis', help='axis: neighbours must be ordered along the axis (matched-filter step beyond the null); fsc: splitbench16 window statistic'); ap.add_argument('--nit', type=int, default=0, help='score-loop iterations: rescore against the chosen samples\' maps (cross-half, score band), new axes, new search')
a = ap.parse_args(); os.makedirs(a.out, exist_ok=True); rng = np.random.default_rng(a.seed)
def log(*s): print(*s, flush=True)
def kv(l): return dict(t.split('=', 1) for t in l.split() if '=' in t)
# ---------------- score space
S = np.load(a.scores + '/scores.npy', mmap_mode='r'); M = np.load(a.scores + '/meta.npz'); names = [str(x) for x in M['names']]; EO = M['eo']; OWN = M['own']; NM = len(OWN)
use = [j for j, nm in enumerate(names) if (not a.use or any(nm.startswith(p) for p in a.use.split(','))) and nm not in a.exclude.split(',')]
K = len(use); rows = np.where(np.isin(OWN, use))[0]; X = np.asarray(S[rows][:, use], np.float64); ok = np.isfinite(X).all(1) & (EO[rows] >= 0); rows, X = rows[ok], X[ok]
Xc = X - X.mean(1, keepdims=True); eo = EO[rows]; N = len(rows)
log('%s: %d particles of %d leaves (%s); score band to %.1f A' % (a.tag, N, K, ' '.join(names[j] for j in use), float(M['em_lim'])))
Gh = [M['gram0'][np.ix_(use, use)], M['gram1'][np.ix_(use, use)]]
w_, U_ = np.linalg.eigh(np.eye(K) - 1. / K); Bc = U_[:, w_ > 0.5]   # orthonormal basis of the centred subspace
def geig(Z, G):
    Cs = np.cov(Z.T); lam, W = eigh(Bc.T @ Cs @ Bc, Bc.T @ G @ Bc); o = np.argsort(-lam); return lam[o], Bc @ W[:, o]
lamh, Vh = [], []
for e in (0, 1): l_, v_ = geig(Xc[eo == e], Gh[e]); lamh.append(l_); Vh.append(v_)
Gp = (Gh[0] * (eo == 0).sum() + Gh[1] * (eo == 1).sum()) / N; lam, V = geig(Xc, Gp); T = Xc @ V   # pooled axes; T = coordinates (noise variance 1 per axis)
def gcos(u, v): return abs(u @ Gp @ v) / np.sqrt((u @ Gp @ u) * (v @ Gp @ v))
na = min(8, K - 1); agree = [[max(gcos(V[:, j], Vh[e][:, i]) for i in range(min(na + 2, K - 1))) for e in (0, 1)] for j in range(na)]
mh = np.array([[gcos(Vh[0][:, i], Vh[1][:, j]) for j in range(na)] for i in range(na)]); ia, ib = linear_sum_assignment(-mh)
sig = np.maximum(lam - 1., 0.); share = sig / max(sig.sum(), 1e-12)
log('axes (generalised eigenvectors of the centred score covariance against the noise covariance; lam - 1 = signal / noise per particle):')
log('  axis  lam-1 pooled | half eo0 / eo1 | loading agreement with the halves | signal share (cum)')
for j in range(na):
    log('  %4d  %7.3f | %7.3f / %7.3f | %.3f / %.3f | %.3f (%.3f)' % (j + 1, lam[j] - 1, lamh[0][j] - 1, lamh[1][j] - 1, agree[j][0], agree[j][1], share[j], share[:j + 1].sum()))
log('  half-to-half matched axes (eo0 axis: eo1 axis cos): %s' % ' '.join('%d:%d %.2f' % (i + 1, j + 1, mh[i, j]) for i, j in zip(ia, ib)))
def subcos(m):   # mean cosine of the principal angles between the two halves' leading m-dim subspaces (G metric)
    Ln = np.linalg.cholesky(Gp + 1e-12 * np.eye(K)).T; A_ = np.linalg.qr(Ln @ Vh[0][:, :m])[0]; B_ = np.linalg.qr(Ln @ Vh[1][:, :m])[0]
    return np.linalg.svd(A_.T @ B_, compute_uv=False)
log('  subspace agreement (principal-angle cosines, halves): %s' % ' | '.join('m=%d: %s' % (m, ' '.join('%.2f' % v for v in subcos(m))) for m in (1, 2, 3) if m < K))
# leaves along the axes (mean coordinate of each leaf's particles, in sd of the particle coordinates)
sdT = T.std(0)
for j in range(min(3, na)):
    mu = sorted(((T[OWN[rows] == u, j].mean() / sdT[j], names[u], int((OWN[rows] == u).sum())) for u in use))
    q = np.percentile(T[:, j] / sdT[j], [0.5, 1, 99, 99.5])
    log('  axis %d: particles 0.5/1/99/99.5 %% at %.2f/%.2f/%.2f/%.2f sd; leaf means %s' % (j + 1, q[0], q[1], q[2], q[3], ' '.join('%s:%.2f' % (nm, m) for m, nm, _ in mu)))
# ---------------- truth (optional)
GT = None
if a.gt:
    g_ = np.loadtxt(a.gt, comments='#', ndmin=2); conf = np.full(NM, -1); clo = np.zeros(NM); conf[g_[:, 0].astype(int) - 1] = g_[:, 1].astype(int); clo[g_[:, 0].astype(int) - 1] = g_[:, 2] if g_.shape[1] > 2 else 0
    GT = dict(conf=conf[rows], closed=clo[rows]); codes = {}
    for cd in [w for w in a.codes.split(',') if w]:
        nm, lo, hi = cd.split(':'); cf = GT['conf']; codes[nm] = np.where((cf >= int(lo)) & (cf <= int(hi)), cf - int(lo), np.nan)
    GT['codes'] = codes
    def nmi(x, y):
        _, xi = np.unique(x, return_inverse=True); _, yi = np.unique(y, return_inverse=True); J = np.zeros((xi.max() + 1, yi.max() + 1)); np.add.at(J, (xi, yi), 1); J /= J.sum()
        px, py = J.sum(1), J.sum(0); nz = J > 0; return (J[nz] * np.log(J[nz] / np.outer(px, py)[nz])).sum() / max(np.sqrt((px * np.log(px)).sum() * (py * np.log(py)).sum()), 1e-12)
    for j in range(min(6, na)):
        s = 'TRUTH axis %d:' % (j + 1)
        for nm, x in codes.items():
            m = np.isfinite(x)
            if m.sum() > 100: s += ' |corr(t, %s)| %.3f (%d particles);' % (nm, abs(np.corrcoef(T[m, j], x[m])[0, 1]), m.sum())
        if GT['closed'].any(): s += ' closed share in the lower / upper 5 %%: %.3f / %.3f (all %.3f);' % (GT['closed'][T[:, j] < np.percentile(T[:, j], 5)].mean(), GT['closed'][T[:, j] > np.percentile(T[:, j], 95)].mean(), GT['closed'].mean())
        b10 = np.clip(((T[:, j] - T[:, j].min()) / np.ptp(T[:, j]) * 10).astype(int), 0, 9); s += ' NMI(conformer; 10 bins) %.3f' % nmi(GT['conf'], b10)
        log(s)
np.savez(a.out + '/axes.npz', lam=lam, V=V, lam0=lamh[0], lam1=lamh[1], rows=rows, T=T[:, :min(6, na)].astype(np.float32), use=np.array(use), names=np.array(names))
if a.naxes <= 0: log('axes only (--naxes 0): done'); raise SystemExit
AX = [int(v) - 1 for v in a.axes.split(',')] if a.axes else list(range(a.naxes)); Tz = T[:, AX] / sdT[AX]; nax = len(AX)
# ---------------- particle data for the maps (crop 48): halves = half stacks (A / B) or the eo partition (full stack)
lamb = 12.2643247 / np.sqrt(a.kv * 1e3 * (1. + a.kv * 1e3 * 0.978466e-6))
def rotm(e1, e2, e3):
    a1, a2, a3 = np.radians(e1), np.radians(e2), np.radians(e3); c1, si1, c2, si2, c3, si3 = np.cos(a1), np.sin(a1), np.cos(a2), np.sin(a2), np.cos(a3), np.sin(a3)
    Z1 = np.zeros((len(e1), 3, 3)); Z1[:, 0, 0] = c1; Z1[:, 0, 1] = -si1; Z1[:, 1, 0] = si1; Z1[:, 1, 1] = c1; Z1[:, 2, 2] = 1
    Y2 = np.zeros((len(e1), 3, 3)); Y2[:, 0, 0] = c2; Y2[:, 0, 2] = si2; Y2[:, 2, 0] = -si2; Y2[:, 2, 2] = c2; Y2[:, 1, 1] = 1
    Z3 = np.zeros((len(e1), 3, 3)); Z3[:, 0, 0] = c3; Z3[:, 0, 1] = -si3; Z3[:, 1, 0] = si3; Z3[:, 1, 1] = c3; Z3[:, 2, 2] = 1
    return Z1 @ Y2 @ Z3
SIMS = a.sim.split(','); full = a.mr == 'full'; n = int(np.frombuffer(open(SIMS[0] + '/all.mrcs', 'rb').read(1024), np.int32)[0]); c = a.crop
x_ = np.arange(n) - n // 2; R0 = a.mskdiam / 2 / a.smpd; fc = np.fft.fftfreq(c) * c / n; fy2, fx2 = np.meshgrid(fc, fc, indexing='ij')
kk0 = np.stack([fx2, fy2, np.zeros_like(fx2)], -1).reshape(-1, 3) * n; r2d = np.hypot(kk0[:, 0], kk0[:, 1]); pcols = np.where(r2d <= c // 2 - 0.5)[0]
kb0 = min(c // 2 - 1, int(round(n * a.smpd / a.em_lim))); bcols = np.where((np.round(r2d[pcols]) >= 2) & (np.round(r2d[pcols]) <= kb0))[0]   # score-band pixels (positions within pcols)
f3 = np.fft.fftfreq(c) * c / n; kx, ky, kz = np.meshgrid(f3, f3, f3, indexing='ij'); sh = np.minimum(np.round(np.sqrt(kx ** 2 + ky ** 2 + kz ** 2) * n).astype(int), c // 2).ravel()
g3 = np.arange(c) - c // 2; r3 = np.sqrt(g3[:, None, None] ** 2 + g3[None, :, None] ** 2 + g3[None, None, :] ** 2) * n / c; msk3 = np.clip((R0 + 4 - r3) / 8., 0, 1).astype(np.float32)
if a.focus and a.focus != 'auto':
    from scipy.ndimage import zoom as _zoom
    _h = np.fromfile(a.focus, np.int32, 256); _m = np.fromfile(a.focus, np.float32, offset=1024 + int(_h[23]), count=int(_h[0] * _h[1] * _h[2])).reshape(int(_h[2]), int(_h[1]), int(_h[0]))
    msk3 = msk3 * np.clip(_zoom(_m, c / _m.shape[0], order=1), 0, 1).astype(np.float32)[:c, :c, :c]; log('focus mask %s: %.1f %% of the crop grid above 0.5' % (a.focus, 100. * np.mean(msk3 > 0.5)))
ns_ = c // 2 + 1; khi = c // 2 - 1; kb = min(khi, int(round(n * a.smpd / a.em_lim))); klo = max(int(round(n * a.smpd / a.bandlo)), 3); kt = min(kb + a.rb_margin, khi - 2)
TB = (sh >= kt) & (sh <= khi); FB = (sh >= 2) & (sh <= khi); KV = np.stack([np.arange(c ** 3) % c, (np.arange(c ** 3) // c) % c, np.arange(c ** 3) // (c * c)], 1); KV = np.where(KV >= c // 2, KV - c, KV).astype(np.float64) / n
P_all = np.array([[float(kv(l).get(k, 0)) for k in ('e1', 'e2', 'e3', 'x', 'y', 'dfx', 'dfy', 'angast')] for l in open(a.poses)]); assert len(P_all) == NM
inrow = np.full(NM, -1); inrow[rows] = np.arange(N)
def load_half(tag, crow, mrow):   # crow: cache rows, mrow: their master rows (all in the node)
    Fc = np.load('%s_%s_c%d.npy' % (a.cache, tag, c), mmap_mode='r'); pp = P_all[mrow]; nn = len(crow)
    Wv = np.zeros((nn, len(pcols)), np.complex64); Wd = np.zeros((nn, len(pcols)), np.float32); idx = np.zeros((nn, len(pcols)), np.int32); F2 = np.zeros((nn, len(bcols)), np.float32)
    kf = (np.hypot(fx2, fy2) / a.smpd).ravel()[pcols]; ang = np.arctan2(fy2, fx2).ravel()[pcols]
    for b in range(0, nn, 20000):
        s_ = slice(b, b + 20000); p = pp[s_]; F = np.asarray(Fc[crow[s_]]).reshape(len(p), -1)[:, pcols]
        df = 0.5 * ((p[:, 5] + p[:, 6])[:, None] + (p[:, 5] - p[:, 6])[:, None] * np.cos(2. * (ang[None] - np.radians(p[:, 7])[:, None]))) * 1e4
        chi = np.pi * lamb * df * kf[None] ** 2 - 0.5 * np.pi * a.cs * 1e7 * lamb ** 3 * kf[None] ** 4; C = (np.sqrt(1. - a.fraca ** 2) * np.sin(chi) + a.fraca * np.cos(chi)).astype(np.float32)
        if a.ctfabs: C = np.abs(C)
        ramp = np.exp(2j * np.pi * (fx2.ravel()[pcols][None] * p[:, 3, None] + fy2.ravel()[pcols][None] * p[:, 4, None])).astype(np.complex64)
        Wv[s_] = C * np.conj(ramp) * F; Wd[s_] = C ** 2; F2[s_] = (np.abs(F[:, bcols]) ** 2)
        q = np.round(np.einsum('qj,nij->nqi', kk0[pcols], rotm(p[:, 0], p[:, 1], p[:, 2]))).astype(np.int32); okq = (np.abs(q) <= c // 2 - 1).all(2); q = q % c
        ii = q[..., 2] * c * c + q[..., 1] * c + q[..., 0]; ii[~okq] = c ** 3; idx[s_] = ii
    ph, th = np.radians(pp[:, 0]), np.radians(pp[:, 1]); v = np.c_[np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)]
    i = np.arange(a.nview) + 0.5; z = 1 - 2 * i / a.nview; rr = np.sqrt(np.maximum(1 - z * z, 0)); gg = np.pi * (3 - np.sqrt(5)) * i
    vb = np.argmax(v @ np.c_[rr * np.cos(gg), rr * np.sin(gg), z].T, 1); dfm = 0.5 * (pp[:, 5] + pp[:, 6]); dq = np.searchsorted(np.quantile(dfm, [1 / 3., 2 / 3.]), dfm)
    return dict(Wv=Wv, Wd=Wd, idx=idx, F2=F2, ti=inrow[mrow], st=vb * 3 + dq, n=nn)
if full: HS = [load_half('A', rows[eo == e], rows[eo == e]) for e in (0, 1)]
else:
    HS = []
    for h, tg in enumerate('AB'):
        mr = np.loadtxt(a.mr.split(',')[h]).astype(int) - 1; sel = np.where(inrow[mr] >= 0)[0]; HS.append(load_half(tg, sel, mr[sel]))
log('map data: halves %d / %d particles (%s), crop %d, pixels %d; bands: score 2..%d (%.1f A), tested %d..%d (%.1f..%.1f A), full %d..%d' % (HS[0]['n'], HS[1]['n'], 'eo partition' if full else 'half stacks', c, len(pcols), kb, n * a.smpd / kb, kt, khi, n * a.smpd / kt, n * a.smpd / khi, klo, khi))
def rec_raw(H, sel):
    ii = H['idx'][sel].ravel(); m = c ** 3 + 1; num = np.bincount(ii, H['Wv'][sel].real.ravel(), m)[:-1] + 1j * np.bincount(ii, H['Wv'][sel].imag.ravel(), m)[:-1]
    return num, np.bincount(ii, H['Wd'][sel].ravel(), m)[:-1]
neq0 = 5000
def shell_eps(H):
    den = rec_raw(H, np.arange(H['n']))[1] * (neq0 / H['n']); return np.array([np.median(den[(sh == k) & (den > 0)]) if ((sh == k) & (den > 0)).any() else 1. for k in range(ns_)]) * a.eps
EPS = 0.5 * (shell_eps(HS[0]) + shell_eps(HS[1]))
def rec(H, sel):   # map of particles sel (indices into the half), one fixed regulariser per N (as splitbench: eps x N / N0), masked in real space
    num, den = rec_raw(H, sel); V_ = num / (den + EPS[sh] * max(len(sel), 1) / neq0); v = np.fft.fftshift(np.real(np.fft.ifftn(V_.reshape(c, c, c)))) * msk3
    return np.fft.fftn(np.fft.ifftshift(v)).ravel()
CONS = [rec(H, np.arange(H['n'])) for H in HS]
def fsc(X1, X2): return np.bincount(sh, np.real(X1 * np.conj(X2)), ns_) / np.maximum(np.sqrt(np.bincount(sh, np.abs(X1) ** 2, ns_) * np.bincount(sh, np.abs(X2) ** 2, ns_)), 1e-30)
def dock(Mp, Ct, m_):
    a_ = Mp[m_] * np.conj(Ct[m_]); k_ = KV[m_]; t = np.zeros(3)
    for _ in range(8):
        z = a_ * np.exp(-2j * np.pi * (k_ @ t)); g_ = 2 * np.pi * (k_ * z.imag[:, None]).sum(0); H_ = -(2 * np.pi) ** 2 * (k_.T * z.real) @ k_
        try: st_ = np.linalg.solve(H_, g_)
        except np.linalg.LinAlgError: break
        t = t - st_
        if np.abs(st_).max() < 1e-4: break
    if not np.all(np.isfinite(t)) or np.abs(t).max() > 3.: return Mp
    return Mp * np.exp(-2j * np.pi * (KV @ t))
def tstat(DA, DB): f = fsc(DA, DB); return max(float(f[k:k + 3].mean()) for k in range(kt, khi - 1))
def pairmaps(selp, selq):   # both groups in both halves from the equal N of the smaller one, docked to the OTHER half's consensus on the tested band
    out = []
    for h, H in enumerate(HS):
        nn = min(len(selp[h]), len(selq[h]))
        for s_ in (selp[h], selq[h]): out.append(dock(rec(H, rng.choice(s_, nn, replace=False)), CONS[1 - h], TB))
    return out   # Ap Aq Bp Bq
def members(lab, g): return [np.where(lab[h] == g)[0] for h in range(2)]
def ends_dirs(lab, gs, nper_unused=None):   # per axis: the low / high end cells (extreme cell means) and their full-N half-maps docked on the tested band
    cen = {g: Tz[cell_of_lab(lab, g)].mean(0) for g in gs}; out = []
    for j in range(nax):
        lo_, hi_ = min(gs, key=lambda g: cen[g][j]), max(gs, key=lambda g: cen[g][j])
        E = [{g: dock(rec(HS[h], np.where(lab[h] == g)[0]), CONS[1 - h], TB) for g in (lo_, hi_)} for h in range(2)]
        out.append((lo_, hi_, E[0][hi_] - E[0][lo_], E[1][hi_] - E[1][lo_]))
    return out, cen
def cell_of_lab(lab, g): return evaluate.cellvec == g
def axstat(Ap, Aq, Bp, Bq, DA, DB, band):   # along-axis step from p to q in units of the end-to-end difference, cross-half (q - p projected on the other half's ends)
    den = np.real(np.vdot(DB[band], DA[band])); return 0.5 * (np.real(np.vdot(DB[band], (Aq - Ap)[band])) + np.real(np.vdot(DA[band], (Bq - Bp)[band]))) / max(den, 1e-30)
def null_thr(lab, groups, prs=(), dirs=()):   # stratified label shuffles at the global equal N: max pair statistic (all pairs) and max |along-axis step| (the neighbour pairs)
    nulls = []; nax_ = []; neq = min(min((lab[h] == g).sum() for g in groups) for h in range(2))
    for s in range(a.nshuf):
        pl = []
        for h, H in enumerate(HS):
            l_ = lab[h].copy()
            for k in np.unique(H['st']): ix = np.where(H['st'] == k)[0]; l_[ix] = rng.permutation(l_[ix])
            pl.append(l_)
        Mh = [{g: dock(rec(HS[h], rng.choice(np.where(pl[h] == g)[0], neq, replace=False)), CONS[1 - h], TB) for g in groups} for h in range(2)]
        gs = sorted(groups); nulls += [tstat(Mh[0][p] - Mh[0][q], Mh[1][p] - Mh[1][q]) for i, p in enumerate(gs) for q in gs[i + 1:]]
        for p, q, j in prs: nax_.append(abs(axstat(Mh[0][p], Mh[0][q], Mh[1][p], Mh[1][q], dirs[j][2], dirs[j][3], TB)))
    return max(nulls), float(np.median(nulls)), (max(nax_) if nax_ else 0.)
def auto_focus(cell):   # falloff-prior focus mask from the cells' cross-half variance (arm_mask2.py: per-shell amplitude matched to the mean, lp 15 A, core 20 % of the peak dilated ~4.4 A, falloff sd ~8.8 A)
    from scipy.ndimage import distance_transform_edt, binary_dilation, gaussian_filter
    lab = half_labels(cell); gs = list(np.unique(cell)); kl = int(round(n * a.smpd / 15.)); lp = (sh <= kl).astype(float)
    sph = np.clip((R0 + 4 - r3) / 8., 0, 1); Mh = [[rec(HS[h], np.where(lab[h] == g)[0]) for g in gs] for h in range(2)]; wts = np.array([(cell == g).sum() for g in gs], float); wts /= wts.sum()
    mA = sum(w * m for w, m in zip(wts, Mh[0])); mB = sum(w * m for w, m in zip(wts, Mh[1])); pm = np.maximum(np.bincount(sh, np.real(mA * np.conj(mB)), ns_), 1e-30)
    def real(X): return np.fft.fftshift(np.real(np.fft.ifftn((X * lp).reshape(c, c, c))))
    E, O = [], []
    for k in range(len(gs)):
        sc = np.sqrt(np.maximum(np.bincount(sh, np.real(Mh[0][k] * np.conj(Mh[1][k])), ns_), 1e-30) / pm); sc[:2] = 1.; E.append(real(Mh[0][k] / sc[sh])); O.append(real(Mh[1][k] / sc[sh]))
    Eb = sum(w * e for w, e in zip(wts, E)); Ob = sum(w * o for w, o in zip(wts, O)); v = gaussian_filter(sum(w * (e - Eb) * (o - Ob) for w, e, o in zip(wts, E, O)) * sph, 0.5)
    core = binary_dilation(v > 0.2 * v.max(), iterations=max(1, int(round(4.4 / (n * a.smpd / c))))); d = distance_transform_edt(~core) * (n * a.smpd / c)
    m = np.where(core, 1., np.exp(-d ** 2 / (2 * 8.8 ** 2))) * (sph > 0)
    log('auto focus mask from %d cells: core %.1f %% of the sphere, mask mean %.2f inside the sphere' % (len(gs), 100. * core[sph > 0.5].mean(), m[sph > 0.5].mean()))
    return m.astype(np.float32)
# ---------------- cells
def make_cells(nper):   # per-row cell index on the AX grid; returns labels (0-based cell ids on the full grid) and the per-axis edges
    idx = []
    for j in range(nax):
        t = Tz[:, j]; lo, hi = np.percentile(t, [a.qlo, 100 - a.qlo]); nj = nper[j]
        if a.ends == 'tail' and nj >= 3: e = np.linspace(lo, hi, nj - 1); idx.append(np.searchsorted(e, t, side='right'))
        else: w = (hi - lo) / (nj - 1); idx.append(np.clip(np.round((t - lo) / w), 0, nj - 1).astype(int))
    return np.ravel_multi_index(idx, nper), idx
def half_labels(cell): return [cell[H['ti']] for H in HS]
def cell_fsc(lab, g):
    s = members(lab, g)
    if min(len(s[0]), len(s[1])) < a.minhalf: return -1.
    f = fsc(rec(HS[0], s[0]), rec(HS[1], s[1])); return float(f[max(kb - 1, 2):kb + 2].mean())
def merge_sparse(cell, nper):   # cells that do not reproduce to the score band limit merge into the nearest kept neighbour on the grid
    cell = cell.copy(); merged = {}; lab = half_labels(cell); fs = {g: cell_fsc(lab, g) for g in np.unique(cell)}
    while True:
        gs = list(fs); bad = [g for g in gs if fs[g] < a.minfsc]
        if not bad or len(gs) < 2: return cell, fs, merged
        g = min(bad, key=lambda q: (cell == q).sum()); keep = [q for q in gs if q != g]
        cen = {q: Tz[cell == q].mean(0) for q in gs}; tgt = min(keep, key=lambda q: ((cen[q] - cen[g]) ** 2).sum()); cell[cell == g] = tgt; merged[g] = tgt
        del fs[g]; lab = half_labels(cell); fs[tgt] = cell_fsc(lab, tgt)
def neighbours(cell, nper):   # pairs of kept cells that touch on the grid (cells merged into another count as that cell)
    gs = np.unique(cell); grid_of = {}; full_ids = np.arange(int(np.prod(nper)))
    occ = {}
    for g in gs:
        for q in np.unique(np.ravel_multi_index(tuple(np.array(make_cells.last_idx)[:, cell == g]), nper)): occ[q] = g
    pr = set()
    for q, g in occ.items():
        iv = np.array(np.unravel_index(q, nper))
        for j in range(nax):
            iw = iv.copy(); iw[j] += 1
            if iw[j] < nper[j]:
                q2 = int(np.ravel_multi_index(iw, nper))
                if q2 in occ and occ[q2] != g: pr.add(tuple(sorted((g, occ[q2]))) + (j,))
    return sorted(pr)
def evaluate(nper, full_checks=False):
    cell0, idx = make_cells(nper); make_cells.last_idx = idx; cell, fs, merged = merge_sparse(cell0, nper); lab = half_labels(cell); gs = list(np.unique(cell))
    prs = neighbours(cell, nper); evaluate.cellvec = cell
    if len(gs) < 2 or not prs: return dict(ok=False, n=len(gs), cell=cell)
    dirs, cen = ends_dirs(lab, gs); prs = [(p, q, j) if cen[p][j] <= cen[q][j] else (q, p, j) for p, q, j in prs]   # oriented low -> high along the axis
    thr, med, athr = null_thr(lab, gs, prs, dirs); st = {}; ax = {}
    for p, q, j in prs:
        Ap, Aq, Bp, Bq = pairmaps(members(lab, p), members(lab, q)); st[(p, q, j)] = tstat(Ap - Aq, Bp - Bq); ax[(p, q, j)] = axstat(Ap, Aq, Bp, Bq, dirs[j][2], dirs[j][3], TB)
    okf = all(v > thr for v in st.values()); oka = all(v > athr for v in ax.values()); ok = oka if a.select == 'axis' else okf
    log('  n %s: %d samples (%d merged as not reproducible), %d neighbour pairs; splitbench16 statistic %s (null max %.3f) -> %s; along-axis steps %s (null max |step| %.3f) -> %s' % (
        'x'.join(map(str, nper)), len(gs), len(merged), len(prs), ' '.join('%.2f' % st[k] for k in prs), thr, 'all distinguishable' if okf else '%d not' % sum(v <= thr for v in st.values()),
        ' '.join('%.3f' % ax[k] for k in prs), athr, 'all ordered' if oka else '%d not' % sum(v <= athr for v in ax.values())))
    return dict(ok=ok, okf=okf, oka=oka, n=len(gs), cell=cell, fs=fs, merged=merged, st=st, ax=ax, thr=thr, athr=athr, prs=prs, nper=nper, idx=idx, dirs=dirs)
# ---------------- resolvability search
cands = [int(v) for v in (a.ns if nax == 1 else a.grid_ns).split(',')]
def search(tag):
    res = {}; best = None; nfail = 0
    log('resolvability search on axes %s (%s ends, quantiles %.1f..%.1f %%)%s:' % ('+'.join(str(j + 1) for j in AX), a.ends, a.qlo, 100 - a.qlo, tag))
    for nc in cands:
        r_ = evaluate([nc] * nax); res[nc] = r_
        if r_['ok']: best = nc; nfail = 0
        else: nfail += 1
        if best is not None and nfail >= 2: break   # two failures past the last success
    return res, best
def rescore(cell):   # every node particle against the samples' maps of the OTHER half (score band, per-shell scale and noise as lsamp_score.py)
    gs_ = list(np.unique(cell)); lab_ = half_labels(cell); MP = [[np.r_[rec(HS[h], np.where(lab_[h] == g)[0]), 0.].astype(np.complex64) for g in gs_] for h in range(2)]
    Sc = np.zeros((N, len(gs_))); Gn = [np.zeros((len(gs_), len(gs_))) for _ in range(2)]
    for h, H in enumerate(HS):
        ref = MP[1 - h]; own = np.searchsorted(gs_, lab_[h]); shp = np.round(r2d[pcols][bcols]).astype(int); Wv_ = H['Wv'][:, bcols]; Wd_ = H['Wd'][:, bcols]; ix_ = H['idx'][:, bcols]
        bn = np.zeros(ns_); bd = np.zeros(ns_)
        for b in range(0, H['n'], 20000):
            s_ = slice(b, b + 20000); Po = np.zeros((len(range(*s_.indices(H['n']))), len(bcols)), np.complex64)
            for k, g in enumerate(gs_):
                m = own[s_] == k
                if m.any(): Po[m] = ref[k][ix_[s_][m]]
            bn += np.bincount(np.broadcast_to(shp, Po.shape).ravel(), np.real(Po * np.conj(Wv_[s_])).ravel(), ns_); bd += np.bincount(np.broadcast_to(shp, Po.shape).ravel(), (Wd_[s_] * np.abs(Po) ** 2).ravel(), ns_)
        beta = np.clip(bn / np.maximum(bd, 1e-30), 0., 4.); bp = beta[shp].astype(np.float32); r_ = np.zeros(ns_); cn_ = np.zeros(ns_)
        for b in range(0, H['n'], 20000):
            s_ = slice(b, b + 20000); Po = np.zeros((len(range(*s_.indices(H['n']))), len(bcols)), np.complex64)
            for k, g in enumerate(gs_):
                m = own[s_] == k
                if m.any(): Po[m] = ref[k][ix_[s_][m]]
            res_ = H['F2'][s_] - 2. * np.real(bp * Po * np.conj(Wv_[s_])) + Wd_[s_] * (bp * np.abs(Po)) ** 2; r_ += np.bincount(np.broadcast_to(shp, res_.shape).ravel(), res_.ravel(), ns_); cn_ += np.bincount(np.broadcast_to(shp, res_.shape).ravel(), minlength=ns_)
        w = (1. / np.maximum(r_ / np.maximum(cn_, 1), 1e-12))[shp].astype(np.float32)
        for b in range(0, H['n'], 20000):
            s_ = slice(b, b + 20000); ng = min(1500, len(range(*s_.indices(H['n'])))); PG = np.zeros((ng, len(bcols), len(gs_)), np.complex64)
            for k in range(len(gs_)):
                P_ = ref[k][ix_[s_]] * bp; Sc[H['ti'][s_], k] = 0.5 * ((2. * np.real(P_ * np.conj(Wv_[s_])) - Wd_[s_] * np.abs(P_) ** 2) * w).sum(1); PG[:, :, k] = P_[:ng]
            X_ = PG * np.sqrt(w * Wd_[s_][:ng])[:, :, None]; Gn[h] += (np.einsum('npk,npl->kl', X_.real, X_.real) + np.einsum('npk,npl->kl', X_.imag, X_.imag)) / H['n'] * (H['n'] / ng) / max(1, (H['n'] + 19999) // 20000)
    return Sc, Gn
if a.focus == 'auto':
    msk_auto = auto_focus(make_cells([8] * nax if nax == 1 else [4] * nax)[0]); msk3 = msk3 * msk_auto; CONS = [rec(H, np.arange(H['n'])) for H in HS]
    np.save(a.out + '/focus_auto.npy', msk_auto)
res, best = search('')
for it in range(a.nit):
    if best is None: break
    cell_ = res[best]['cell']; Sc, Gn = rescore(cell_); Kc = Sc.shape[1]; Xr = Sc - Sc.mean(1, keepdims=True)
    w2, U2 = np.linalg.eigh(np.eye(Kc) - 1. / Kc); Bc2 = U2[:, w2 > 0.5]; Ghalf = Gn[0] if HS[0]['n'] else Gn[1]; Gp2 = 0.5 * (Gn[0] + Gn[1])
    lam2, W2 = eigh(Bc2.T @ np.cov(Xr.T) @ Bc2, Bc2.T @ Gp2 @ Bc2); o = np.argsort(-lam2); lam2 = lam2[o]; V2 = Bc2 @ W2[:, o]; T2 = Xr @ V2
    sgn = np.sign([np.corrcoef(T2[:, j], Tz[:, min(j, nax - 1)])[0, 1] for j in range(nax)]); Tz = T2[:, :nax] * sgn / T2[:, :nax].std(0)
    msg = 'score loop %d: rescored against %d sample maps; lam - 1 of the new axes %s' % (it + 1, Kc, ' '.join('%.2f' % v for v in lam2[:4] - 1))
    if GT is not None:
        for nm, x in GT['codes'].items():
            m = np.isfinite(x)
            if m.sum() > 100: msg += '; |corr(t, %s)| %s' % (nm, ' '.join('%.3f' % abs(np.corrcoef(Tz[m, j], x[m])[0, 1]) for j in range(nax)))
    log(msg); res, best = search(' after score loop %d' % (it + 1))
if best is None: log('no candidate has every neighbour distinguishable: the axes carry no resolvable landscape'); raise SystemExit
R = res[best]; cell = R['cell']; nper = R['nper']; lab = half_labels(cell); gs = list(np.unique(cell)); make_cells.last_idx = R['idx']
order = sorted(gs, key=lambda g: tuple(np.round(Tz[cell == g].mean(0), 3)))   # sample numbering by position (axis 1 first)
sid = {g: i + 1 for i, g in enumerate(order)}; out = np.zeros(NM, int); out[rows] = [sid[g] for g in cell]; np.savetxt(a.out + '/samples.txt', out, fmt='%d')
log('CHOSEN: %s per axis -> %d samples; labels %s/samples.txt (master rows)' % ('x'.join(map(str, nper)), len(gs), a.out))
# ---------------- checks
lines = []
for g in order:
    s = members(lab, g); f = fsc(rec(HS[0], s[0]), rec(HS[1], s[1])); r143 = next((n * a.smpd / k for k in range(2, khi + 1) if f[k] < 0.143), n * a.smpd / khi)
    lines.append('  sample %2d: %6d particles, centre %s sd, FSC at %.0f A %.2f, FSC 0.143 at %.1f A' % (sid[g], (cell == g).sum(), '/'.join('%+.2f' % v for v in Tz[cell == g].mean(0)), a.em_lim, R['fs'].get(g, np.nan), r143))
log('check 1, reproducibility:'); [log(l_) for l_ in lines]
log('check 2, neighbours distinguishable: splitbench16 statistic (null max %.3f) %s' % (R['thr'], ' '.join('%d-%d(ax%d):%.2f' % (sid[p], sid[q], j + 1, v) for (p, q, j), v in sorted(R['st'].items(), key=lambda kv: (kv[0][2], sid[kv[0][0]])))))
log('         along-axis steps (null max |step| %.3f) %s' % (R['athr'], ' '.join('%d-%d(ax%d):%.3f' % (sid[p], sid[q], j + 1, v) for (p, q, j), v in sorted(R['ax'].items(), key=lambda kv: (kv[0][2], sid[kv[0][0]])))))
evaluate.cellvec = cell; MAc = {g: dock(rec(HS[0], members(lab, g)[0]), CONS[1], TB) for g in gs}; MBc = {g: dock(rec(HS[1], members(lab, g)[1]), CONS[0], TB) for g in gs}
log('check 4a, positions along each axis (cross-half on the tested band, 0 / 1 = the low / high end sample) and the steps between neighbours:')
for j in range(nax):
    lo_, hi_, DA, DB = R['dirs'][j]; pos = {g: axstat(MAc[lo_], MAc[g], MBc[lo_], MBc[g], DA, DB, TB) for g in gs}; srt = sorted(gs, key=lambda g: Tz[cell == g, j].mean())
    stp = np.array([R['ax'][(p, q, jj)] for (p, q, jj) in R['prs'] if jj == j]); log('  axis %d: %s; neighbour steps mean %.3f CV %.2f' % (AX[j] + 1, ' '.join('%d:%.2f' % (sid[g], pos[g]) for g in srt), stp.mean(), stp.std() / max(abs(stp.mean()), 1e-12)))
def d2(Ap, Aq, Bp, Bq, band):
    num = np.real(np.vdot((Bp - Bq)[band], (Ap - Aq)[band])); sp = np.real(np.vdot(Bp[band], Ap[band])); sq = np.real(np.vdot(Bq[band], Aq[band])); return num / np.sqrt(max(sp * sq, 1e-30))
log('check 4, neighbour map distances d2 (full band %.0f-%.0f A | tested band):' % (n * a.smpd / 2, n * a.smpd / khi))
for j in range(nax):
    dd = [(sid[p], sid[q], d2(*pairmaps(members(lab, p), members(lab, q)), FB), d2(*pairmaps(members(lab, p), members(lab, q)), TB)) for (p, q, jj) in R['prs'] if jj == j]
    dd.sort()
    if dd:
        f1 = np.array([v[2] for v in dd]); f2 = np.array([v[3] for v in dd])
        log('  axis %d: %s; CV %.2f | %.2f' % (AX[j] + 1, ' '.join('%d-%d:%.3f|%.3f' % v for v in dd), f1.std() / max(abs(f1.mean()), 1e-12), f2.std() / max(abs(f2.mean()), 1e-12)))
# check 3: end samples vs the tails, position along the end-to-end direction (cross-half, tested band; 0 = low end, 1 = high end)
log('check 3, end samples vs direct reconstructions of the tails (position along the end-to-end direction, cross-half; 0 / 1 = the end samples):')
for j in range(nax):
    t = Tz[:, j]; lo, hi = np.percentile(t, [a.qlo, 100 - a.qlo]); cen = {g: Tz[cell == g].mean(0) for g in gs}
    other = [k for k in range(nax) if k != j]
    def endof(sign):   # end sample on axis j: the extreme cell among those central on the other axes
        cand = sorted(gs, key=lambda g: sign * cen[g][j]); return cand[-1]
    glo, ghi = endof(-1), endof(1); sets = {'end_lo': (cell == glo), 'end_hi': (cell == ghi), 'tail_lo': t < lo, 'tail_hi': t > hi, 'tail_lo0.5': t < np.percentile(t, a.qlo / 2), 'tail_hi0.5': t > np.percentile(t, 100 - a.qlo / 2)}
    MA = {k: rec(HS[0], np.where(v[HS[0]['ti']])[0]) for k, v in sets.items()}; MB = {k: rec(HS[1], np.where(v[HS[1]['ti']])[0]) for k, v in sets.items()}
    MA = {k: dock(v, CONS[1], TB) for k, v in MA.items()}; MB = {k: dock(v, CONS[0], TB) for k, v in MB.items()}
    for band, bn in ((TB, 'tested'), (FB, 'full')):
        DA, DB = MA['end_hi'] - MA['end_lo'], MB['end_hi'] - MB['end_lo']; den = np.real(np.vdot(DB[band], DA[band]))
        def pos(k): return 0.5 * (np.real(np.vdot(DB[band], (MA[k] - MA['end_lo'])[band])) + np.real(np.vdot(DA[band], (MB[k] - MB['end_lo'])[band]))) / den, (np.real(np.vdot(DB[band], (MA[k] - MA['end_lo'])[band])) - np.real(np.vdot(DA[band], (MB[k] - MB['end_lo'])[band]))) / den
        log('  axis %d (%s band): %s' % (AX[j] + 1, bn, ' '.join('%s %d ptcl %.2f (+-%.2f)' % (k, int(v.sum()), *pos(k)) for k, v in sets.items())))
# truth per sample
if GT is not None:
    log('TRUTH per sample:')
    for g in order:
        m = cell == g; s = '  sample %2d (%6d):' % (sid[g], m.sum())
        for nm, x in GT['codes'].items():
            xs = x[m]; f = np.isfinite(xs); s += ' %s %.1f+-%.1f (%.0f %% coded)' % (nm, np.nanmean(xs) if f.any() else np.nan, np.nanstd(xs) if f.any() else np.nan, 100 * f.mean())
        if GT['closed'].any(): s += ' closed %.3f (%d)' % (GT['closed'][m].mean(), int(GT['closed'][m].sum()))
        cc = np.bincount(GT['conf'][m].clip(0)); s += ' top conformers %s' % ' '.join('%d:%.2f' % (k, cc[k] / m.sum()) for k in np.argsort(-cc)[:3] if cc[k] > 0)
        log(s)
    for nm, x in GT['codes'].items():
        mus = []
        for g in order:
            m = (cell == g) & np.isfinite(x)
            if m.sum() >= 100 and m.sum() >= 0.5 * (cell == g).sum(): mus.append(x[m].mean())
        if len(mus) >= 2:
            st_ = np.diff(mus) if nax == 1 else np.array([0.])
            f = np.isfinite(x); W = sum(((x[(cell == g) & f] - x[(cell == g) & f].mean()) ** 2).sum() for g in gs if ((cell == g) & f).any()); e2 = 1 - W / ((x[f] - x[f].mean()) ** 2).sum()
            log('TRUTH %s: samples >= 50 %% coded span %.1f .. %.1f; eta2 %.3f%s' % (nm, min(mus), max(mus), e2, '' if nax > 1 else '; steps %s (CV %.2f)' % (' '.join('%.1f' % v for v in st_), np.std(st_) / max(abs(np.mean(st_)), 1e-12))))
    if GT['closed'].any():
        k = max(order, key=lambda g: GT['closed'][cell == g].mean()); log('TRUTH closed: best sample %d holds %d particles, closed share %.3f, %.1f %% of all closed' % (sid[k], (cell == k).sum(), GT['closed'][cell == k].mean(), 100 * GT['closed'][cell == k].sum() / max(GT['closed'].sum(), 1)))
    log('TRUTH NMI(conformer; sample) %.3f' % nmi(GT['conf'], cell))
log('DONE')
