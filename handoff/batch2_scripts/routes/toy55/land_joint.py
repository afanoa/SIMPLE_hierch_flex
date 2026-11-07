#!/usr/bin/env python3
"""land_joint.py: JOINT map + particle fit of an occupancy landscape p(t) on a physical basis grid (t_k, maps V_k), the coordinate's calibration
lam(t) (monotone, knots) and its noise scale c. Per fold (windows on lam of pixel set S, maps from the other pixel set) and half stack:
  particle term: lam_S,i ~ sum_k p_k N(lam(t_k), c / I_S,i)
  map term:      window map M_w ~ a(shell) * sum_k qhat_w(k) V_k, qhat_w(k) = mean over the window's particles of their model posterior
                 P(k | lam_S,i); per-shell noise from the two half stacks' difference; a(shell) profiled out
  prior:         lam_pen ||D2 G||^2, G = -ln p
The window maps measure how much the window compositions overlap, i.e. the noise the deconvolution needs: they constrain c and lam(t),
which the particle term alone cannot separate from the width of p. Input: a land_e1b.py (synthetic, --gt) or land_real.py run directory.
usage: land_joint.py --run DIR --kind e1b|real [data args as the run] --basis LIST|--swing DIR --out OUT [--nw 12 --lam_pen 1e-3,1e-2 --nboot 10]"""
import os, sys, glob, argparse, json, numpy as np
from scipy.optimize import minimize
from scipy.special import ndtr
ap = argparse.ArgumentParser()
for k in ('run', 'kind', 'sim', 'mr', 'cache', 'poses', 'out'): ap.add_argument('--' + k, required=True)
ap.add_argument('--basis', default=''); ap.add_argument('--swing', default=''); ap.add_argument('--gt', default=''); ap.add_argument('--codes', default='swing:200:230')
ap.add_argument('--nw', type=int, default=12); ap.add_argument('--lam_pen', default='1e-4,1e-3,1e-2'); ap.add_argument('--nboot', type=int, default=0); ap.add_argument('--nknots', type=int, default=7)
ap.add_argument('--crop', type=int, default=48); ap.add_argument('--mskdiam', type=float, default=230.); ap.add_argument('--smpd', type=float, default=1.464); ap.add_argument('--kv', type=float, default=300.)
ap.add_argument('--cs', type=float, default=2.7); ap.add_argument('--fraca', type=float, default=0.1); ap.add_argument('--ctfabs', type=int, default=0); ap.add_argument('--eps', type=float, default=0.05)
ap.add_argument('--lp', type=float, default=10.); ap.add_argument('--seed', type=int, default=0); ap.add_argument('--tmax', type=float, default=-1., help='real: use basis steps with travel <= TMAX')
a = ap.parse_args(); os.makedirs(a.out, exist_ok=True); rng = np.random.default_rng(a.seed)
def log(*s): print(*s, flush=True)
def kv(l): return dict(t.split('=', 1) for t in l.split() if '=' in t)
def rd(p):
    h = np.fromfile(p, np.int32, 256); n_ = int(h[0]); return np.fromfile(p, np.float32, offset=1024 + int(h[23]), count=n_ ** 3).reshape(n_, n_, n_).astype(np.float64)
# ---------------- run data: particles (master rows, half, lam per pixel set, Fisher), truth if synthetic
if a.kind == 'e1b':
    Z = np.load(a.run + '/e1b_particles.npz'); mrow = Z['rows']; half = Z['eo'].astype(int); tt = Z['t']
else:
    Z = np.load(a.run + '/real_particles.npz'); mrow = Z['mrow']; half = Z['half'].astype(int); tt = np.full(len(mrow), -1)
lamA, IA, lamB, IB = Z['lamA'], Z['IA'], Z['lamB'], Z['IB']; N = len(mrow)
# ---------------- basis grid
if a.basis:
    BL = [l.split() for l in open(a.basis) if l.strip()]; BF = [l[1] for l in BL]; tk = np.arange(len(BF), dtype=float)
else:
    tr = [l.split() for l in open(a.swing + '/truth_swing.txt') if l.strip() and not l.startswith('#')]; BF = [glob.glob('%s/conf%s_*.mrc' % (a.swing, l[1]))[0] for l in tr]; tk = np.array([float(l[3]) for l in tr])
    if a.tmax > 0: keep = tk <= a.tmax + 1e-6; BF = [f for f, k_ in zip(BF, keep) if k_]; tk = tk[keep]
K = len(BF)
# ---------------- Fourier data (crop c), halves = eo (one stack) or the two half stacks
lamb = 12.2643247 / np.sqrt(a.kv * 1e3 * (1. + a.kv * 1e3 * 0.978466e-6))
def rotm(e1, e2, e3):
    a1, a2, a3 = np.radians(e1), np.radians(e2), np.radians(e3); c1, si1, c2, si2, c3, si3 = np.cos(a1), np.sin(a1), np.cos(a2), np.sin(a2), np.cos(a3), np.sin(a3)
    Z1 = np.zeros((len(e1), 3, 3)); Z1[:, 0, 0] = c1; Z1[:, 0, 1] = -si1; Z1[:, 1, 0] = si1; Z1[:, 1, 1] = c1; Z1[:, 2, 2] = 1
    Y2 = np.zeros((len(e1), 3, 3)); Y2[:, 0, 0] = c2; Y2[:, 0, 2] = si2; Y2[:, 2, 0] = -si2; Y2[:, 2, 2] = c2; Y2[:, 1, 1] = 1
    Z3 = np.zeros((len(e1), 3, 3)); Z3[:, 0, 0] = c3; Z3[:, 0, 1] = -si3; Z3[:, 1, 0] = si3; Z3[:, 1, 1] = c3; Z3[:, 2, 2] = 1
    return Z1 @ Y2 @ Z3
SIMS = a.sim.split(','); n = int(np.frombuffer(open(SIMS[0] + '/all.mrcs', 'rb').read(1024), np.int32)[0]); c = a.crop
R0 = a.mskdiam / 2 / a.smpd; fc = np.fft.fftfreq(c) * c / n; fy2, fx2 = np.meshgrid(fc, fc, indexing='ij'); kk0 = np.stack([fx2, fy2, np.zeros_like(fx2)], -1).reshape(-1, 3) * n
r2d = np.hypot(kk0[:, 0], kk0[:, 1]); pcols = np.where(r2d <= c // 2 - 0.5)[0]; psh = np.round(r2d[pcols]).astype(int)
ppar = np.mod(np.round(kk0[pcols, 0]).astype(int) + np.round(kk0[pcols, 1]).astype(int), 2) == 1
klp = int(round(n * a.smpd / a.lp))
f3 = np.fft.fftfreq(c) * c / n; kx, ky, kz = np.meshgrid(f3, f3, f3, indexing='ij'); sh = np.minimum(np.round(np.sqrt(kx ** 2 + ky ** 2 + kz ** 2) * n).astype(int), c // 2).ravel()
g3 = np.arange(c) - c // 2; r3 = np.sqrt(g3[:, None, None] ** 2 + g3[None, :, None] ** 2 + g3[None, None, :] ** 2) * n / c; msk3 = np.clip((R0 + 4 - r3) / 8., 0, 1).astype(np.float32)
ns_ = c // 2 + 1; khi = c // 2 - 1; KV = np.stack([np.arange(c ** 3) % c, (np.arange(c ** 3) // c) % c, np.arange(c ** 3) // (c * c)], 1); KV = np.where(KV >= c // 2, KV - c, KV).astype(np.float64) / n
FB = (sh >= 2) & (sh <= min(klp, khi)); COLS = {'A': np.where(~ppar)[0], 'B': np.where(ppar)[0]}
if a.kind == 'e1b':
    NMp = sum(1 for _ in open(a.poses)); need = np.zeros(NMp, bool); need[mrow] = True; P_all = np.zeros((NMp, 8))
    for r, l in enumerate(open(a.poses)):
        if need[r]: d = kv(l); P_all[r] = [float(d.get(k, 0)) for k in ('e1', 'e2', 'e3', 'x', 'y', 'dfx', 'dfy', 'angast')]
    HR = [(np.where(half == h)[0], None) for h in (0, 1)]   # (particle indices, cache rows = master rows)
else:
    NMp = int(mrow.max()) + 1; need = np.zeros(NMp, bool); need[mrow] = True; P_all = np.zeros((NMp, 8))
    for r, l in enumerate(open(a.poses)):
        if r < NMp and need[r]: d = kv(l); P_all[r] = [float(d.get(k, 0)) for k in ('e1', 'e2', 'e3', 'x', 'y', 'dfx', 'dfy', 'angast')]
    HR = []
    for h in (0, 1):
        mr = np.loadtxt(a.mr.split(',')[h]).astype(int) - 1; pos = {m: i for i, m in enumerate(mr)}; pi_ = np.where(half == h)[0]; HR.append((pi_, np.array([pos[m] for m in mrow[pi_]])))
def load_half(h):
    pi_, crow = HR[h]; mr_ = mrow[pi_]; crow = mr_ if crow is None else crow; Fc = np.load('%s_%s_c%d.npy' % (a.cache, 'A' if a.kind == 'e1b' else 'AB'[h], c), mmap_mode='r')
    pp = P_all[mr_]; nn = len(pi_); Wv = np.zeros((nn, len(pcols)), np.complex64); Wd = np.zeros((nn, len(pcols)), np.float32); idx = np.zeros((nn, len(pcols)), np.int32)
    kf = (np.hypot(fx2, fy2) / a.smpd).ravel()[pcols]; ang = np.arctan2(fy2, fx2).ravel()[pcols]; o = np.argsort(crow); inv = np.empty_like(o); inv[o] = np.arange(nn)
    for b in range(0, nn, 20000):
        s_ = slice(b, b + 20000); p = pp[s_]; F = np.asarray(Fc[crow[s_]]).reshape(len(p), -1)[:, pcols]
        df = 0.5 * ((p[:, 5] + p[:, 6])[:, None] + (p[:, 5] - p[:, 6])[:, None] * np.cos(2. * (ang[None] - np.radians(p[:, 7])[:, None]))) * 1e4
        chi = np.pi * lamb * df * kf[None] ** 2 - 0.5 * np.pi * a.cs * 1e7 * lamb ** 3 * kf[None] ** 4; C = (np.sqrt(1. - a.fraca ** 2) * np.sin(chi) + a.fraca * np.cos(chi)).astype(np.float32)
        if a.ctfabs: C = np.abs(C)
        ramp = np.exp(2j * np.pi * (fx2.ravel()[pcols][None] * p[:, 3, None] + fy2.ravel()[pcols][None] * p[:, 4, None])).astype(np.complex64)
        Wv[s_] = C * np.conj(ramp) * F; Wd[s_] = C ** 2
        q = np.round(np.einsum('qj,nij->nqi', kk0[pcols], rotm(p[:, 0], p[:, 1], p[:, 2]))).astype(np.int32); okq = (np.abs(q) <= c // 2 - 1).all(2); q = q % c
        ii = q[..., 2] * c * c + q[..., 1] * c + q[..., 0]; ii[~okq] = c ** 3; idx[s_] = ii
    return dict(Wv=Wv, Wd=Wd, idx=idx, pi=pi_, n=nn)
HS = [load_half(h) for h in (0, 1)]; log('%s: %d particles (halves %d / %d), %d basis maps (t %.0f..%.0f)' % (a.run, N, HS[0]['n'], HS[1]['n'], K, tk[0], tk[-1]))
neq0 = 5000
def rec_raw(H, s_, cs):
    cc = COLS[cs]; ii = H['idx'][np.ix_(s_, cc)].ravel(); m = c ** 3 + 1; Wv_ = H['Wv'][np.ix_(s_, cc)]
    return np.bincount(ii, Wv_.real.ravel(), m)[:-1] + 1j * np.bincount(ii, Wv_.imag.ravel(), m)[:-1], np.bincount(ii, H['Wd'][np.ix_(s_, cc)].ravel(), m)[:-1]
EPS = None
def rec(H, s_, cs):
    global EPS
    num, den = rec_raw(H, s_, cs)
    if EPS is None:
        d0 = rec_raw(H, np.arange(H['n']), cs)[1] * (neq0 / H['n']); EPS = np.array([np.median(d0[(sh == k) & (d0 > 0)]) if ((sh == k) & (d0 > 0)).any() else 1. for k in range(ns_)]) * a.eps
    V_ = num / (den + EPS[sh] * max(len(s_), 1) / neq0); v = np.fft.fftshift(np.real(np.fft.ifftn(V_.reshape(c, c, c)))) * msk3; return np.fft.fftn(np.fft.ifftshift(v)).ravel()
def dock(Mp, Ct, m_=FB):
    a_ = Mp[m_] * np.conj(Ct[m_]); k_ = KV[m_]; t = np.zeros(3)
    for _ in range(8):
        z = a_ * np.exp(-2j * np.pi * (k_ @ t)); g_ = 2 * np.pi * (k_ * z.imag[:, None]).sum(0); H_ = -(2 * np.pi) ** 2 * (k_.T * z.real) @ k_
        try: st_ = np.linalg.solve(H_, g_)
        except np.linalg.LinAlgError: break
        t = t - st_
        if np.abs(st_).max() < 1e-4: break
    if not np.all(np.isfinite(t)) or np.abs(t).max() > 3.: return Mp
    return Mp * np.exp(-2j * np.pi * (KV @ t))
def crop_four(v):
    nb = v.shape[0]; sb = np.r_[0:c // 2, nb - c // 2:nb]; Fv = np.fft.fftn(np.fft.ifftshift(v)); vc = np.fft.fftshift(np.real(np.fft.ifftn(Fv[np.ix_(sb, sb, sb)]))) * msk3
    return np.fft.fftn(np.fft.ifftshift(vc)).ravel()
BAS = np.array([crop_four(rd(f)) for f in BF]); BM = BAS.mean(0); shF = sh[FB]; BF_ = BAS[:, FB]
Bs = np.array([np.real(BF_[:, shF == s] @ np.conj(BF_[:, shF == s]).T) for s in range(ns_)])   # per shell Gram of the basis
NVOX = np.bincount(shF, minlength=ns_).astype(float)   # voxels per shell (noise per voxel = shell noise power / NVOX)
# ---------------- held-out windows and their maps
FOLDS = (('A', 'B'), ('B', 'A')); LS = {'A': (lamA, IA), 'B': (lamB, IB)}; FW = []
for fs, ft in FOLDS:
    ls_, Is_ = LS[fs]; qs = np.percentile(ls_, np.linspace(0, 100, a.nw + 1)); wl = np.clip(np.searchsorted(qs, ls_, side='right') - 1, 0, a.nw - 1)
    Cw = np.zeros((2, a.nw, ns_, K)); Ew = np.zeros((2, a.nw, ns_)); Dw = np.zeros((a.nw, ns_)); Mh = []
    for w in range(a.nw):
        mm = []
        for h, H in enumerate(HS):
            loc = np.where(wl[H['pi']] == w)[0]; Mw = dock(rec(H, loc, ft), BM)[FB]; mm.append(Mw)
            for s in range(ns_):
                ms = shF == s; Cw[h, w, s] = np.real(BF_[:, ms] @ np.conj(Mw[ms])); Ew[h, w, s] = np.sum(np.abs(Mw[ms]) ** 2)
        for s in range(ns_): ms = shF == s; Dw[w, s] = 0.5 * np.sum(np.abs(mm[0][ms] - mm[1][ms]) ** 2)   # per-map noise power per shell
    grp = [(h, w, np.where((half == h) & (wl == w))[0]) for h in (0, 1) for w in range(a.nw)]
    FW.append(dict(ls=ls_, Is=Is_, wl=wl, Cw=Cw, Ew=Ew, sig2=np.maximum(Dw.mean(0), 1e-30), grp=grp))
log('windows: %d per fold, maps from the test pixels, %d shells in the band' % (a.nw, int(FB.sum() and len(np.unique(shF)))))
# ---------------- model
kn = np.linspace(tk[0], tk[-1], a.nknots); Hk = np.array([[np.interp(t_, kn, np.eye(a.nknots)[j]) for j in range(a.nknots)] for t_ in tk])
D2 = np.diff(np.eye(K), 2, axis=0); PEN = D2.T @ D2; svalid = [s for s in range(ns_) if Bs[s].trace() > 0]
def unpack(th):
    G = th[:K]; p = np.exp(-(G - G.min())); p /= p.sum(); l0, dl = th[K], np.logaddexp(0., th[K + 1:K + a.nknots]); lknot = l0 + np.r_[0., np.cumsum(dl)]; lt = Hk @ lknot; cc = np.exp(th[-1]); return G, p, lt, cc
def post_full(ig): raise SystemExit('subset fits not implemented')
def objective(th, lam_pen, sel=None, parts=False):
    G, p, lt, cc = unpack(th); tot_p = 0.; tot_m = 0.
    for F in FW:
        ls_, Is_ = F['ls'], F['Is']; m_ = slice(None) if sel is None else sel
        d = ls_[m_, None] - lt[None]; ivar = (Is_[m_] / cc)[:, None]; lg = -0.5 * d * d * ivar + 0.5 * np.log(ivar) + np.log(np.maximum(p, 1e-300))[None]
        mx = lg.max(1, keepdims=True); den = np.exp(lg - mx).sum(1, keepdims=True); post = np.exp(lg - mx) / den; tot_p += -(mx[:, 0] + np.log(den[:, 0])).sum()
        for h, w, ig in F['grp']:
            if sel is not None: ig = np.intersect1d(ig, np.where(sel)[0] if sel.dtype == bool else sel)
            if len(ig) == 0: continue
            q = post[ig].mean(0) if sel is None else post_full(ig)
            num = (F['Cw'][h, w][svalid] @ q); qbq = np.einsum('k,skl,l->s', q, Bs[svalid], q)
            tot_m += np.sum((F['Ew'][h, w][svalid] - num ** 2 / np.maximum(qbq, 1e-30)) * NVOX[svalid] / (2. * F['sig2'][svalid]))   # a(shell) profiled out
    val = tot_p / len(FW) + tot_m / len(FW) + lam_pen * N * (G @ PEN @ G)
    return (val, tot_p, tot_m) if parts else val
def fit(lam_pen, th0, sel=None, maxiter=300):
    f = lambda th: objective(th, lam_pen, sel)
    def fg(th):
        v = f(th); g = np.zeros_like(th); e = 1e-4
        for j in range(len(th)): tj = th.copy(); tj[j] += e; g[j] = (f(tj) - v) / e
        return v, g
    r = minimize(fg, th0, jac=True, method='L-BFGS-B', options=dict(maxiter=maxiter)); return r.x, r.fun
# start: flat p, lam(t) from the windows' mean lam vs equal spacing, c = 1
lam_all = 0.5 * (lamA + lamB); th0 = np.r_[np.zeros(K), np.percentile(lam_all, 2), np.log(np.expm1(np.maximum(np.diff(np.percentile(lam_all, np.linspace(2, 98, a.nknots))), 1e-4))), 0.]
res = dict(K=K, tk=tk.tolist(), fits={})
pt = np.bincount(tt[tt >= 0], minlength=K) / max((tt >= 0).sum(), 1) if (tt >= 0).any() else None
def sm3(p): return np.convolve(p, np.ones(3), 'same') / np.convolve(np.ones_like(p), np.ones(3), 'same')
for lp_ in [float(x) for x in a.lam_pen.split(',')]:
    th, fv = fit(lp_, th0); G, p, lt, cc = unpack(th); v_, tp, tm = objective(th, lp_, parts=True)
    out = dict(p=p.tolist(), lam_t=lt.tolist(), c=float(cc), obj=float(fv), part=float(tp), map=float(tm))
    if pt is not None:
        Gs = -np.log(np.maximum(sm3(p), 1e-12)); out['tv'] = 0.5 * float(np.abs(p - pt).sum()); out['low'] = float(p[:15].sum() + 0.5 * p[15]); out['minor'] = float(p[:4].sum())
        out['barrier'] = float(Gs[15] - min(Gs[5], Gs[25]))
    res['fits'][str(lp_)] = out
    log('RESULT lam_pen %.0e: c %.3f, lam(t) %+.3f..%+.3f, obj %.1f (particles %.1f, maps %.1f)%s' % (lp_, cc, lt[0], lt[-1], fv, tp, tm,
        '' if pt is None else ' | tv %.3f low %.3f minor %.3f barrier %+.2f' % (out['tv'], out['low'], out['minor'], out['barrier'])))
    th0 = th
if pt is not None:
    Gt = -np.log(np.maximum(sm3(pt), 1e-12)); res['truth'] = dict(p=pt.tolist(), low=float(pt[:15].sum() + 0.5 * pt[15]), minor=float(pt[:4].sum()), barrier=float(Gt[15] - min(Gt[5], Gt[25])))
    log('TRUTH low %.3f minor %.3f barrier %+.2f' % (res['truth']['low'], res['truth']['minor'], res['truth']['barrier']))
json.dump(res, open(a.out + '/joint.json', 'w'), indent=1); log('DONE %s' % a.out)
