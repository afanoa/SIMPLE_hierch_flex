#!/usr/bin/env python3
"""Size a site class from its posterior RANKING (site EM of minor_em3.py --post 1), not from the mixture weight. For a grid of N, the top-N
rows by P(class 2) (N/2 per half, ranked within each half) give half-maps (in-memory direct Fourier insertion at the half's poses, crop
48, as minor_em3.py) and the deviation D_h(N) = top-N_h - all_h inside the site mask; the cross-half correlation of D_A and D_B is measured
(1) in the EM band (shells 2..kb: selection imprint - particles chosen against the other half's class-2 reference reproduce it there even
as noise) and (2) in the HOLD-OUT band beyond it (shells kb + MARGIN .. crop Nyquist), where the selection is blind. For true members
the hold-out correlation grows with N and falls once the ranking dilutes (p(N)^2 N |D|^2 / (p(N)^2 N |D|^2 + s^2)), so N* = argmax of the
hold-out curve. Writes OUT/post_master.txt (P(class 2) per master row, -1 outside), OUT/curve.txt, OUT/top_star.txt (master rows in the
top N*) and, with --truth / --site, the true precision per N, the truth's state count and the AUC of the posterior.
usage: site_size.py --em DIR --start S_A,S_B --sim H_A,H_B --poses P_A,P_B --mr MRA,MRB --nmaster N --mask M.mrc --out DIR [--cache C
       --ctfabs 0 --em_lim 15 --margin 2 --grid 500,1000,... --truth sites_truth.txt --site 2 --state 1]"""
import os, argparse, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('--em', required=True); ap.add_argument('--start', required=True); ap.add_argument('--sim', required=True); ap.add_argument('--poses', required=True)
ap.add_argument('--mr', required=True); ap.add_argument('--nmaster', type=int, required=True); ap.add_argument('--mask', required=True); ap.add_argument('--out', required=True); ap.add_argument('--cache', default='')
ap.add_argument('--ctfabs', type=int, default=0); ap.add_argument('--em_lim', type=float, default=15.); ap.add_argument('--margin', type=int, default=2); ap.add_argument('--grid', default='500,1000,1500,2500,4000,6000,10000,15000,25000,40000,60000,90000')
ap.add_argument('--truth', default=''); ap.add_argument('--site', type=int, default=-1); ap.add_argument('--state', type=int, default=1); ap.add_argument('--crop', type=int, default=48); ap.add_argument('--mskdiam', type=float, default=230.)
ap.add_argument('--smpd', type=float, default=1.464); ap.add_argument('--kv', type=float, default=300.); ap.add_argument('--cs', type=float, default=2.7); ap.add_argument('--fraca', type=float, default=0.1); ap.add_argument('--conv', default='1,1,1,1,1,1')
ap.add_argument('--eps', type=float, default=1.); a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
cv = [float(v) for v in a.conv.split(',')]; convT, s1, yf, s3, ssh, sctf = bool(cv[0]), cv[1], cv[2], cv[3], cv[4], cv[5]
def log(*s): print(*s, flush=True)
def kv(l): return dict(t.split('=', 1) for t in l.split() if '=' in t)
def rd_mrc(p):
    h = np.fromfile(p, np.int32, 256); n = int(h[0]); return np.fromfile(p, np.float32, offset=1024 + int(h[23]), count=n ** 3).reshape(n, n, n).astype(np.float64)
lam = 12.2643247 / np.sqrt(a.kv * 1e3 * (1. + a.kv * 1e3 * 0.978466e-6))
def rotm(e1, e2, e3):
    a1, a2, a3 = np.radians(e1), np.radians(e2), np.radians(e3); c1, si1, c2, si2, c3, si3 = np.cos(a1), np.sin(a1), np.cos(a2), np.sin(a2), np.cos(a3), np.sin(a3)
    Z1 = np.zeros((len(e1), 3, 3)); Z1[:, 0, 0] = c1; Z1[:, 0, 1] = -si1; Z1[:, 1, 0] = si1; Z1[:, 1, 1] = c1; Z1[:, 2, 2] = 1
    Y2 = np.zeros((len(e1), 3, 3)); Y2[:, 0, 0] = c2; Y2[:, 0, 2] = si2; Y2[:, 2, 0] = -si2; Y2[:, 2, 2] = c2; Y2[:, 1, 1] = 1
    Z3 = np.zeros((len(e1), 3, 3)); Z3[:, 0, 0] = c3; Z3[:, 0, 1] = -si3; Z3[:, 1, 0] = si3; Z3[:, 1, 1] = c3; Z3[:, 2, 2] = 1
    return Z1 @ Y2 @ Z3
SIMS = a.sim.split(','); h0 = np.frombuffer(open(SIMS[0] + '/all.mrcs', 'rb').read(1024), np.int32); n = int(h0[0]); c = a.crop
x = np.arange(n) - n // 2; R0 = a.mskdiam / 2 / a.smpd; msk2 = np.clip((R0 + 4 - np.hypot(x[:, None], x[None, :])) / 8., 0, 1).astype(np.float32); sl = np.r_[0:c // 2, n - c // 2:n]
fc = np.fft.fftfreq(c) * c / n; fy2, fx2 = np.meshgrid(fc, fc, indexing='ij'); kk0 = np.stack([fx2, fy2, np.zeros_like(fx2)], -1).reshape(-1, 3) * n
f3 = np.fft.fftfreq(c) * c / n; kx, ky, kz = np.meshgrid(f3, f3, f3, indexing='ij'); sh = np.minimum(np.round(np.sqrt(kx ** 2 + ky ** 2 + kz ** 2) * n).astype(int), c // 2).ravel()
ns = c // 2 + 1; m3 = c ** 3 + 1; kb = min(c // 2 - 1, int(round(n * a.smpd / a.em_lim))); kt = c // 2 - 1
r2 = np.round(np.hypot(kk0[:, 0], kk0[:, 1])).astype(int); cols = np.where((r2 >= 2) & (r2 <= kt))[0]
kf = (np.hypot(fx2, fy2) / a.smpd).ravel()[cols]; ang = np.arctan2(fy2, fx2).ravel()[cols]
g3 = np.arange(c) - c // 2; r3 = np.sqrt(g3[:, None, None] ** 2 + g3[None, :, None] ** 2 + g3[None, None, :] ** 2) * n / c; msk3 = np.clip((R0 + 4 - r3) / 8., 0, 1)
log('crop %d (Nyquist %.1f A); EM band shells 2..%d (to %.1f A), hold-out band %d..%d (%.1f-%.1f A)' % (c, n * a.smpd / (c // 2), kb, n * a.smpd / kb, kb + a.margin, kt, n * a.smpd / (kb + a.margin), n * a.smpd / kt))
def load_half(simdir, posef, rows, tag):
    Pz = np.array([[float(kv(l).get(k, 0)) for k in ('e1', 'e2', 'e3', 'x', 'y', 'dfx', 'dfy', 'angast')] for l in open(posef) if l.strip()])[rows]
    hh = np.frombuffer(open(simdir + '/all.mrcs', 'rb').read(1024), np.int32); nz, nsym = int(hh[2]), int(hh[23])
    cache = '%s_%s_c%d.npy' % (a.cache, tag, c) if a.cache else ''; Fc = np.load(cache, mmap_mode='r') if cache and os.path.isfile(cache) else None
    I = None if Fc is not None else np.memmap(simdir + '/all.mrcs', np.float32, 'r', offset=1024 + nsym, shape=(nz, n, n))
    N = len(rows); Wv = np.zeros((N, len(cols)), np.complex64); Wd = np.zeros((N, len(cols)), np.float32); ix = np.zeros((N, len(cols)), np.int32)
    for b in range(0, N, 10000):
        rr = rows[b:b + 10000]; pp = Pz[b:b + 10000]
        if Fc is not None: F = np.asarray(Fc[rr]).reshape(len(rr), -1)[:, cols]
        else: blk = np.array(I[rr], np.float32) * msk2; F = np.fft.fft2(np.fft.ifftshift(blk, axes=(1, 2)))[:, sl][:, :, sl].reshape(len(rr), -1)[:, cols]
        df = 0.5 * ((pp[:, 5] + pp[:, 6])[:, None] + (pp[:, 5] - pp[:, 6])[:, None] * np.cos(2. * (ang[None] - np.radians(pp[:, 7])[:, None]))) * 1e4
        chi = np.pi * lam * df * kf[None] ** 2 - 0.5 * np.pi * a.cs * 1e7 * lam ** 3 * kf[None] ** 4; C = (sctf * (np.sqrt(1. - a.fraca ** 2) * np.sin(chi) + a.fraca * np.cos(chi))).astype(np.float32)
        if a.ctfabs: C = np.abs(C)
        ramp = np.exp(ssh * 2j * np.pi * (fx2.ravel()[cols][None] * pp[:, 3, None] + fy2.ravel()[cols][None] * pp[:, 4, None])).astype(np.complex64)
        Wv[b:b + len(rr)] = C * np.conj(ramp) * F; Wd[b:b + len(rr)] = C ** 2
        Rm = rotm(s1 * pp[:, 0], pp[:, 1], s3 * pp[:, 2]); kk = kk0[cols] * np.array([1., yf, 1.])
        p3 = np.einsum('qj,nij->nqi', kk, Rm) if convT else np.einsum('qj,nji->nqi', kk, Rm); q = np.round(p3).astype(np.int32)
        ok = (np.abs(q) <= c // 2 - 1).all(2); q = q % c; ixb = q[..., 2] * c * c + q[..., 1] * c + q[..., 0]; ixb[~ok] = c ** 3; ix[b:b + len(rr)] = ixb
    log('%s: %d rows loaded' % (tag, N)); return dict(Wv=Wv, Wd=Wd, ix=ix)
ST = [np.loadtxt(f, dtype=int) for f in a.start.split(',')]; PO = [np.loadtxt('%s/post_soft_%s.txt' % (a.em, h))[:, 1] for h in 'AB']; MR = [np.loadtxt(f).astype(int) - 1 for f in a.mr.split(',')]
rows = [np.where(s > 0)[0] for s in ST]; D = [load_half(SIMS[h], a.poses.split(',')[h], rows[h], 'AB'[h]) for h in (0, 1)]
post = np.full(a.nmaster, -1.)
for h in (0, 1): post[MR[h][rows[h]]] = PO[h][rows[h]]
np.savetxt(a.out + '/post_master.txt', post, fmt='%.4f')
neq = 5000.
def shell_eps(d):
    den = np.bincount(d['ix'].ravel(), d['Wd'].ravel(), m3)[:-1] * (neq / len(d['Wd'])); return np.array([np.median(den[(sh == k) & (den > 0)]) if ((sh == k) & (den > 0)).any() else 1. for k in range(ns)]) * a.eps
eps = 0.5 * (shell_eps(D[0]) + shell_eps(D[1])); epsv = np.r_[eps[sh], eps[-1]]
def wmap(d, sel):
    ii = d['ix'][sel].ravel(); num = np.bincount(ii, d['Wv'][sel].real.ravel(), m3) + 1j * np.bincount(ii, d['Wv'][sel].imag.ravel(), m3)
    V = num / (np.bincount(ii, d['Wd'][sel].ravel(), m3) + epsv * len(sel) / neq); V[-1] = 0.; return V
def to_real(V): return np.fft.fftshift(np.real(np.fft.ifftn(V[:-1].reshape(c, c, c))))
def crop_real(v):
    F = np.fft.fftn(np.fft.ifftshift(v)); return np.fft.fftshift(np.real(np.fft.ifftn(F[np.ix_(sl, sl, sl)]))) * (c / v.shape[0]) ** 3
M = np.clip(crop_real(rd_mrc(a.mask)), 0, 1) * msk3
allm = [to_real(wmap(D[h], np.arange(len(rows[h])))) for h in (0, 1)]
def band_corr(x, y, lo, hi):
    X = np.fft.fftn(np.fft.ifftshift(x)).ravel(); Y = np.fft.fftn(np.fft.ifftshift(y)).ravel(); m = (sh >= lo) & (sh <= hi)
    return float(np.real(np.vdot(X[m], Y[m])) / max(np.linalg.norm(X[m]) * np.linalg.norm(Y[m]), 1e-30))
def band_amp(x, y, lo, hi):   # cross-half deviation amplitude in a band (noise-free in expectation)
    X = np.fft.fftn(np.fft.ifftshift(x)).ravel(); Y = np.fft.fftn(np.fft.ifftshift(y)).ravel(); m = (sh >= lo) & (sh <= hi); return float(np.sqrt(max(np.real(np.vdot(X[m], Y[m])), 0.)) / c ** 1.5)
truth = None
if a.truth and a.site >= 0:
    T = np.genfromtxt(a.truth, comments='#'); tv = T[:, 3 + a.site]; truth = ((tv > 0) == bool(a.state)) & (T[:, 2] == 0)   # discrete 0/1 or continuous (sign)
    inem = post >= 0; nt = int((truth & inem).sum()); y = truth[inem]; p = post[inem]
    o = np.argsort(p); rk = np.empty(len(o)); rk[o] = np.arange(len(o)); auc = (rk[y].sum() - y.sum() * (y.sum() - 1) / 2) / max(y.sum() * (~y).sum(), 1)
    log('truth: %d rows of %d in the EM hold the site state %d (%.3f); posterior AUC %.3f; mixture weight (sum of P) %.0f' % (nt, inem.sum(), a.state, nt / inem.sum(), auc, p.sum()))
grid = [int(v) for v in a.grid.split(',') if int(v) <= 0.7 * sum(len(r) for r in rows)]; res = []
amp_true = None
if truth is not None:   # deviation of the TRUE members: reference amplitudes for the imprint check (EM band / hold-out band)
    Dt = [(to_real(wmap(D[h], np.where(truth[MR[h][rows[h]]])[0])) - allm[h]) * M for h in (0, 1)]
    amp_true = (band_amp(Dt[0], Dt[1], 2, kb), band_amp(Dt[0], Dt[1], kb + a.margin, kt))
    log('true members: %d / %d per half; deviation amplitude EM band %.3g, hold-out band %.3g' % (truth[MR[0][rows[0]]].sum(), truth[MR[1][rows[1]]].sum(), amp_true[0], amp_true[1]))
order = [np.argsort(-PO[h][rows[h]]) for h in (0, 1)]
for N in grid:
    Dm = []
    for h in (0, 1):
        sel = order[h][:N // 2]; Dm.append((to_real(wmap(D[h], sel)) - allm[h]) * M)
    ci = band_corr(Dm[0], Dm[1], 2, kb); co = band_corr(Dm[0], Dm[1], kb + a.margin, kt); occ = float(np.sqrt(max((Dm[0] * Dm[1]).sum(), 0.)))
    aho = band_amp(Dm[0], Dm[1], kb + a.margin, kt)
    line = 'N %6d: P(class 2) >= %.3f | cross-half corr of the site deviation: EM band %.3f, hold-out band %.3f | deviation amplitude %.3g, hold-out %.3g' % (N, min(PO[h][rows[h]][order[h][N // 2 - 1]] for h in (0, 1)), ci, co, occ, aho)
    prec = np.nan
    if truth is not None:
        top = np.zeros(a.nmaster, bool)
        for h in (0, 1): top[MR[h][rows[h][order[h][:N // 2]]]] = True
        prec = truth[top].mean(); line += ' | TRUE precision %.3f, recall %.3f' % (prec, (truth & top).sum() / max(nt, 1))
        line += ' | apparent occupancy (amplitude / true members\' amplitude): EM band %.2f, hold-out %.2f' % (band_amp(Dm[0], Dm[1], 2, kb) / max(amp_true[0], 1e-30), aho / max(amp_true[1], 1e-30))
    res.append((N, ci, co, occ, prec, aho)); log(line)
r_ = np.array(res); ns_ = int(r_[np.argmax(r_[:, 2]), 0]); log('N* = %d (argmax of the hold-out cross-half correlation %.3f)' % (ns_, r_[:, 2].max()))
eq = r_[:, 5] * r_[:, 0] / max(r_[0, 5], 1e-30)   # equivalent pure members: hold-out amplitude x N over the smallest-N amplitude (saturates at the true count for a binary site)
log('equivalent pure members (hold-out amplitude x N / amplitude at N=%d): %s -> %s' % (r_[0, 0], ' '.join('%d:%.0f' % (n_, e_) for n_, e_ in zip(r_[:, 0], eq)), 'saturates (max %.0f at N %d)' % (eq.max(), r_[np.argmax(eq), 0]) if np.argmax(eq) < len(eq) - 1 else 'still rising at the largest N (continuum or weak ranking)'))
np.savetxt(a.out + '/curve.txt', np.c_[r_, eq], fmt='%.5f', header='N corr_em_band corr_holdout amplitude true_precision holdout_amplitude equivalent_pure_members')
top = np.zeros(a.nmaster, int)
for h in (0, 1): top[MR[h][rows[h][order[h][:ns_ // 2]]]] = 1
np.savetxt(a.out + '/top_star.txt', top, fmt='%d'); open(a.out + '/nstar.txt', 'w').write('%d\n' % ns_)
if truth is not None: log('top N* = %d: TRUE precision %.3f, recall %.3f (truth count %d)' % (ns_, truth[top > 0].mean(), (truth & (top > 0)).sum() / max(nt, 1), nt))
