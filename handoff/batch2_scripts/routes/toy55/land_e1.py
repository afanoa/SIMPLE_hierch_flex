#!/usr/bin/env python3
"""land_e1.py: TRUTH-FREE occupancy landscape of one tree node along its leading reproducible axis, with the truth used only to score it.
(1) COORDINATE: start = the node's leading score axis (lsamp2: per-particle-normalised scores vs the node leaves, generalised eigenvector
    against the noise covariance); then NIT passes of a per-particle LEAST-SQUARES position lam_i along the end-to-end direction
    D = M_hi - M_lo of the OTHER eo half (cell maps of the coordinate's ends, docked), images modelled as C beta (M0 + lam D) + noise
    (per-shell scale beta from the consensus, per-shell noise, per-particle noise level), computed on pixel set A and on set B (parity of
    kx + ky: independent noise) -> lam = Fisher-weighted mean, I_i = its Fisher information; the split z = (lam_A - lam_B) / sqrt(1/I_A +
    1/I_B) calibrates the noise model (var z = c, per view x defocus stratum).
(2) WINDOWS: NW equal-count windows of lam -> half-maps -> NNLS onto a PHYSICAL basis (here the true swing maps, per-shell amplitude
    matched, docked): window composition q_w(t). Map-level occupancy p_map(t) = sum_w P(w) q_w(t); calibration lam(t) (monotone piecewise
    linear) from mean lam of each window = sum_t q_w(t) lam(t).
(3) DECONVOLUTION: NPMLE of p(t) on the t grid with lam_i ~ N(lam(t), c / I_i) (truth-free: windows' lam(t), split c; oracle: the true
    E[lam | t] and the true residual c), bootstrap bands; naive = histogram of lam mapped through lam(t).
Known landscapes are imposed by subsampling the node's coded particles by their true t (--shape; contaminants keep their share).
usage: land_e1.py --scores DIR --use r22_ --sim S --mr full --cache CP --poses P --gt GT --codes swing:200:230 --basis LEAVES --out OUT
       [--shape sim|flat|broad|dw2|dw4|shoulder --family 0|1 (1: coded particles only) --nit 2 --nw 10]"""
import os, argparse, json, numpy as np
from scipy.linalg import eigh
from scipy.optimize import nnls
ap = argparse.ArgumentParser()
for k in ('scores', 'sim', 'mr', 'cache', 'poses', 'gt', 'out', 'basis'): ap.add_argument('--' + k, required=True)
ap.add_argument('--use', default=''); ap.add_argument('--codes', default='swing:200:230'); ap.add_argument('--shape', default='sim'); ap.add_argument('--family', type=int, default=0)
ap.add_argument('--nit', type=int, default=2); ap.add_argument('--nw', type=int, default=10); ap.add_argument('--qend', type=float, default=15.); ap.add_argument('--nboot', type=int, default=20)
ap.add_argument('--crop', type=int, default=48); ap.add_argument('--mskdiam', type=float, default=230.); ap.add_argument('--smpd', type=float, default=1.464); ap.add_argument('--kv', type=float, default=300.)
ap.add_argument('--cs', type=float, default=2.7); ap.add_argument('--fraca', type=float, default=0.1); ap.add_argument('--ctfabs', type=int, default=0); ap.add_argument('--eps', type=float, default=0.05)
ap.add_argument('--nview', type=int, default=40); ap.add_argument('--seed', type=int, default=0); ap.add_argument('--bstep', type=int, default=3, help='physical basis every BSTEP grid points')
ap.add_argument('--lp', type=float, default=10., help='resolution limit of the lam statistic and of the window unmixing (A)')
a = ap.parse_args(); os.makedirs(a.out, exist_ok=True); rng = np.random.default_rng(a.seed)
def log(*s): print(*s, flush=True)
def kv(l): return dict(t.split('=', 1) for t in l.split() if '=' in t)
def rd(p):
    h = np.fromfile(p, np.int32, 256); n_ = int(h[0]); return np.fromfile(p, np.float32, offset=1024 + int(h[23]), count=n_ ** 3).reshape(n_, n_, n_).astype(np.float64)
# ---------------- node particles and the starting axis (lsamp2's score space)
M = np.load(a.scores + '/meta.npz'); names = [str(x) for x in M['names']]; EO = M['eo']; OWN = M['own']; NM = len(OWN)
use = [j for j, nm in enumerate(names) if (not a.use or any(nm.startswith(p) for p in a.use.split(',')))]; K = len(use)
g_ = np.loadtxt(a.gt, comments='#', ndmin=2); conf = np.full(NM, -1); conf[g_[:, 0].astype(int) - 1] = g_[:, 1].astype(int)
cn, lo_, hi_ = a.codes.split(':'); lo_, hi_ = int(lo_), int(hi_); NG = hi_ - lo_ + 1; tg = np.arange(NG, dtype=float)
tcode = np.where((conf >= lo_) & (conf <= hi_), conf - lo_, -1)
rows = np.where(np.isin(OWN, use) & (EO >= 0))[0]
if a.family: rows = rows[tcode[rows] >= 0]
def target(name):
    if name == 'flat': g = np.zeros(NG)
    elif name.startswith('dw'): h = float(name[2:]); g = h * (((tg - 15.) / 10.) ** 2 - 1.) ** 2
    elif name == 'broad': g = 0.5 * ((tg - 15.) / 7.) ** 2
    elif name == 'shoulder':
        m_ = np.exp(-0.5 * ((tg - 20.) / 4.) ** 2); p = 0.95 * m_ / m_.sum() + 0.05 * (tg <= 3) / (tg <= 3).sum(); return p / p.sum()
    else: raise SystemExit('unknown shape ' + name)
    p = np.exp(-g); return p / p.sum()
if a.shape != 'sim':   # impose p(t) on the coded particles; contaminants keep their share of the node
    tc = tcode[rows]; Nt = np.bincount(tc[tc >= 0], minlength=NG); p = target(a.shape); s = np.min(Nt / np.maximum(p, 1e-12)); nk = np.floor(s * p).astype(int)
    keep = [rng.choice(rows[tc == k], nk[k], replace=False) for k in range(NG)]; ncod = int(nk.sum()); con = rows[tc < 0]
    ncon = int(round(len(con) * ncod / max((tc >= 0).sum(), 1))); keep.append(rng.choice(con, min(ncon, len(con)), replace=False)); rows = np.sort(np.concatenate(keep))
SA = np.asarray(np.load(a.scores + '/scores_pA.npy', mmap_mode='r')[rows][:, use], np.float64); SB = np.asarray(np.load(a.scores + '/scores_pB.npy', mmap_mode='r')[rows][:, use], np.float64)
ok = np.isfinite(SA).all(1) & np.isfinite(SB).all(1); rows, SA, SB = rows[ok], SA[ok], SB[ok]; N = len(rows); eo = EO[rows]; tt = tcode[rows]
ppn = np.asarray(M['ppnoise'])[rows].astype(np.float64); ppn = np.where(np.isfinite(ppn) & (ppn > 0), ppn, np.nanmedian(ppn)); SA /= ppn[:, None]; SB /= ppn[:, None]
X = SA + SB; X -= X.mean(1, keepdims=True)
G = {(e, q): M['gram%s%d' % (q, e)][np.ix_(use, use)] * np.mean(1. / ppn[eo == e]) for e in (0, 1) for q in 'AB'}
Gp = sum((G[(e, 'A')] + G[(e, 'B')]) * (eo == e).sum() for e in (0, 1)) / N
w_, U_ = np.linalg.eigh(np.eye(K) - 1. / K); Bc = U_[:, w_ > 0.5]
lamv, W_ = eigh(Bc.T @ np.cov(X.T) @ Bc, Bc.T @ Gp @ Bc); v1 = Bc @ W_[:, np.argmax(lamv)]; z0 = X @ v1; z0 = (z0 - z0.mean()) / z0.std()
cod = tt >= 0
log('%s shape %s: %d particles (%d coded = %.1f %%), %d leaves; start axis lam-1 %.3f, |corr| with the truth %.3f' % (a.use or 'all', a.shape, N, cod.sum(), 100 * cod.mean(), K,
    lamv.max() - 1, abs(np.corrcoef(z0[cod], tt[cod])[0, 1])))
# ---------------- particle Fourier data at crop c (lsamp2 layout): eo halves of one stack
lamb = 12.2643247 / np.sqrt(a.kv * 1e3 * (1. + a.kv * 1e3 * 0.978466e-6))
def rotm(e1, e2, e3):
    a1, a2, a3 = np.radians(e1), np.radians(e2), np.radians(e3); c1, si1, c2, si2, c3, si3 = np.cos(a1), np.sin(a1), np.cos(a2), np.sin(a2), np.cos(a3), np.sin(a3)
    Z1 = np.zeros((len(e1), 3, 3)); Z1[:, 0, 0] = c1; Z1[:, 0, 1] = -si1; Z1[:, 1, 0] = si1; Z1[:, 1, 1] = c1; Z1[:, 2, 2] = 1
    Y2 = np.zeros((len(e1), 3, 3)); Y2[:, 0, 0] = c2; Y2[:, 0, 2] = si2; Y2[:, 2, 0] = -si2; Y2[:, 2, 2] = c2; Y2[:, 1, 1] = 1
    Z3 = np.zeros((len(e1), 3, 3)); Z3[:, 0, 0] = c3; Z3[:, 0, 1] = -si3; Z3[:, 1, 0] = si3; Z3[:, 1, 1] = c3; Z3[:, 2, 2] = 1
    return Z1 @ Y2 @ Z3
n = int(np.frombuffer(open(a.sim.split(',')[0] + '/all.mrcs', 'rb').read(1024), np.int32)[0]); c = a.crop
R0 = a.mskdiam / 2 / a.smpd; fc = np.fft.fftfreq(c) * c / n; fy2, fx2 = np.meshgrid(fc, fc, indexing='ij'); kk0 = np.stack([fx2, fy2, np.zeros_like(fx2)], -1).reshape(-1, 3) * n
r2d = np.hypot(kk0[:, 0], kk0[:, 1]); pcols = np.where(r2d <= c // 2 - 0.5)[0]; pr = r2d[pcols]; psh = np.round(pr).astype(int)
ppar = np.mod(np.round(kk0[pcols, 0]).astype(int) + np.round(kk0[pcols, 1]).astype(int), 2) == 1   # True = pixel set B
klp = int(round(n * a.smpd / a.lp)); pin = (psh >= 2) & (psh <= klp)
f3 = np.fft.fftfreq(c) * c / n; kx, ky, kz = np.meshgrid(f3, f3, f3, indexing='ij'); sh = np.minimum(np.round(np.sqrt(kx ** 2 + ky ** 2 + kz ** 2) * n).astype(int), c // 2).ravel()
g3 = np.arange(c) - c // 2; r3 = np.sqrt(g3[:, None, None] ** 2 + g3[None, :, None] ** 2 + g3[None, None, :] ** 2) * n / c; msk3 = np.clip((R0 + 4 - r3) / 8., 0, 1).astype(np.float32)
ns_ = c // 2 + 1; khi = c // 2 - 1; KV = np.stack([np.arange(c ** 3) % c, (np.arange(c ** 3) // c) % c, np.arange(c ** 3) // (c * c)], 1); KV = np.where(KV >= c // 2, KV - c, KV).astype(np.float64) / n
FB = (sh >= 2) & (sh <= min(klp, khi))
P_all = np.array([[float(kv(l).get(k, 0)) for k in ('e1', 'e2', 'e3', 'x', 'y', 'dfx', 'dfy', 'angast')] for l in open(a.poses)]); assert len(P_all) == NM
def load_half(sel):
    Fc = np.load('%s_A_c%d.npy' % (a.cache, c), mmap_mode='r'); mrow = rows[sel]; pp = P_all[mrow]; nn = len(sel)
    Wv = np.zeros((nn, len(pcols)), np.complex64); Wd = np.zeros((nn, len(pcols)), np.float32); idx = np.zeros((nn, len(pcols)), np.int32); pw = np.zeros((2, ns_)); pc = np.zeros((2, ns_))
    kf = (np.hypot(fx2, fy2) / a.smpd).ravel()[pcols]; ang = np.arctan2(fy2, fx2).ravel()[pcols]
    for b in range(0, nn, 20000):
        s_ = slice(b, b + 20000); p = pp[s_]; F = np.asarray(Fc[mrow[s_]]).reshape(len(p), -1)[:, pcols]
        df = 0.5 * ((p[:, 5] + p[:, 6])[:, None] + (p[:, 5] - p[:, 6])[:, None] * np.cos(2. * (ang[None] - np.radians(p[:, 7])[:, None]))) * 1e4
        chi = np.pi * lamb * df * kf[None] ** 2 - 0.5 * np.pi * a.cs * 1e7 * lamb ** 3 * kf[None] ** 4; C = (np.sqrt(1. - a.fraca ** 2) * np.sin(chi) + a.fraca * np.cos(chi)).astype(np.float32)
        if a.ctfabs: C = np.abs(C)
        ramp = np.exp(2j * np.pi * (fx2.ravel()[pcols][None] * p[:, 3, None] + fy2.ravel()[pcols][None] * p[:, 4, None])).astype(np.complex64)
        Wv[s_] = C * np.conj(ramp) * F; Wd[s_] = C ** 2; F2 = np.abs(F) ** 2
        for q_ in (0, 1): m_ = ppar == bool(q_); pw[q_] += np.bincount(psh[m_], F2[:, m_].sum(0), ns_); pc[q_] += np.bincount(psh[m_], minlength=ns_) * len(p)
        q = np.round(np.einsum('qj,nij->nqi', kk0[pcols], rotm(p[:, 0], p[:, 1], p[:, 2]))).astype(np.int32); okq = (np.abs(q) <= c // 2 - 1).all(2); q = q % c
        ii = q[..., 2] * c * c + q[..., 1] * c + q[..., 0]; ii[~okq] = c ** 3; idx[s_] = ii
    ph, th = np.radians(pp[:, 0]), np.radians(pp[:, 1]); v = np.c_[np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)]
    i = np.arange(a.nview) + 0.5; z = 1 - 2 * i / a.nview; rr = np.sqrt(np.maximum(1 - z * z, 0)); gg = np.pi * (3 - np.sqrt(5)) * i
    vb = np.argmax(v @ np.c_[rr * np.cos(gg), rr * np.sin(gg), z].T, 1); dfm = 0.5 * (pp[:, 5] + pp[:, 6]); dq = np.searchsorted(np.quantile(dfm, [1 / 3., 2 / 3.]), dfm)
    wpx = (1. / np.maximum(pw / np.maximum(pc, 1), 1e-30))   # per pixel set and shell: 1 / image power (noise at this SNR)
    return dict(Wv=Wv, Wd=Wd, idx=idx, sel=sel, st=vb * 3 + dq, n=nn, wpix=np.where(ppar, wpx[1][psh], wpx[0][psh]).astype(np.float32), ppn=ppn[sel].astype(np.float32))
HS = [load_half(np.where(eo == e)[0]) for e in (0, 1)]
log('map data: halves %d / %d particles, crop %d, %d pixels per image (%d A / %d B), lam band shells 2..%d (%.1f A)' % (HS[0]['n'], HS[1]['n'], c, len(pcols), (~ppar).sum(), ppar.sum(), klp, n * a.smpd / klp))
neq0 = 5000
def rec_raw(H, s_):
    ii = H['idx'][s_].ravel(); m = c ** 3 + 1
    return np.bincount(ii, H['Wv'][s_].real.ravel(), m)[:-1] + 1j * np.bincount(ii, H['Wv'][s_].imag.ravel(), m)[:-1], np.bincount(ii, H['Wd'][s_].ravel(), m)[:-1]
def shell_eps(H):
    den = rec_raw(H, np.arange(H['n']))[1] * (neq0 / H['n']); return np.array([np.median(den[(sh == k) & (den > 0)]) if ((sh == k) & (den > 0)).any() else 1. for k in range(ns_)]) * a.eps
EPS = 0.5 * (shell_eps(HS[0]) + shell_eps(HS[1]))
def rec(H, s_):   # map of the particles s_ (indices into the half), fixed regulariser x N / N0, masked in real space
    num, den = rec_raw(H, s_); V_ = num / (den + EPS[sh] * max(len(s_), 1) / neq0)
    v = np.fft.fftshift(np.real(np.fft.ifftn(V_.reshape(c, c, c)))) * msk3; return np.fft.fftn(np.fft.ifftshift(v)).ravel()
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
CONS = [rec(H, np.arange(H['n'])) for H in HS]
def beta_of(H, M0):   # per-shell scale of the model C beta P M0 to this half's images (M0 = the other half's consensus)
    bn = np.zeros(ns_); bd = np.zeros(ns_); M0x = np.r_[M0, 0]
    for b in range(0, H['n'], 20000):
        s_ = slice(b, b + 20000); P0 = M0x[H['idx'][s_]]
        bn += np.bincount(psh, np.real(P0 * np.conj(H['Wv'][s_])).sum(0), ns_); bd += np.bincount(psh, (H['Wd'][s_] * np.abs(P0) ** 2).sum(0), ns_)
    return np.clip(bn / np.maximum(bd, 1e-30), 0., None)
BETA = [beta_of(HS[e], CONS[1 - e]) for e in (0, 1)]
def lam_of(H, e, D, M0):   # per-particle LS position along D (other half's maps) on pixel sets A and B, and its Fisher information
    Dx = np.r_[D, 0]; M0x = np.r_[M0, 0]; b = BETA[e][psh][None]; out = np.zeros((H['n'], 4))
    for b0 in range(0, H['n'], 20000):
        s_ = slice(b0, b0 + 20000); PD = b * Dx[H['idx'][s_]]; P0 = b * M0x[H['idx'][s_]]; w = H['wpix'][None] / H['ppn'][s_, None] * pin[None]
        r = H['Wv'][s_] - H['Wd'][s_] * P0; num = w * np.real(np.conj(PD) * r); den = w * H['Wd'][s_] * np.abs(PD) ** 2
        nA, dA, nB, dB = num[:, ~ppar].sum(1), den[:, ~ppar].sum(1), num[:, ppar].sum(1), den[:, ppar].sum(1); out[s_] = np.c_[nA / np.maximum(dA, 1e-30), dA, nB / np.maximum(dB, 1e-30), dB]
    return out
def coordinate(z, it):
    qlo, qhi = np.percentile(z, [a.qend, 100 - a.qend]); LH = []
    for e, H in enumerate(HS):
        zz = z[H['sel']]; lo = rec(H, np.where(zz <= qlo)[0]); hi = rec(H, np.where(zz >= qhi)[0]); LH.append((dock(lo, CONS[e]), dock(hi, CONS[e])))
    S = np.zeros((N, 4))
    for e, H in enumerate(HS): S[H['sel']] = lam_of(H, e, LH[1 - e][1] - LH[1 - e][0], CONS[1 - e])
    I = S[:, 1] + S[:, 3]; lam = (S[:, 0] * S[:, 1] + S[:, 2] * S[:, 3]) / np.maximum(I, 1e-30); zsp = (S[:, 0] - S[:, 2]) / np.sqrt(1. / np.maximum(S[:, 1], 1e-30) + 1. / np.maximum(S[:, 3], 1e-30))
    log('  pass %d: |corr(lam, t)| %.3f (coded), split z var %.3f, Fisher sd median %.3f, lam sd %.3f, end-to-end units' % (it, abs(np.corrcoef(lam[cod], tt[cod])[0, 1]), zsp.var(),
        np.median(1. / np.sqrt(np.maximum(I, 1e-30))), lam.std()))
    return lam, I, zsp, LH
z = z0.copy()
for it in range(1, a.nit + 1):
    lam, I, zsp, LH = coordinate(z, it); z = (lam - lam.mean()) / lam.std()
# ---------------- windows: physical composition by NNLS onto the basis, map-level occupancy, calibration lam(t)
BL = [l.split() for l in open(a.basis) if l.strip()]; bidx = list(range(0, len(BL), a.bstep)); bidx = bidx if bidx[-1] == len(BL) - 1 else bidx + [len(BL) - 1]
def crop_four(v):
    nb = v.shape[0]; sb = np.r_[0:c // 2, nb - c // 2:nb]; Fv = np.fft.fftn(np.fft.ifftshift(v)); vc = np.fft.fftshift(np.real(np.fft.ifftn(Fv[np.ix_(sb, sb, sb)]))) * msk3
    return np.fft.fftn(np.fft.ifftshift(vc)).ravel()
BAS = np.array([crop_four(rd(BL[j][1])) for j in bidx]); tb = tg[bidx]; NB = len(bidx); BM = BAS.mean(0)
def unmix(Mw):   # Mw ~ a(shell) * sum_k q_k BAS_k over the band, q >= 0: alternate q (NNLS) and a (per-shell LS)
    Mw = dock(Mw, BM); m_ = FB; q = np.full(NB, 1. / NB)
    for _ in range(6):
        Bq = q @ BAS; aa = np.bincount(sh[m_], np.real(Mw[m_] * np.conj(Bq[m_])), ns_) / np.maximum(np.bincount(sh[m_], np.abs(Bq[m_]) ** 2, ns_), 1e-30)
        Amat = (BAS[:, m_] * aa[sh[m_]][None]).T; q, _ = nnls(np.r_[Amat.real, Amat.imag], np.r_[Mw[m_].real, Mw[m_].imag])
        if q.sum() <= 0: q = np.full(NB, 1. / NB)
    return q / q.sum()
qs = np.percentile(lam, np.linspace(0, 100, a.nw + 1)); wlab = np.clip(np.searchsorted(qs, lam, side='right') - 1, 0, a.nw - 1)
Q = np.zeros((a.nw, NB)); QH = np.zeros((2, a.nw, NB)); LW = np.zeros(a.nw); PW = np.zeros(a.nw); TW = np.zeros(a.nw)
for w in range(a.nw):
    for e, H in enumerate(HS): QH[e, w] = unmix(rec(H, np.where(wlab[H['sel']] == w)[0]))
    Q[w] = QH[:, w].mean(0); LW[w] = lam[wlab == w].mean(); PW[w] = (wlab == w).mean(); m_ = (wlab == w) & cod; TW[w] = tt[m_].mean() if m_.any() else np.nan
if np.corrcoef(LW, Q @ tb)[0, 1] < 0: lam, zsp, LW = -lam, -zsp, -LW   # truth-free orientation: lam increases with the unmixed swing
for w in range(a.nw): log('  window %2d: %6d ptcls, mean lam %+.3f, true t %.1f, unmixed mean t %.1f (halves %.1f / %.1f)' % (w + 1, (wlab == w).sum(), LW[w], TW[w], Q[w] @ tb, QH[0, w] @ tb, QH[1, w] @ tb))
c_all = float(zsp.var()); st = np.zeros(N, int)
for H in HS: st[H['sel']] = H['st']
cs_ = {s: (zsp[st == s].var() if (st == s).sum() > 30 else c_all) for s in np.unique(st)}; c_st = np.array([cs_[s] for s in st])
Et = np.array([lam[tt == k].mean() if (tt == k).any() else np.nan for k in range(NG)])   # oracle E[lam | t]
c_or = float(np.mean(((lam[cod] - Et[tt[cod]]) ** 2) * I[cod]))
log('noise model: split c %.3f (strata 5/50/95 %% %.2f/%.2f/%.2f), oracle c %.3f (residual about E[lam|t] x Fisher); |corr(lam, t)| %.3f' % (c_all, *np.percentile(c_st, [5, 50, 95]), c_or,
    abs(np.corrcoef(lam[cod], tt[cod])[0, 1])))
def expand(qb):   # basis weights -> the t grid (hat functions between basis points)
    p = np.zeros(NG)
    for k in range(NB): p += qb[k] * np.interp(tg, tb, np.eye(NB)[k])
    return p / p.sum()
p_map = expand(PW @ Q)
knots = tb; nk_ = len(knots); Bk = np.array([[np.interp(t_, knots, np.eye(nk_)[j]) for j in range(nk_)] for t_ in tg])   # lam on the grid = Bk @ lam_knots
Aw = np.array([expand(Q[w]) for w in range(a.nw)]) @ Bk   # window w's mean lam = Aw[w] @ lam_knots
Tm = np.tril(np.ones((nk_, nk_)))[:, 1:]   # lam_knots = l0 + Tm @ d with d >= 0 (monotone increasing)
x_, _ = nnls(np.c_[Aw.sum(1), -Aw.sum(1), Aw @ Tm], LW); lknots = (x_[0] - x_[1]) + Tm @ x_[2:]; lam_t = Bk @ lknots
log('calibration lam(t) (truth-free, from the windows): ' + ' '.join('%.0f:%+.3f' % (t_, l_) for t_, l_ in zip(knots, lknots)))
log('oracle E[lam | t] at the knots:                      ' + ' '.join('%.0f:%+.3f' % (t_, Et[int(t_)]) for t_ in knots))
# ---------------- deconvolution on the t grid
def F_(L, p): d = L @ p; return p * (L.T @ (1. / d)) / L.shape[0]
def obj(L, p): return np.log(L @ p).mean()
def npmle(L, p0=None, nit=3000, tol=1e-10):
    p = np.full(L.shape[1], 1. / L.shape[1]) if p0 is None else p0.copy()
    for _ in range(nit):
        p1 = F_(L, p); p2 = F_(L, p1); r = p1 - p; v = p2 - 2. * p1 + p; nv = np.linalg.norm(v)
        if nv < 1e-14: p = p2; break
        al = min(-np.linalg.norm(r) / nv, -1.); pq = np.maximum(p - 2. * al * r + al * al * v, 0.); pq /= pq.sum(); pq = F_(L, pq)
        if not np.isfinite(obj(L, pq)) or obj(L, pq) < obj(L, p2): pq = p2
        done = np.abs(pq - p).max() < tol; p = pq
        if done: break
    return p
def lik(lt, cvec, sel=slice(None)):
    dl = lam[sel, None] - lt[None]; L = np.exp(-0.5 * dl ** 2 * (I[sel] / cvec[sel])[:, None]); return L / np.maximum(L.max(1, keepdims=True), 1e-300)
pt = np.bincount(tt[cod], minlength=NG) / cod.sum()
def sm3(p): return np.convolve(p, np.ones(3), 'same') / np.convolve(np.ones_like(p), np.ones(3), 'same')
def metrics(p):
    G_ = -np.log(np.maximum(sm3(p), 1e-12)); return dict(tv=0.5 * float(np.abs(p - pt).sum()), low=float(p[:15].sum() + 0.5 * p[15]), minor=float(p[:4].sum()), barrier=float(G_[15] - min(G_[5], G_[25])))
EtF = np.where(np.isfinite(Et), Et, np.interp(tg, tg[np.isfinite(Et)], Et[np.isfinite(Et)]))
variants = {'truthfree': (lam_t, c_st), 'oracle': (EtF, np.full(N, c_or)), 'calib_oracle_c_split': (EtF, c_st)}
res = dict(n=int(N), coded=int(cod.sum()), shape=a.shape, c_split=c_all, c_oracle=c_or, truth=metrics(pt), map=metrics(p_map), curves=dict(truth=pt.tolist(), map=p_map.tolist()))
naive_t = np.clip(np.rint(np.interp(lam, lam_t, tg) if np.all(np.diff(lam_t) > 0) else np.interp(lam, np.sort(lam_t), tg)), 0, NG - 1).astype(int); p_naive = np.bincount(naive_t, minlength=NG) / N
res['naive'] = metrics(p_naive); res['curves']['naive'] = p_naive.tolist()
for vn, (lt, cv) in variants.items():
    L = lik(lt, cv); pe = npmle(L); boots = []
    for b in range(a.nboot):
        ib = rng.integers(0, N, N); boots.append(npmle(L[ib], p0=pe, nit=400, tol=1e-8))
    boots = np.array(boots); res[vn] = dict(m=metrics(pe), boot={k: [float(np.percentile([metrics(bb)[k] for bb in boots], q)) for q in (10, 50, 90)] for k in ('tv', 'low', 'minor', 'barrier')})
    res['curves'][vn] = pe.tolist(); res['curves'][vn + '_lo'] = np.percentile(boots, 10, 0).tolist(); res['curves'][vn + '_hi'] = np.percentile(boots, 90, 0).tolist()
log('RESULT %s n %d (coded %d): truth low %.3f minor %.3f barrier %+.2f | map-level tv %.3f low %.3f minor %.3f barrier %+.2f | naive tv %.3f barrier %+.2f minor %.3f' % (
    a.shape, N, cod.sum(), res['truth']['low'], res['truth']['minor'], res['truth']['barrier'], res['map']['tv'], res['map']['low'], res['map']['minor'], res['map']['barrier'],
    res['naive']['tv'], res['naive']['barrier'], res['naive']['minor']))
for vn in variants:
    m_, b_ = res[vn]['m'], res[vn]['boot']; log('RESULT %s %-22s tv %.3f low %.3f minor %.3f [%.3f %.3f] barrier %+.2f [%+.2f %+.2f]' % (a.shape, vn, m_['tv'], m_['low'], m_['minor'], b_['minor'][0], b_['minor'][2],
        m_['barrier'], b_['barrier'][0], b_['barrier'][2]))
json.dump(res, open(a.out + '/e1.json', 'w'), indent=1)
np.savez(a.out + '/e1_particles.npz', rows=rows, lam=lam, I=I, zsplit=zsp, c_st=c_st, t=tt, wlab=wlab, Q=Q, tb=tb, LW=LW, PW=PW, lam_t=lam_t, Et=Et)
try:
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(13, 3.8)); C_ = res['curves']
    for k_, col in (('truth', 'k'), ('map', 'C3'), ('naive', 'C2'), ('truthfree', 'b'), ('oracle', 'C1')):
        ax[0].plot(tg, C_[k_], color=col, lw=2 if k_ == 'truth' else 1.2, label=k_); ax[1].plot(tg, -np.log(np.maximum(sm3(np.array(C_[k_])), 1e-6)) + np.log(sm3(pt).max()), color=col, label=k_)
    ax[0].fill_between(tg, C_['truthfree_lo'], C_['truthfree_hi'], color='b', alpha=0.2); ax[0].legend(fontsize=7); ax[0].set_xlabel('t (A)'); ax[0].set_ylabel('occupancy')
    ax[1].set_ylim(-0.5, 6); ax[1].set_xlabel('t (A)'); ax[1].set_ylabel('G (kT)')
    ax[2].plot(tg, EtF, 'k.-', label='oracle E[lam|t]'); ax[2].plot(tg, lam_t, 'b-', label='windows lam(t)'); ax[2].set_xlabel('t (A)'); ax[2].set_ylabel('lam'); ax[2].legend(fontsize=7)
    fig.suptitle('%s %s' % (a.use or 'all', a.shape)); fig.tight_layout(); fig.savefig(a.out + '/e1.png', dpi=110)
except Exception as e: log('plot skipped: %s' % e)
log('DONE %s' % a.out)
