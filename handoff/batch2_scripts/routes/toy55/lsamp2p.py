#!/usr/bin/env python3
"""lsamp2p.py = lsamp2h.py + --coords: an external 1-D coordinate (coords.npz rows / tA / tB / tfull, e.g. lsamp_poleaxis.py's pole axis)
replaces the score axes wherever the stage takes a coordinate (use with --axes 1). lsamp2h.py = lsamp2.py + held-out sample maps: per fold, the fold's samples (cells of the fold coordinate at the chosen n) reconstructed from
the TEST pixels only (both eo halves pooled), zero-padded from the crop to box 160 and written as OUT/heldout_<fold>/state_NN.mrc: maps
without the selection imprint of the binning, for readouts in the coordinate's band. lsamp2.py (landscape sampling with pixel-split cross-validation): samples of a tree node's reproducible variability, tail to tail at equal
steps, the number of samples set by resolvability.
(1) AXES: centred score vectors of the node's particles against the node's leaf maps (lsamp_score2.py: full score = pixel set A + B) against
    the noise covariance G: generalised eigenvectors C v = lam G v (lam - 1 = map-space signal variance per particle, in noise units), per eo
    half and pooled; reproducible = loadings agree between the halves. NAX leading axes (--naxes; 0 = the leading axis plus every axis with
    lam - 1 >= half of the leading one, at most 3).
(2) SAMPLES on a coordinate t (scores . v, standardised): per axis n samples at equal steps from the QLO to the 100 - QLO percentile (sample k
    centred at q_lo + k (q_hi - q_lo) / (n - 1), end samples open-ended); grids for 2-3 axes; a cell whose half-maps do not reproduce to
    20 A (FSC < MINFSC) merges into its nearest kept neighbour (only such cells merge).
(3) RESOLVABILITY by pixel-split cross-validation (no selection imprint): fold A takes the coordinate from pixel set A only and tests on the
    pixels of set B (independent noise: unmasked images, Fourier pixels split by the parity of kx + ky), fold B the reverse. Per axis the
    end-to-end direction D = high-end minus low-end cell map of the OTHER eo half built from the test pixels (docked); every particle's
    matched-filter projection u on D over the test pixels (all shells, per-shell noise weights), centred per view x defocus stratum and scaled
    so the end cells differ by 1. A neighbour step p -> q is resolved when the mean u increases in both halves and the pooled z exceeds the
    one-sided Bonferroni threshold (alpha ALPHA over the neighbour pairs). n per fold = the largest candidate with every step resolved; the
    node's n = the smaller of the two folds.
(4) FINAL samples on the full coordinate (A + B) at that n. CHECKS: reproducibility (FSC at 20 A, resolution), the fold tests' steps
    (distinguishable neighbours, even spacing = CV of the steps in end-to-end units), end samples vs the tail particles beyond the quantiles
    (fold positions on the test pixels; 0 / 1 = the end samples). With --gt: truth per sample.
usage: lsamp2.py --scores DIR --sim .. --mr .. --cache CP --poses POSES --out OUT [--use r22_] [--naxes 1] [--gt GT --codes swing:200:230]"""
import os, argparse, numpy as np
from scipy.linalg import eigh
from scipy.stats import norm
ap = argparse.ArgumentParser()
ap.add_argument('--scores', required=True); ap.add_argument('--sim', required=True); ap.add_argument('--mr', required=True); ap.add_argument('--cache', required=True)
ap.add_argument('--poses', required=True); ap.add_argument('--out', required=True); ap.add_argument('--use', default=''); ap.add_argument('--exclude', default='')
ap.add_argument('--naxes', type=int, default=1); ap.add_argument('--axes', default=''); ap.add_argument('--ns', default='4,5,6,7,8,10,12,14,16,20,24'); ap.add_argument('--grid_ns', default='2,3,4,5,6,7,8')
ap.add_argument('--qlo', type=float, default=1.); ap.add_argument('--leaf_heldout', type=int, default=0, help='1: only the leaf held-out readout (the tree\'s leaves re-assigned on one pixel set, maps from the other) and exit'); ap.add_argument('--min_per_axis', type=int, default=3, help='grids only while every axis resolves at least this many samples'); ap.add_argument('--tail_tol', type=float, default=0.1, help='grids only while check 3 holds: no tail more than this (end-to-end units) beyond its end sample'); ap.add_argument('--exclude_rows', default='', help='file of master-row flags (1 = drop): e.g. an ORACLE junk removal for synthetic diagnostics'); ap.add_argument('--ends', default='centre', help='centre: end samples centred at the quantiles (open-ended); tail: end samples = the particles beyond the quantiles, interior samples between'); ap.add_argument('--alpha', type=float, default=0.01); ap.add_argument('--minfsc', type=float, default=0.143); ap.add_argument('--fsc_res', type=float, default=20.)
ap.add_argument('--minhalf', type=int, default=150); ap.add_argument('--gt', default=''); ap.add_argument('--codes', default=''); ap.add_argument('--crop', type=int, default=48)
ap.add_argument('--mskdiam', type=float, default=230.); ap.add_argument('--smpd', type=float, default=1.464); ap.add_argument('--kv', type=float, default=300.); ap.add_argument('--cs', type=float, default=2.7)
ap.add_argument('--fraca', type=float, default=0.1); ap.add_argument('--ctfabs', type=int, default=0); ap.add_argument('--eps', type=float, default=0.05); ap.add_argument('--nview', type=int, default=40)
ap.add_argument('--coords', default='', help='external coordinate (coords.npz) instead of the score axes; with --axes 1'); ap.add_argument('--seed', type=int, default=0); ap.add_argument('--tag', default='lsamp2'); ap.add_argument('--axes_only', type=int, default=0)
ap.add_argument('--nit', type=int, default=0, help='score-loop iterations per coordinate (fold A, fold B, full): rescore on the coordinate\'s own pixels against NREF cells\' maps of the other eo half (amplitude-matched per shell), new axes'); ap.add_argument('--nref', type=int, default=12); ap.add_argument('--junk_ssnr', type=float, default=0.35, help='leaves whose per-particle SSNR at 20-10 A relative to the node is below this are junk: dropped before the axes (splitbench16 junk rule); 0 = off'); ap.add_argument('--ppnorm', type=int, default=1, help='1: every particle\'s scores divided by its own noise level (meta ppnoise): the per-particle likelihood'); ap.add_argument('--spacing', default='struct', help="coord: equal steps of the coordinate; struct: equal steps of the held-out structural position E[u | t] (isotonic), ends still at the quantiles"); ap.add_argument('--em_lim', type=float, default=0., help='score band of the loop (default: the scores\' band)')
a = ap.parse_args(); os.makedirs(a.out, exist_ok=True); rng = np.random.default_rng(a.seed)
def log(*s): print(*s, flush=True)
def kv(l): return dict(t.split('=', 1) for t in l.split() if '=' in t)
# ---------------- score space
M = np.load(a.scores + '/meta.npz'); names = [str(x) for x in M['names']]; EO = M['eo']; OWN = M['own']; NM = len(OWN)
use = [j for j, nm in enumerate(names) if (not a.use or any(nm.startswith(p) for p in a.use.split(','))) and nm not in a.exclude.split(',')]
K = len(use); USE0 = list(use); rows = np.where(np.isin(OWN, use) & (EO >= 0))[0]
if a.exclude_rows: drop = np.loadtxt(a.exclude_rows).astype(int) > 0; log('excluded rows (%s): %d of the node\'s %d' % (a.exclude_rows, int(drop[rows].sum()), len(rows))); rows = rows[~drop[rows]]
SA = np.asarray(np.load(a.scores + '/scores_pA.npy', mmap_mode='r')[rows][:, use], np.float64); SB = np.asarray(np.load(a.scores + '/scores_pB.npy', mmap_mode='r')[rows][:, use], np.float64)
ok = np.isfinite(SA).all(1) & np.isfinite(SB).all(1); rows, SA, SB = rows[ok], SA[ok], SB[ok]; N = len(rows); eo = EO[rows]
if a.ppnorm and 'ppnoise' in M.files:
    lam_i = np.asarray(M['ppnoise'])[rows].astype(np.float64); lam_i = np.where(np.isfinite(lam_i) & (lam_i > 0), lam_i, np.nanmedian(lam_i)); SA /= lam_i[:, None]; SB /= lam_i[:, None]
    log('per-particle noise normalisation: noise level 5/50/95 %% = %.2f/%.2f/%.2f' % tuple(np.percentile(lam_i, [5, 50, 95])))
XA = SA - SA.mean(1, keepdims=True); XB = SB - SB.mean(1, keepdims=True); X = XA + XB
log('%s: %d particles of %d leaves (%s); score band to %.1f A' % (a.tag, N, K, ' '.join(names[j] for j in use), float(M['em_lim'])))
G = {(e, q): M['gram%s%d' % (q, e)][np.ix_(use, use)] for e in (0, 1) for q in 'AB'}
if a.ppnorm and 'ppnoise' in M.files: G = {k: v * np.mean(1. / lam_i[eo == k[0]]) for k, v in G.items()}
Gh = [G[(e, 'A')] + G[(e, 'B')] for e in (0, 1)]
w_, U_ = np.linalg.eigh(np.eye(K) - 1. / K); Bc = U_[:, w_ > 0.5]
def geig(Z, Gm):
    lam, W = eigh(Bc.T @ np.cov(Z.T) @ Bc, Bc.T @ Gm @ Bc); o = np.argsort(-lam); return lam[o], Bc @ W[:, o]
lamh, Vh = zip(*[geig(X[eo == e], Gh[e]) for e in (0, 1)]); Gp = (Gh[0] * (eo == 0).sum() + Gh[1] * (eo == 1).sum()) / N; lam, V = geig(X, Gp)
def gcos(u, v): return abs(u @ Gp @ v) / np.sqrt((u @ Gp @ u) * (v @ Gp @ v))
na = min(8, K - 1); agree = [[max(gcos(V[:, j], Vh[e][:, i]) for i in range(min(na + 2, K - 1))) for e in (0, 1)] for j in range(na)]; sig = np.maximum(lam - 1., 0.); share = sig / max(sig.sum(), 1e-12)
log('axes (lam - 1 = signal / noise per particle): axis lam-1 | eo0 / eo1 | loading agreement | share (cum)')
for j in range(na): log('  %4d %7.3f | %7.3f / %7.3f | %.3f / %.3f | %.3f (%.3f)' % (j + 1, lam[j] - 1, lamh[0][j] - 1, lamh[1][j] - 1, agree[j][0], agree[j][1], share[j], share[:j + 1].sum()))
if a.axes: AX = [int(v) - 1 for v in a.axes.split(',')]
elif a.naxes > 0: AX = list(range(a.naxes))
else: AX = [j for j in range(min(3, na)) if j == 0 or (lam[j] - 1 >= 0.5 * (lam[0] - 1) and min(agree[j]) >= 0.9)]
nax = len(AX); log('sampled axes: %s' % ' '.join(str(j + 1) for j in AX))
def coord(Xf): T_ = Xf @ V[:, AX]; return (T_ - T_.mean(0)) / T_.std(0)   # V / AX are module globals (the score loop replaces V)
TZ = {'A': coord(XA), 'B': coord(XB), 'full': coord(X)}
if a.coords:   # the external coordinate stands in for every later coord() call (the junk screen and the axis gate recompute TZ through it)
    C_ = np.load(a.coords); ix_ = np.minimum(np.searchsorted(C_['rows'], rows), len(C_['rows']) - 1); assert np.array_equal(C_['rows'][ix_], rows), 'coords rows differ'
    EXT = {'A': C_['tA'][ix_], 'B': C_['tB'][ix_], 'full': C_['tfull'][ix_]}
    def coord(Xf): return EXT['A'] if Xf is XA else EXT['B'] if Xf is XB else EXT['full']
    TZ = {'A': coord(XA), 'B': coord(XB), 'full': coord(X)}; log('coordinate from %s replaces the score axes' % a.coords)
GT = None
if a.gt:
    g_ = np.loadtxt(a.gt, comments='#', ndmin=2); conf = np.full(NM, -1); clo = np.zeros(NM); conf[g_[:, 0].astype(int) - 1] = g_[:, 1].astype(int); clo[g_[:, 0].astype(int) - 1] = g_[:, 2] if g_.shape[1] > 2 else 0
    GT = dict(conf=conf[rows], closed=clo[rows], codes={})
    for cd in [w for w in a.codes.split(',') if w]: nm, lo, hi = cd.split(':'); cf = GT['conf']; GT['codes'][nm] = np.where((cf >= int(lo)) & (cf <= int(hi)), cf - int(lo), np.nan)
    def nmi(x, y):
        _, xi = np.unique(x, return_inverse=True); _, yi = np.unique(y, return_inverse=True); J = np.zeros((xi.max() + 1, yi.max() + 1)); np.add.at(J, (xi, yi), 1); J /= J.sum()
        px, py = J.sum(1), J.sum(0); nz = J > 0; return (J[nz] * np.log(J[nz] / np.outer(px, py)[nz])).sum() / max(np.sqrt((px * np.log(px)).sum() * (py * np.log(py)).sum()), 1e-12)
    for nm, x in GT['codes'].items():
        m = np.isfinite(x)
        if m.sum() > 100: log('TRUTH |corr(t, %s)| per sampled axis: full %s | pixel set A %s | B %s' % (nm, ' '.join('%.3f' % abs(np.corrcoef(TZ['full'][m, j], x[m])[0, 1]) for j in range(nax)),
                              ' '.join('%.3f' % abs(np.corrcoef(TZ['A'][m, j], x[m])[0, 1]) for j in range(nax)), ' '.join('%.3f' % abs(np.corrcoef(TZ['B'][m, j], x[m])[0, 1]) for j in range(nax))))
    for nm, x in GT['codes'].items():
        m = np.isfinite(x)
        if m.sum() > 100:
            Tall = X @ V; bd = []
            for mm in (1, 2, 3, min(7, K - 1)):
                A_ = np.c_[Tall[m, :mm], np.ones(m.sum())]; cf = np.linalg.lstsq(A_, x[m], rcond=None)[0]; bd.append('%d axes %.3f' % (mm, np.corrcoef(A_ @ cf, x[m])[0, 1]))
            log('TRUTH best linear direction of %s within the leading axes: %s' % (nm, ', '.join(bd)))
if a.axes_only: raise SystemExit
# ---------------- particle data (crop 48, the scoring cache): halves = half stacks or the eo partition
lamb = 12.2643247 / np.sqrt(a.kv * 1e3 * (1. + a.kv * 1e3 * 0.978466e-6))
def rotm(e1, e2, e3):
    a1, a2, a3 = np.radians(e1), np.radians(e2), np.radians(e3); c1, si1, c2, si2, c3, si3 = np.cos(a1), np.sin(a1), np.cos(a2), np.sin(a2), np.cos(a3), np.sin(a3)
    Z1 = np.zeros((len(e1), 3, 3)); Z1[:, 0, 0] = c1; Z1[:, 0, 1] = -si1; Z1[:, 1, 0] = si1; Z1[:, 1, 1] = c1; Z1[:, 2, 2] = 1
    Y2 = np.zeros((len(e1), 3, 3)); Y2[:, 0, 0] = c2; Y2[:, 0, 2] = si2; Y2[:, 2, 0] = -si2; Y2[:, 2, 2] = c2; Y2[:, 1, 1] = 1
    Z3 = np.zeros((len(e1), 3, 3)); Z3[:, 0, 0] = c3; Z3[:, 0, 1] = -si3; Z3[:, 1, 0] = si3; Z3[:, 1, 1] = c3; Z3[:, 2, 2] = 1
    return Z1 @ Y2 @ Z3
SIMS = a.sim.split(','); full = a.mr == 'full'; n = int(np.frombuffer(open(SIMS[0] + '/all.mrcs', 'rb').read(1024), np.int32)[0]); c = a.crop
R0 = a.mskdiam / 2 / a.smpd; fc = np.fft.fftfreq(c) * c / n; fy2, fx2 = np.meshgrid(fc, fc, indexing='ij'); kk0 = np.stack([fx2, fy2, np.zeros_like(fx2)], -1).reshape(-1, 3) * n
r2d = np.hypot(kk0[:, 0], kk0[:, 1]); pcols = np.where(r2d <= c // 2 - 0.5)[0]; pr = r2d[pcols]; psh = np.round(pr).astype(int)
ppar = np.mod(np.round(kk0[pcols, 0]).astype(int) + np.round(kk0[pcols, 1]).astype(int), 2) == 1   # True = pixel set B
f3 = np.fft.fftfreq(c) * c / n; kx, ky, kz = np.meshgrid(f3, f3, f3, indexing='ij'); sh = np.minimum(np.round(np.sqrt(kx ** 2 + ky ** 2 + kz ** 2) * n).astype(int), c // 2).ravel()
g3 = np.arange(c) - c // 2; r3 = np.sqrt(g3[:, None, None] ** 2 + g3[None, :, None] ** 2 + g3[None, None, :] ** 2) * n / c; msk3 = np.clip((R0 + 4 - r3) / 8., 0, 1).astype(np.float32)
ns_ = c // 2 + 1; khi = c // 2 - 1; kres = int(round(n * a.smpd / a.fsc_res)); KV = np.stack([np.arange(c ** 3) % c, (np.arange(c ** 3) // c) % c, np.arange(c ** 3) // (c * c)], 1); KV = np.where(KV >= c // 2, KV - c, KV).astype(np.float64) / n
FB = (sh >= 2) & (sh <= khi)
P_all = np.array([[float(kv(l).get(k, 0)) for k in ('e1', 'e2', 'e3', 'x', 'y', 'dfx', 'dfy', 'angast')] for l in open(a.poses)]); assert len(P_all) == NM
inrow = np.full(NM, -1); inrow[rows] = np.arange(N)
def load_half(tag, crow, mrow):
    Fc = np.load('%s_%s_c%d.npy' % (a.cache, tag, c), mmap_mode='r'); pp = P_all[mrow]; nn = len(crow)
    Wv = np.zeros((nn, len(pcols)), np.complex64); Wd = np.zeros((nn, len(pcols)), np.float32); idx = np.zeros((nn, len(pcols)), np.int32); pw = np.zeros((2, ns_)); pc = np.zeros((2, ns_))
    kf = (np.hypot(fx2, fy2) / a.smpd).ravel()[pcols]; ang = np.arctan2(fy2, fx2).ravel()[pcols]
    for b in range(0, nn, 20000):
        s_ = slice(b, b + 20000); p = pp[s_]; F = np.asarray(Fc[crow[s_]]).reshape(len(p), -1)[:, pcols]
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
    return dict(Wv=Wv, Wd=Wd, idx=idx, ti=inrow[mrow], st=vb * 3 + dq, n=nn, w=(1. / np.maximum(pw / np.maximum(pc, 1), 1e-30)).astype(np.float32))
if full: HS = [load_half('A', rows[eo == e], rows[eo == e]) for e in (0, 1)]
else:
    HS = []
    for h, tg in enumerate('AB'):
        mr = np.loadtxt(a.mr.split(',')[h]).astype(int) - 1; sel = np.where(inrow[mr] >= 0)[0]; HS.append(load_half(tg, sel, mr[sel]))
log('map data: halves %d / %d particles (%s), crop %d, %d pixels per image (%d set A / %d set B)' % (HS[0]['n'], HS[1]['n'], 'eo partition' if full else 'half stacks', c, len(pcols), (~ppar).sum(), ppar.sum()))
COLS = {'all': np.arange(len(pcols)), 'A': np.where(~ppar)[0], 'B': np.where(ppar)[0]}; neq0 = 5000
def rec_raw(H, sel, cs):
    ii = H['idx'][np.ix_(sel, cs)].ravel(); m = c ** 3 + 1; Wv_ = H['Wv'][np.ix_(sel, cs)]
    return np.bincount(ii, Wv_.real.ravel(), m)[:-1] + 1j * np.bincount(ii, Wv_.imag.ravel(), m)[:-1], np.bincount(ii, H['Wd'][np.ix_(sel, cs)].ravel(), m)[:-1]
def shell_eps(H):
    den = rec_raw(H, np.arange(H['n']), COLS['all'])[1] * (neq0 / H['n']); return np.array([np.median(den[(sh == k) & (den > 0)]) if ((sh == k) & (den > 0)).any() else 1. for k in range(ns_)]) * a.eps
EPS = 0.5 * (shell_eps(HS[0]) + shell_eps(HS[1]))
def rec(H, sel, cs='all'):   # map of particles sel (indices into the half) from the pixels of set cs, fixed regulariser x N / N0, masked in real space
    cs_ = COLS[cs]; num, den = rec_raw(H, sel, cs_); f_ = len(cs_) / len(pcols); V_ = num / (den + EPS[sh] * f_ * max(len(sel), 1) / neq0)
    v = np.fft.fftshift(np.real(np.fft.ifftn(V_.reshape(c, c, c)))) * msk3; return np.fft.fftn(np.fft.ifftshift(v)).ravel()
CONS = {cs: [rec(H, np.arange(H['n']), cs) for H in HS] for cs in ('all', 'A', 'B')}
def fsc(X1, X2): return np.bincount(sh, np.real(X1 * np.conj(X2)), ns_) / np.maximum(np.sqrt(np.bincount(sh, np.abs(X1) ** 2, ns_) * np.bincount(sh, np.abs(X2) ** 2, ns_)), 1e-30)
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
def write_mrc(fn, v, apix):
    v = np.ascontiguousarray(v, np.float32); nz, ny, nx = v.shape; h = np.zeros(256, np.int32); h[0:3] = [nx, ny, nz]; h[3] = 2; h[7:10] = [nx, ny, nz]
    hf = h.view(np.float32); hf[10:13] = [nx * apix, ny * apix, nz * apix]; hf[13:16] = 90.; h[16:19] = [1, 2, 3]; hf[19], hf[20], hf[21] = v.min(), v.max(), v.mean()
    h[52] = np.frombuffer(b'MAP ', np.int32)[0]; h[53] = np.frombuffer(bytes([0x44, 0x44, 0, 0]), np.int32)[0]
    with open(fn, 'wb') as f: f.write(h.tobytes()); f.write(v.tobytes())
NB = 160
def padmap(M_):
    Fc = M_.reshape(c, c, c); Fb = np.zeros((NB, NB, NB), np.complex128); i_ = np.r_[0:c // 2, NB - c // 2:NB]; Fb[np.ix_(i_, i_, i_)] = Fc
    return np.fft.fftshift(np.real(np.fft.ifftn(Fb))) * (NB / c) ** 3
if a.leaf_heldout:   # the tree's own leaves, re-assigned on one pixel set (argmax of the cross-half scores over the node's leaves) and read on the other
    orig = np.searchsorted(use, OWN[rows]); lab_ = {'A': SA.argmax(1), 'B': SB.argmax(1)}
    log('leaf held-out mode: %d particles of %d leaves; re-assignment agrees with the tree\'s labels for %.1f %% (set A) / %.1f %% (set B); A vs B %.1f %%'
        % (N, K, 100 * np.mean(lab_['A'] == orig), 100 * np.mean(lab_['B'] == orig), 100 * np.mean(lab_['A'] == lab_['B'])))
    for k in range(K):
        mo = orig == k; log('  leaf %d = %s: tree %d particles; re-assigned A %d (%.0f %% of them from the leaf), B %d (%.0f %%)' % (k + 1, names[use[k]], mo.sum(),
            (lab_['A'] == k).sum(), 100 * np.mean(mo[lab_['A'] == k]) if (lab_['A'] == k).any() else 0, (lab_['B'] == k).sum(), 100 * np.mean(mo[lab_['B'] == k]) if (lab_['B'] == k).any() else 0))
    for sub, labv, pix in (('heldout_leafA', lab_['A'], 'B'), ('own_leafA', lab_['A'], 'A'), ('heldout_leafB', lab_['B'], 'A'), ('own_leafB', lab_['B'], 'B'), ('tree_all', orig, 'all')):
        D = os.path.join(a.out, sub); os.makedirs(D, exist_ok=True); lo_ = np.zeros(NM, int); lo_[rows] = labv + 1; np.savetxt(os.path.join(D, 'labels.txt'), lo_, fmt='%d')
        for k in range(K):
            sl = [np.where(labv[H['ti']] == k)[0] for H in HS]
            if min(len(sl[0]), len(sl[1])) < 20: continue
            write_mrc(os.path.join(D, 'state_%02d.mrc' % (k + 1)), padmap(0.5 * (rec(HS[0], sl[0], pix) + rec(HS[1], sl[1], pix))), n * a.smpd / NB)
        log('%s: leaf labels from %s, maps from pixels %s -> %s' % (sub, 'the tree' if sub == 'tree_all' else 'pixel set ' + sub[-1], pix, D))
    log('DONE'); raise SystemExit
# ---------------- junk screen (splitbench16 rule): leaves without the node's high-resolution signal are set aside before the axes
if a.junk_ssnr > 0:
    def relssnr(sA, sB, cA, cB, nc):
        b_ = (sh > kres) & (sh <= khi); out = []
        for X1, X2, nn in ((rec(HS[0], sA), rec(HS[1], sB), 0.5 * (len(sA) + len(sB))), (cA, cB, nc)):
            f = np.real(np.vdot(X2[b_], X1[b_])) / max(np.sqrt(np.real(np.vdot(X1[b_], X1[b_])) * np.real(np.vdot(X2[b_], X2[b_]))), 1e-30); out.append(max(f, 0.) / max(1 - f, 1e-3) / nn)
        return out[0] / max(out[1], 1e-12)
    own_h = [OWN[rows][H['ti']] for H in HS]; leafs = sorted(set(OWN[rows])); nc = 0.5 * (HS[0]['n'] + HS[1]['n'])
    rs = {u: relssnr(np.where(own_h[0] == u)[0], np.where(own_h[1] == u)[0], CONS['all'][0], CONS['all'][1], nc) for u in leafs}; junk = [u for u in leafs if rs[u] < a.junk_ssnr]
    log('junk screen (per-particle SSNR at %.0f-%.0f A relative to the node, threshold %.2f): %s' % (a.fsc_res, n * a.smpd / khi, a.junk_ssnr, ' '.join('%s:%.2f%s' % (names[u], rs[u], '*' if u in junk else '') for u in leafs)))
    if junk and len(junk) < len(leafs):
        keep = ~np.isin(OWN[rows], junk); newi = np.full(N, -1); newi[keep] = np.arange(keep.sum()); use = [u for u in use if u not in junk]; K = len(use);
        colk = [j for j, u in enumerate(USE0) if u not in junk]; XA, XB = XA[keep][:, colk], XB[keep][:, colk]; XA -= XA.mean(1, keepdims=True); XB -= XB.mean(1, keepdims=True); X = XA + XB
        rows = rows[keep]; eo = eo[keep]; N = len(rows); G = {k: v[np.ix_(colk, colk)] for k, v in G.items()}; Gh = [G[(e, 'A')] + G[(e, 'B')] for e in (0, 1)]
        if GT is not None: GT['conf'], GT['closed'] = GT['conf'][keep], GT['closed'][keep]; GT['codes'] = {k: v[keep] for k, v in GT['codes'].items()}
        for H in HS:
            m_ = keep[H['ti']]
            for k_ in ('Wv', 'Wd', 'idx', 'st'): H[k_] = H[k_][m_]
            H['ti'] = newi[H['ti'][m_]]; H['n'] = int(m_.sum())
        K = len(colk); w_, U_ = np.linalg.eigh(np.eye(K) - 1. / K); Bc = U_[:, w_ > 0.5]
        lamh, Vh = zip(*[geig(X[eo == e], Gh[e]) for e in (0, 1)]); Gp = (Gh[0] * (eo == 0).sum() + Gh[1] * (eo == 1).sum()) / N; lam, V = geig(X, Gp); sig = np.maximum(lam - 1., 0.); share = sig / max(sig.sum(), 1e-12)
        na = min(8, K - 1); agree = [[max(gcos(V[:, j], Vh[e][:, i]) for i in range(min(na + 2, K - 1))) for e in (0, 1)] for j in range(na)]
        log('after the junk screen: %d particles of %d leaves; axes lam - 1 %s (shares %s)' % (N, K, ' '.join('%.2f' % v for v in lam[:6] - 1), ' '.join('%.2f' % v for v in share[:6])))
        if not a.axes: AX = list(range(a.naxes)) if a.naxes > 0 else [j for j in range(min(3, na)) if j == 0 or (lam[j] - 1 >= 0.5 * (lam[0] - 1) and min(agree[j]) >= 0.9)]
        nax = len(AX); TZ = {'A': coord(XA), 'B': coord(XB), 'full': coord(X)}; CONS = {cs: [rec(H, np.arange(H['n']), cs) for H in HS] for cs in ('all', 'A', 'B')}
        log('sampled axes after the junk screen: %s' % ' '.join(str(j + 1) for j in AX))
        if GT is not None:
            for nm, x in GT['codes'].items():
                mk = np.isfinite(x)
                if mk.sum() > 100: log('TRUTH |corr(t, %s)| after the junk screen: %s' % (nm, ' '.join('%.3f' % abs(np.corrcoef(TZ['full'][mk, j], x[mk])[0, 1]) for j in range(nax))))
# ---------------- cells
CENT = {}   # per coordinate set and axis: function n -> sample centres in t (structural spacing); absent = equal coordinate steps
def make_cells(Tz, nper, cs=None):
    idx = []
    for j in range(nax):
        t = Tz[:, j]; lo, hi = np.percentile(t, [a.qlo, 100 - a.qlo])
        if a.ends == 'tail' and nper[j] >= 3:   # ends = the tails beyond the quantiles; interior samples between them
            if cs is not None and (cs, j) in CENT: ce = CENT[(cs, j)](nper[j])[1:-1]; mid = 0.5 * (ce[1:] + ce[:-1]); ii = 1 + np.searchsorted(mid, t)
            else: ii = 1 + np.clip(((t - lo) / (hi - lo) * (nper[j] - 2)).astype(int), 0, nper[j] - 3)
            idx.append(np.where(t < lo, 0, np.where(t > hi, nper[j] - 1, ii)))
        elif cs is not None and (cs, j) in CENT: ce = CENT[(cs, j)](nper[j]); mid = 0.5 * (ce[1:] + ce[:-1]); idx.append(np.searchsorted(mid, t))
        else: w = (hi - lo) / (nper[j] - 1); idx.append(np.clip(np.round((t - lo) / w), 0, nper[j] - 1).astype(int))
    return np.ravel_multi_index(idx, nper), np.array(idx)
def pava(y, w):   # isotonic (non-decreasing) least squares
    v, wt, cnt = [], [], []
    for yi, wi in zip(y, w):
        v.append(yi); wt.append(wi); cnt.append(1)
        while len(v) > 1 and v[-2] > v[-1]: w_ = wt[-2] + wt[-1]; v[-2] = (v[-2] * wt[-2] + v[-1] * wt[-1]) / w_; wt[-2] = w_; cnt[-2] += cnt[-1]; v.pop(); wt.pop(); cnt.pop()
    return np.repeat(v, cnt)
def struct_centres(Tz, cs, test):   # E[u | t] per axis from 4 equal-step cells' end-to-end direction on the test pixels; centres at equal steps of it
    cell, idx = make_cells(Tz, [4] * nax); U = ustat(cell, idx, [4] * nax, test, Tz); out = {}
    for j in range(nax):
        t = Tz[:, j]; u = np.zeros(N)
        for h, H in enumerate(HS): u[H['ti']] = U[j][h]
        qb = np.percentile(t, np.linspace(0, 100, 41)); b = np.clip(np.searchsorted(qb, t, side='right') - 1, 0, 39); tm = np.array([t[b == k].mean() for k in range(40)]); um = np.array([u[b == k].mean() for k in range(40)])
        f = pava(um, np.bincount(b, minlength=40).astype(float)) + 1e-6 * np.arange(40); lo, hi = np.percentile(t, [a.qlo, 100 - a.qlo]); flo, fhi = np.interp([lo, hi], tm, f)
        if a.ends == 'tail': out[j] = (lambda nn, tm=tm, f=f, lo=lo, hi=hi, flo=flo, fhi=fhi: np.r_[lo, np.interp(flo + (fhi - flo) * (np.arange(1, nn - 1) - 0.5) / max(nn - 2, 1), f, tm), hi])
        else: out[j] = (lambda nn, tm=tm, f=f, lo=lo, hi=hi, flo=flo, fhi=fhi: np.r_[lo, np.interp(flo + (fhi - flo) * np.arange(1, nn - 1) / (nn - 1), f, tm), hi])
        log('  structural position E[u | t] along axis %d (coordinate %s, test pixels %s), 10 quantile steps: %s' % (AX[j] + 1, cs, test, ' '.join('%.2f' % v for v in np.interp(np.percentile(t, np.linspace(5, 95, 10)), tm, f))))
    return out
def hl(cell): return [cell[H['ti']] for H in HS]
def cell_fsc(cell, g, cs):
    s = [np.where(l == g)[0] for l in hl(cell)]
    if min(len(s[0]), len(s[1])) < a.minhalf: return -1.
    f = fsc(rec(HS[0], s[0], cs), rec(HS[1], s[1], cs)); return float(f[max(kres - 1, 2):kres + 2].mean())
def merge_sparse(Tz, cell, cs):
    cell = cell.copy(); merged = {}; fs = {g: cell_fsc(cell, g, cs) for g in np.unique(cell)}
    while True:
        gs = list(fs); bad = [g for g in gs if fs[g] < a.minfsc]
        if not bad or len(gs) < 2: return cell, fs, merged
        g = min(bad, key=lambda q: (cell == q).sum()); keep = [q for q in gs if q != g]; cen = {q: Tz[cell == q].mean(0) for q in gs}
        tgt = min(keep, key=lambda q: ((cen[q] - cen[g]) ** 2).sum()); cell[cell == g] = tgt; merged[g] = tgt; del fs[g]; fs[tgt] = cell_fsc(cell, tgt, cs)
def neighbours(cell, idx, nper):   # pairs of kept cells adjacent on the grid, oriented low -> high along axis j
    occ = {}
    for g in np.unique(cell):
        for q in np.unique(np.ravel_multi_index(tuple(idx[:, cell == g]), nper)): occ[int(q)] = g
    pr = set()
    for q, g in occ.items():
        iv = np.array(np.unravel_index(q, nper))
        for j in range(nax):
            iw = iv.copy(); iw[j] += 1
            if iw[j] < nper[j]:
                q2 = int(np.ravel_multi_index(iw, nper))
                if q2 in occ and occ[q2] != g: pr.add((g, occ[q2], j))
    return sorted(pr)
def ustat(cell, idx, nper, test, Tz):   # per axis: the test-pixel projection u of every particle on the other half's end-to-end direction, centred per stratum, end layers 1 apart
    U = []
    for j in range(nax):
        lo_m, hi_m = idx[j] == 0, idx[j] == nper[j] - 1
        D = []
        for h, H in enumerate(HS):
            ti = H['ti']; E = [dock(rec(H, np.where(m_[ti])[0], test), CONS[test][1 - h]) for m_ in (lo_m, hi_m)]; D.append(E[1] - E[0])
        uh = []
        for h, H in enumerate(HS):
            Dp = np.r_[D[1 - h], 0.].astype(np.complex64); cs_ = COLS[test]; w_ = H['w'][1 if test == 'B' else 0][psh[cs_]] * (psh[cs_] >= 2); u = np.zeros(H['n'])
            for b in range(0, H['n'], 20000):
                s_ = slice(b, b + 20000); u[s_] = (np.real(Dp[H['idx'][s_][:, cs_]] * np.conj(H['Wv'][s_][:, cs_])) * w_).sum(1)
            for k in np.unique(H['st']): m_ = H['st'] == k; u[m_] -= u[m_].mean()
            ti = H['ti']; sc = u[hi_m[ti]].mean() - u[lo_m[ti]].mean(); uh.append(u / sc if abs(sc) > 1e-30 else u)
        U.append(uh)
    return U
def steps(cell, prs, U):   # per neighbour pair: pooled step (end-to-end units), z, both halves' signs
    out = {}
    for p, q, j in prs:
        ds, vs = [], []
        for h, H in enumerate(HS):
            lab = cell[H['ti']]; up, uq = U[j][h][lab == p], U[j][h][lab == q]
            if len(up) < 2 or len(uq) < 2: ds.append(0.); vs.append(1e30); continue
            ds.append(uq.mean() - up.mean()); vs.append(up.var() / len(up) + uq.var() / len(uq))
        wts = 1. / np.array(vs); d = float((wts * np.array(ds)).sum() / wts.sum()); out[(p, q, j)] = (d, d * np.sqrt(wts.sum()), ds[0], ds[1])
    return out
def evaluate(Tz, nper, coordset, test):
    cell0, idx = make_cells(Tz, nper, coordset); cell, fs, merged = merge_sparse(Tz, cell0, test); prs = neighbours(cell, idx, nper)
    if len(np.unique(cell)) < 2 or not prs: return dict(ok=False, cell=cell)
    U = ustat(cell, idx, nper, test, Tz); st = steps(cell, prs, U); zc = norm.ppf(1 - a.alpha / len(prs))
    ok = all(v[1] > zc and v[2] > 0 and v[3] > 0 for v in st.values())
    log('  fold %s->%s n %s: %d samples (%d merged), steps %s; z %s (threshold %.2f) -> %s' % (coordset, test, 'x'.join(map(str, nper)), len(np.unique(cell)), len(merged),
        ' '.join('%.3f' % st[k][0] for k in prs), ' '.join('%.1f' % st[k][1] for k in prs), zc, 'ALL RESOLVED' if ok else '%d not resolved' % sum(not (v[1] > zc and v[2] > 0 and v[3] > 0) for v in st.values())))
    return dict(ok=ok, cell=cell, idx=idx, fs=fs, merged=merged, prs=prs, st=st, U=U, nper=nper, zc=zc)
kbl = int(round(n * a.smpd / (a.em_lim if a.em_lim > 0 else float(M['em_lim'])))); BCOL = {cs: np.where((psh >= 2) & (psh <= kbl) & (np.ones_like(ppar) if cs == 'all' else (ppar == (cs == 'B'))))[0] for cs in ('all', 'A', 'B')}
def kmeans(Z, k, reps=8):
    best = None
    for r in range(reps):
        C = Z[rng.choice(len(Z), k, replace=False)]
        for it in range(60):
            lab = ((Z[:, None] - C[None]) ** 2).sum(-1).argmin(1); C2 = np.array([Z[lab == j].mean(0) if (lab == j).any() else C[j] for j in range(k)])
            if np.allclose(C2, C): break
            C = C2
        sse = ((Z - C[lab]) ** 2).sum(); best = (sse, lab) if best is None or sse < best[0] else best
    return best[1]
def unmix(Zs):   # score-loop pass (lrec_pc2 / rflex_coord): NREF k-means cells in the leading whitened subspace -> their maps of the other eo half (amplitude matched
    # per shell to the mean) -> every particle rescored on pixel set A and on set B -> new axes from the sum; fold coordinates use the same loadings
    cell = kmeans(Zs, a.nref); gs_ = list(np.unique(cell)); lab_ = hl(cell); MP = [[rec(HS[h], np.where(lab_[h] == g)[0]) for g in gs_] for h in range(2)]
    wsz = np.array([(cell == g).sum() for g in gs_], float); mA = sum(w * m for w, m in zip(wsz, MP[0])) / wsz.sum(); mB = sum(w * m for w, m in zip(wsz, MP[1])) / wsz.sum()
    pm = np.maximum(np.bincount(sh, np.real(mA * np.conj(mB)), ns_), 1e-30)
    for k in range(len(gs_)):
        sc = np.sqrt(np.maximum(np.bincount(sh, np.real(MP[0][k] * np.conj(MP[1][k])), ns_), 1e-30) / pm); sc[:2] = 1.; sc = np.clip(sc, 0.25, 4.)
        MP[0][k] = np.r_[MP[0][k] / sc[sh], 0.].astype(np.complex64); MP[1][k] = np.r_[MP[1][k] / sc[sh], 0.].astype(np.complex64)
    Kc = len(gs_); S2 = {q: np.zeros((N, Kc)) for q in 'AB'}; G2 = {(e, q): np.zeros((Kc, Kc)) for e in (0, 1) for q in 'AB'}
    for h, H in enumerate(HS):
        ref = MP[1 - h]
        for q in 'AB':
            cs_ = BCOL[q]; shp = psh[cs_]; w = H['w'][0 if q == 'A' else 1][shp]; nb = 0
            for b in range(0, H['n'], 20000):
                s_ = slice(b, b + 20000); Wv_ = H['Wv'][s_][:, cs_]; Wd_ = H['Wd'][s_][:, cs_]; ix_ = H['idx'][s_][:, cs_]; ng = min(1500, Wv_.shape[0]); PG = np.zeros((ng, len(cs_), Kc), np.complex64)
                for k in range(Kc): P_ = ref[k][ix_]; S2[q][H['ti'][s_], k] = ((np.real(P_ * np.conj(Wv_)) - 0.5 * Wd_ * np.abs(P_) ** 2) * w).sum(1); PG[:, :, k] = P_[:ng]
                X_ = PG * np.sqrt(w * Wd_[:ng])[:, :, None]; G2[(h, q)] += (np.einsum('npk,npl->kl', X_.real, X_.real) + np.einsum('npk,npl->kl', X_.imag, X_.imag)) / ng; nb += 1
            G2[(h, q)] /= max(nb, 1)
    XA2 = S2['A'] - S2['A'].mean(1, keepdims=True); XB2 = S2['B'] - S2['B'].mean(1, keepdims=True); X2 = XA2 + XB2
    nh = [HS[h]['n'] for h in range(2)]; Gp2 = sum((G2[(h, 'A')] + G2[(h, 'B')]) * nh[h] for h in range(2)) / sum(nh)
    w2, U2 = np.linalg.eigh(np.eye(Kc) - 1. / Kc); Bc2 = U2[:, w2 > 0.5]; lam2, W2 = eigh(Bc2.T @ np.cov(X2.T) @ Bc2, Bc2.T @ Gp2 @ Bc2); o = np.argsort(-lam2); V2 = Bc2 @ W2[:, o]
    return XA2, XB2, X2, V2, lam2[o]
for it in range(a.nit):
    m_ = max(2, nax); T0 = X @ V[:, :m_]; Zs = (T0 - T0.mean(0)) / T0.std(0)
    XA, XB, X, V, lam = unmix(Zs); TZ = {'A': coord(XA), 'B': coord(XB), 'full': coord(X)}
    msg = 'score loop %d (k-means %d in the leading %d-D subspace): new axes lam - 1 %s' % (it + 1, a.nref, m_, ' '.join('%.2f' % v for v in lam[:4] - 1))
    if GT is not None:
        for nm, x in GT['codes'].items():
            mk = np.isfinite(x)
            if mk.sum() > 100:
                Tall = X @ V[:, :4]; msg += '; |corr(axis 1..4, %s)| %s; sampled axes: full %s | A %s | B %s' % (nm, ' '.join('%.3f' % abs(np.corrcoef(Tall[mk, j], x[mk])[0, 1]) for j in range(min(4, Tall.shape[1]))),
                    ' '.join('%.3f' % abs(np.corrcoef(TZ['full'][mk, j], x[mk])[0, 1]) for j in range(nax)), ' '.join('%.3f' % abs(np.corrcoef(TZ['A'][mk, j], x[mk])[0, 1]) for j in range(nax)), ' '.join('%.3f' % abs(np.corrcoef(TZ['B'][mk, j], x[mk])[0, 1]) for j in range(nax)))
    log(msg)
AX_all = list(AX)
for m_ax in range(len(AX_all), 0, -1):   # a grid only while every axis still resolves >= MIN_PER_AXIS samples; otherwise drop the weakest axis
    AX = AX_all[:m_ax]; nax = m_ax; TZ = {'A': coord(XA), 'B': coord(XB), 'full': coord(X)}; CENT.clear()
    if a.spacing == 'struct':
        for cs_, te_ in (('A', 'B'), ('B', 'A')):
            for j, fn in struct_centres(TZ[cs_], cs_, te_).items(): CENT[(cs_, j)] = fn
        def full_centres(j):   # the final coordinate takes the folds' centre quantile levels
            def fn(nn):
                lv = np.mean([np.searchsorted(np.sort(TZ[cs_][:, j]), CENT[(cs_, j)](nn)) / N for cs_ in ('A', 'B')], 0); return np.percentile(TZ['full'][:, j], 100 * np.clip(lv, 0, 1))
            return fn
        for j in range(nax): CENT[('full', j)] = full_centres(j)
    cands = [int(v) for v in (a.ns if nax == 1 else a.grid_ns).split(',')]; best = {}; RES = {}
    for coordset, test in (('A', 'B'), ('B', 'A')):
        b_ = None; nf = 0
        log('resolvability, fold %s -> %s (coordinate from pixel set %s, tests on set %s; axes %s):' % (coordset, test, coordset, test, '+'.join(str(j + 1) for j in AX)))
        for nc in cands:
            r_ = evaluate(TZ[coordset], [nc] * nax, coordset, test); RES[(coordset, nc)] = r_
            if r_['ok']: b_ = nc; nf = 0
            else: nf += 1
            if b_ is not None and nf >= 2: break
        best[coordset] = b_
    log('resolvable samples per axis: fold A %s, fold B %s' % (best['A'], best['B']))
    if nax == 1: break
    if best['A'] is None or best['B'] is None or min(best['A'], best['B']) < a.min_per_axis:
        log('a grid on %d axes resolves fewer than %d samples per axis: dropping axis %d' % (nax, a.min_per_axis, AX[-1] + 1)); continue
    nb_ = min(best['A'], best['B']); short = 0.   # check 3 as the gate: the tails beyond the quantiles must not lie beyond the end samples
    for cs_ in ('A', 'B'):
        r_ = RES[(cs_, nb_)]
        for j in range(nax):
            t = TZ[cs_][:, j]; lo, hi = np.percentile(t, [a.qlo, 100 - a.qlo]); U_ = r_['U'][j]
            def pos_(m): return np.mean(np.r_[U_[0][m[HS[0]['ti']]], U_[1][m[HS[1]['ti']]]])
            b0 = pos_(r_['idx'][j] == 0); short = max(short, b0 - pos_(t < lo), pos_(t > hi) - b0 - 1.)
    if short <= a.tail_tol: break
    log('a grid on %d axes leaves its end samples %.2f (end-to-end units) short of the tails (check 3, tolerance %.2f): dropping axis %d' % (nax, short, a.tail_tol, AX[-1] + 1))
if best['A'] is None or best['B'] is None: log('a fold resolves no candidate: no landscape on these axes'); raise SystemExit
nbest = min(best['A'], best['B']); nper = [nbest] * nax
# ---------------- final samples on the full coordinate
Tz = TZ['full']; cell0, idx = make_cells(Tz, nper, 'full'); cell, fs, merged = merge_sparse(Tz, cell0, 'all'); gs = list(np.unique(cell))
order = sorted(gs, key=lambda g: tuple(np.round(Tz[cell == g].mean(0), 3))); sid = {g: i + 1 for i, g in enumerate(order)}
out = np.zeros(NM, int); out[rows] = [sid[g] for g in cell]; np.savetxt(a.out + '/samples.txt', out, fmt='%d')
log('CHOSEN %s per axis -> %d samples (%d merged as not reproducible); labels %s/samples.txt' % ('x'.join(map(str, nper)), len(gs), len(merged), a.out))
log('check 1, reproducibility (full pixels):')
for g in order:
    s = [np.where(l == g)[0] for l in hl(cell)]; f = fsc(rec(HS[0], s[0]), rec(HS[1], s[1])); r143 = next((n * a.smpd / k for k in range(2, khi + 1) if f[k] < 0.143), n * a.smpd / khi)
    log('  sample %2d: %6d particles, centre %s sd, FSC at %.0f A %.2f, FSC 0.143 at %.1f A' % (sid[g], (cell == g).sum(), '/'.join('%+.2f' % v for v in Tz[cell == g].mean(0)), a.fsc_res, fs.get(g, np.nan), r143))
for coordset in ('A', 'B'):
    r_ = RES[(coordset, nbest)]; st = r_['st']
    log('check 2 + 4, fold %s at n %d: neighbour steps (end-to-end units) and z (threshold %.2f): %s' % (coordset, nbest, r_['zc'], ' '.join('%d-%d/ax%d:%.3f(z %.1f)' % (p, q, j + 1, *st[(p, q, j)][:2]) for p, q, j in r_['prs'])))
    for j in range(nax):
        d = np.array([st[k][0] for k in r_['prs'] if k[2] == j]); log('   axis %d: steps CV %.2f (mean %.3f, min %.3f, max %.3f)' % (AX[j] + 1, d.std() / max(abs(d.mean()), 1e-12), d.mean(), d.min(), d.max()))
    Tf = TZ[coordset]; test = 'B' if coordset == 'A' else 'A'
    for j in range(nax):   # check 3: tails vs end cells on this fold's test pixels
        t = Tf[:, j]; lo, hi = np.percentile(t, [a.qlo, 100 - a.qlo]); idxf = r_['idx']; U = r_['U'][j]
        def posof(m):
            v_ = [U[h][m[H['ti']]] for h, H in enumerate(HS)]; return np.mean(np.r_[v_[0], v_[1]]), np.sqrt(np.var(np.r_[v_[0], v_[1]]) / max(m.sum(), 1))
        base = posof(idxf[j] == 0)[0]
        sets = [('end_lo', idxf[j] == 0), ('end_hi', idxf[j] == nper[j] - 1), ('tail_lo(<%.1f%%)' % a.qlo, t < lo), ('tail_hi(>%.1f%%)' % (100 - a.qlo), t > hi), ('tail_lo0.5', t < np.percentile(t, a.qlo / 2)), ('tail_hi0.5', t > np.percentile(t, 100 - a.qlo / 2))]
        log('check 3, fold %s axis %d (positions on the test pixels, end cells at 0 / 1): %s' % (coordset, AX[j] + 1, ' '.join('%s %d %.2f+-%.2f' % (nm, m.sum(), posof(m)[0] - base, posof(m)[1]) for nm, m in sets)))
if GT is not None:
    log('TRUTH per final sample:')
    for g in order:
        m = cell == g; s = '  sample %2d (%6d):' % (sid[g], m.sum())
        for nm, x in GT['codes'].items(): xs = x[m]; f = np.isfinite(xs); s += ' %s %.1f+-%.1f (%.0f %% coded)' % (nm, np.nanmean(xs) if f.any() else np.nan, np.nanstd(xs) if f.any() else np.nan, 100 * f.mean())
        if GT['closed'].any(): s += ' closed %.3f (%d)' % (GT['closed'][m].mean(), int(GT['closed'][m].sum()))
        cc = np.bincount(GT['conf'][m].clip(0)); s += ' top %s' % ' '.join('%d:%.2f' % (k, cc[k] / m.sum()) for k in np.argsort(-cc)[:3] if cc[k] > 0); log(s)
    for nm, x in GT['codes'].items():
        f = np.isfinite(x); mus = [x[(cell == g) & f].mean() for g in order if ((cell == g) & f).sum() >= 100 and ((cell == g) & f).sum() >= 0.5 * (cell == g).sum()]
        W = sum(((x[(cell == g) & f] - x[(cell == g) & f].mean()) ** 2).sum() for g in gs if ((cell == g) & f).any()); e2 = 1 - W / ((x[f] - x[f].mean()) ** 2).sum()
        if len(mus) >= 2: st_ = np.diff(mus); log('TRUTH %s: samples >= 50 %% coded span %.1f .. %.1f, eta2 %.3f%s' % (nm, min(mus), max(mus), e2, '' if nax > 1 else ', steps %s (CV %.2f)' % (' '.join('%.1f' % v for v in st_), np.std(st_) / max(abs(np.mean(st_)), 1e-12))))
    if GT['closed'].any():
        k = max(order, key=lambda g: GT['closed'][cell == g].mean()); log('TRUTH closed: best sample %d: %d particles, closed share %.3f, %.1f %% of the node\'s closed' % (sid[k], (cell == k).sum(), GT['closed'][cell == k].mean(), 100 * GT['closed'][cell == k].sum() / max(GT['closed'].sum(), 1)))
    log('TRUTH NMI(conformer; sample) %.3f' % nmi(GT['conf'], cell))
np.savez(os.path.join(a.out, 'coords.npz'), rows=rows, tfull=TZ['full'], tA=TZ['A'], tB=TZ['B'], ax=np.array(AX), lam=lam, nper=np.array(nper))   # per-particle coordinates for cross-run comparisons
jobs = [(cs_, te_, 'heldout_%s' % cs_) for cs_, te_ in (('A', 'B'), ('B', 'A'))] + [(cs_, cs_, 'own_%s' % cs_) for cs_ in ('A', 'B')] + [('full', 'all', 'final_all')]
for coordset, test, sub in jobs:   # heldout_X: fold X's cells from the test pixels (no imprint); own_X: the same cells from the coordinate's own pixels; final_all: the final samples, all pixels
    if coordset == 'full': cellf = cell; gsf = order
    else: r_ = RES[(coordset, nbest)]; cellf = r_['cell']; gsf = sorted(np.unique(cellf), key=lambda g: tuple(np.round(TZ[coordset][cellf == g].mean(0), 3)))
    labh = hl(cellf); D = os.path.join(a.out, sub); os.makedirs(D, exist_ok=True); lab_out = np.zeros(NM, int); sid_ = {g: k for k, g in enumerate(gsf, 1)}; lab_out[rows] = [sid_[g] for g in cellf]
    np.savetxt(os.path.join(D, 'labels.txt'), lab_out, fmt='%d')
    for k, g in enumerate(gsf, 1):
        M_ = 0.5 * (rec(HS[0], np.where(labh[0] == g)[0], test) + rec(HS[1], np.where(labh[1] == g)[0], test)); write_mrc(os.path.join(D, 'state_%02d.mrc' % k), padmap(M_), n * a.smpd / NB)
    log('%s maps (cells of coordinate %s, pixels %s): %d samples in %s (sizes %s)' % (sub, coordset, test, len(gsf), D, ' '.join(str((cellf == g).sum()) for g in gsf)))
log('DONE')
