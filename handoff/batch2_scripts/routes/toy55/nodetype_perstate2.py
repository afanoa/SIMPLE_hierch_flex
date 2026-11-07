#!/usr/bin/env python3
"""Per-state stop gate (nodetype_perstate.py) plus two readouts. --tmaps: with --gt, each state's TRUTH deviation share in the band
(members' conformer-weighted truth maps, Fourier-cropped to the data box, per-shell regression on the leave-one-out truth mean, relative
deviation energy |dX_k|^2 / |X_k|^2; junk members carry no truth map) = how much the state really differs from the others there.
--cmp LABELS: where each state's members came from in another labelling of the same rows (top two sources) and the NMI between them.
Gate: every hard state's halves vs the other half's leave-one-out mean (Newton-refined 3-D shift fit first), per-shell regression,
R_k, SE_k from NSEC equal-solid-angle Fourier sectors, z_k, one-sided p_k from Student t (NSEC - 1 dof); GO if K * min_k p_k < ALPHA.
usage: nodetype_perstate2.py REC LABELS [--bands 17:10,13:10] [--mskdiam 220] [--nsec 32] [--alpha 0.00135] [--gt gt [--tmaps DIR]] [--cmp LABELS] [--tag NAME]"""
import argparse, glob, os, numpy as np
from scipy import stats
ap = argparse.ArgumentParser(); ap.add_argument('rec'); ap.add_argument('labels'); ap.add_argument('--bands', default='17:10,13:10'); ap.add_argument('--mskdiam', type=float, default=220.)
ap.add_argument('--nsec', type=int, default=32); ap.add_argument('--alpha', type=float, default=0.00135); ap.add_argument('--gt', default=''); ap.add_argument('--tmaps', default='')
ap.add_argument('--cmp', default=''); ap.add_argument('--tag', default=''); a = ap.parse_args()
def rd(p):
    h = np.fromfile(p, np.int32, 256); n = int(h[0]); return np.fromfile(p, np.float32, offset=1024 + int(h[23]), count=n ** 3).reshape(n, n, n).astype(np.float64), float(np.fromfile(p, np.float32, 13)[10]) / n
def nmi(x, y):
    _, xi = np.unique(x, return_inverse=True); _, yi = np.unique(y, return_inverse=True); J = np.zeros((xi.max() + 1, yi.max() + 1)); np.add.at(J, (xi, yi), 1); J /= J.sum()
    px, py = J.sum(1), J.sum(0); nz = J > 0; I = (J[nz] * np.log(J[nz] / np.outer(px, py)[nz])).sum(); h = lambda p: -(p[p > 0] * np.log(p[p > 0])).sum(); return I / max(np.sqrt(h(px) * h(py)), 1e-12)
lab = np.loadtxt(a.labels).astype(int); ks = [k for k in range(1, lab.max() + 1) if (lab == k).any() and os.path.exists(os.path.join(a.rec, 'recvol_state%02d_even_unfil.mrc' % k))]
nk = np.array([(lab == k).sum() for k in ks]); w = nk / nk.sum(); K = len(ks)
v0, smpd = rd(os.path.join(a.rec, 'recvol_state%02d_even_unfil.mrc' % ks[0])); n = v0.shape[0]; c = n // 2
g = np.arange(n) - c; R3 = np.sqrt(g[:, None, None] ** 2 + g[None, :, None] ** 2 + g[None, None, :] ** 2); r0 = a.mskdiam / 2 / smpd
M = np.where(R3 <= r0, 1., np.where(R3 <= r0 + 6, 0.5 + 0.5 * np.cos(np.pi * (R3 - r0) / 6), 0.))
f = np.fft.fftfreq(n); fz = np.fft.rfftfreq(n); KX, KY, KZ = np.meshgrid(f, f, fz, indexing='ij'); kr = np.sqrt(KX ** 2 + KY ** 2 + KZ ** 2); res = np.where(kr > 0, smpd / np.maximum(kr, 1e-9), 1e9)
wb_full = np.where((KZ == 0) | (np.abs(KZ) == 0.5), 1., 2.); sh_full = np.round(kr * n).astype(int); nsh = sh_full.max() + 1; tag = (a.tag + ' ') if a.tag else ''
npol = max(a.nsec // 8, 1); nsec = 8 * npol   # cos(polar angle) of the half-space in equal bins = equal solid angle
sec_full = (np.floor((np.arctan2(KY, KX) + np.pi) / (2 * np.pi) * 8).astype(int) % 8) + 8 * np.minimum(np.floor(KZ / np.maximum(kr, 1e-12) * npol).astype(int), npol - 1)
E = np.array([np.fft.rfftn(M * rd(os.path.join(a.rec, 'recvol_state%02d_even_unfil.mrc' % k))[0]) for k in ks]); O = np.array([np.fft.rfftn(M * rd(os.path.join(a.rec, 'recvol_state%02d_odd_unfil.mrc' % k))[0]) for k in ks])
def fitshift(X, Ref, band, newton=True):   # sub-pixel 3-D translation maximising the band-limited cross-correlation of X with Ref (within +-4 px)
    P = np.zeros_like(X); P[band] = X[band] * Ref[band].conj(); cc = np.fft.irfftn(P, s=(n, n, n)); cc = np.roll(cc, (4, 4, 4), axis=(0, 1, 2))[:9, :9, :9]
    i = np.array(np.unravel_index(np.argmax(cc), cc.shape)); t = i - 4.
    for d in range(3):
        if 0 < i[d] < 8:
            lo_ = cc[tuple(np.where(np.arange(3) == d, i - 1, i))]; mi = cc[tuple(i)]; hi_ = cc[tuple(np.where(np.arange(3) == d, i + 1, i))]; den = lo_ - 2 * mi + hi_
            if den < 0: t[d] += 0.5 * (lo_ - hi_) / den
    if newton:   # the parabolic peak of a band-pass correlation is biased the same way in both halves: Newton on the exact correlation
        q0 = wb_full[band] * P[band]; kk = np.stack([KX[band], KY[band], KZ[band]])
        for _ in range(8):
            q = q0 * np.exp(2j * np.pi * (kk.T @ t)); gr = -2 * np.pi * (kk * q.imag).sum(1); H = -4 * np.pi ** 2 * (kk * q.real) @ kk.T
            if np.linalg.eigvalsh(H).max() >= 0: break
            dt = -np.linalg.solve(H, gr); t = t + dt
            if np.linalg.norm(dt) < 1e-4: break
    return t
def shifted(X, t): return X * np.exp(2j * np.pi * (KX * t[0] + KY * t[1] + KZ * t[2]))
def stat_global(Ex, Ox, band):   # nodetype_stop3's all-state statistic (deviations from the weighted mean of all states, 8 sectors)
    Em, Om = w @ Ex.reshape(K, -1), w @ Ox.reshape(K, -1); Em = Em.reshape(Ex.shape[1:]); Om = Om.reshape(Ox.shape[1:])
    m = band; sh = sh_full[m]; wb = wb_full[m]; ss = lambda z: np.bincount(sh, weights=wb * z, minlength=nsh)
    den = ss(np.real(Em[m] * Om[m].conj())); ok = np.abs(den) > 1e-9 * np.abs(den).max(); den = np.where(ok, den, 1.)
    dE = np.array([e[m] - np.where(ok, ss(np.real(e[m] * Om[m].conj())) / den, 1.)[sh] * Em[m] for e in Ex]); dO = np.array([o[m] - np.where(ok, ss(np.real(o[m] * Em[m].conj())) / den, 1.)[sh] * Om[m] for o in Ox])
    x = (w[:, None] * np.real(dE * dO.conj())).sum(0) * wb; pe = (w[:, None] * np.abs(dE) ** 2).sum(0) * wb; po = (w[:, None] * np.abs(dO) ** 2).sum(0) * wb
    sec = ((KX[m] >= 0).astype(int) + 2 * (KY[m] >= 0) + 4 * (np.abs(KZ[m]) >= np.abs(KZ[m]).max() / 2)).astype(int)
    Rv = x.sum() / np.sqrt(pe.sum() * po.sum()); Rs = np.array([x[sec == s].sum() / np.sqrt(pe[sec == s].sum() * po[sec == s].sum()) for s in range(8)]); return Rv, Rv / (Rs.std(ddof=1) / np.sqrt(8))
def loo(i, X, ww=None):
    wl = (w if ww is None else ww).copy(); wl[i] = 0.; wl /= max(wl.sum(), 1e-300); return np.tensordot(wl, X, 1)
def stat_state(i, Es, Os, band):
    m = band; e = Es[i][m]; o = Os[i][m]; eb = loo(i, Es)[m]; ob = loo(i, Os)[m]
    sh = sh_full[m]; wb = wb_full[m]; sec = sec_full[m]; ss = lambda z: np.bincount(sh, weights=wb * z, minlength=nsh)
    den = ss(np.real(eb * ob.conj())); ok = np.abs(den) > 1e-9 * np.abs(den).max(); den = np.where(ok, den, 1.)
    dE = e - np.where(ok, ss(np.real(e * ob.conj())) / den, 1.)[sh] * eb; dO = o - np.where(ok, ss(np.real(o * eb.conj())) / den, 1.)[sh] * ob
    x = np.real(dE * dO.conj()) * wb; pe = np.abs(dE) ** 2 * wb; po = np.abs(dO) ** 2 * wb
    Rv = x.sum() / np.sqrt(pe.sum() * po.sum()); Rs = np.array([x[sec == s].sum() / np.sqrt(pe[sec == s].sum() * po[sec == s].sum()) for s in range(nsec)])
    z = Rv / (Rs.std(ddof=1) / np.sqrt(nsec)); return Rv, z, stats.t.sf(z, nsec - 1)
truth = [''] * K; src = [''] * K; X = None; extra = ''
if a.gt:
    Gt = np.loadtxt(a.gt, comments='#', ndmin=2); conf = np.zeros(max(len(lab), int(Gt[:, 0].max())), int); conf[Gt[:, 0].astype(int) - 1] = Gt[:, 1].astype(int); conf = conf[:len(lab)]
    for i, k in enumerate(ks):
        cm = conf[lab == k]; ok = cm[cm > 0]; sw = (ok >= 200) & (ok <= 230); u, cn = np.unique(ok, return_counts=True)
        top = '%d:%.0f%%' % (u[np.argmax(cn)], 100 * cn.max() / len(ok)) if len(ok) else '-'
        truth[i] = '| truth junk %4.1f%% top %-8s swing %3.0f%% (mean %s)' % (100 * (cm == 0).mean(), top, 100 * sw.mean() if len(ok) else 0, ('%.1f' % (ok[sw] - 200).mean()) if sw.any() else '-')
    if a.tmaps:
        def crop(v):
            m_ = v.shape[0]
            if m_ == n: return v
            Fv = np.fft.fftshift(np.fft.fftn(np.fft.ifftshift(v))); o_ = m_ // 2 - n // 2; return np.real(np.fft.fftshift(np.fft.ifftn(np.fft.ifftshift(Fv[o_:o_ + n, o_:o_ + n, o_:o_ + n]))))
        cache = {}; X = []; wt = []
        for k in ks:
            cs, cn = np.unique(conf[(lab == k) & (conf > 0)], return_counts=True); acc = np.zeros_like(E[0]); ws = 0
            for ci, ni in zip(cs, cn):
                if ci not in cache:
                    fl = glob.glob(os.path.join(a.tmaps, ('conf%02d_*_m160.mrc' % ci) if ci < 100 else ('conf%d_*_m160.mrc' % ci)))
                    if fl and not extra: tv, ts_ = rd(fl[0]); extra = ', truth maps box %d @ %.2f A' % (tv.shape[0], ts_)
                    cache[ci] = np.fft.rfftn(M * crop(rd(fl[0])[0])) if fl else None
                if cache[ci] is not None: acc = acc + ni * cache[ci]; ws += ni
            X.append(acc / max(ws, 1)); wt.append(float(ws))
        X = np.array(X); wt = np.array(wt)
def tdev(i, band):   # relative deviation energy of state i's truth map from the leave-one-out truth mean after per-shell regression
    if wt[i] == 0 or wt.sum() == wt[i]: return np.nan
    xb = loo(i, X, wt)[band]; x = X[i][band]; sh = sh_full[band]; wb = wb_full[band]; ss = lambda z: np.bincount(sh, weights=wb * z, minlength=nsh)
    den = ss(np.abs(xb) ** 2); ok = den > 1e-12 * den.max(); den = np.where(ok, den, 1.); d = x - np.where(ok, ss(np.real(x * xb.conj())) / den, 0.)[sh] * xb
    return (np.abs(d) ** 2 * wb).sum() / max((np.abs(x) ** 2 * wb).sum(), 1e-300)
if a.cmp:
    lab2 = np.loadtxt(a.cmp).astype(int)[:len(lab)]; mm = (lab > 0) & (lab2 > 0); extra += ', NMI vs %s %.3f' % (os.path.basename(os.path.dirname(os.path.abspath(a.cmp))) + '/' + os.path.basename(a.cmp), nmi(lab[mm], lab2[mm]))
    for i, k in enumerate(ks):
        u, cn = np.unique(lab2[(lab == k) & (lab2 > 0)], return_counts=True); o_ = np.argsort(-cn)[:2]; src[i] = '| from ' + ' '.join('%d:%.0f%%' % (u[j], 100 * cn[j] / cn.sum()) for j in o_)
print('%s%d states (members %s), box %d @ %.2f A, %d sectors, alpha %g%s' % (tag, K, ' '.join(str(x) for x in nk), n, smpd, nsec, a.alpha, extra))
for b in a.bands.split(','):
    lo, hi = [float(t) for t in b.split(':')]; band = (res <= lo) & (res >= hi)
    tE = [fitshift(E[i], (w[:, None, None, None] * O).sum(0), band, False) for i in range(K)]; tO = [fitshift(O[i], (w[:, None, None, None] * E).sum(0), band, False) for i in range(K)]
    Rg, zg = stat_global(np.array([shifted(E[i], tE[i]) for i in range(K)]), np.array([shifted(O[i], tO[i]) for i in range(K)]), band)
    Es, Os = E, O
    for _ in range(2):   # every state's halves against the other half's leave-one-out mean, then again against the shift-corrected means
        tE = [fitshift(E[i], loo(i, Os), band) for i in range(K)]; tO = [fitshift(O[i], loo(i, Es), band) for i in range(K)]
        Es = np.array([shifted(E[i], tE[i]) for i in range(K)]); Os = np.array([shifted(O[i], tO[i]) for i in range(K)])
    rows = [stat_state(i, Es, Os, band) for i in range(K)]; p = np.array([r[2] for r in rows]); zthr = stats.t.isf(a.alpha / K, nsec - 1); go = K * p.min() < a.alpha; j = int(np.argmin(p))
    for i, k in enumerate(ks):
        Rv, z, pk = rows[i]; ts = max(np.linalg.norm(tE[i]), np.linalg.norm(tO[i])); td = (' tdev %.4f' % tdev(i, band)) if X is not None else ''
        print('%s  band %g-%g A state %2d n %6d  R %7.4f  z %5.1f  p %.1e  shift %.2f px %s %s%s %s' % (tag, lo, hi, k, nk[i], Rv, z, pk, ts, '*' if K * pk < a.alpha else ' ', truth[i], td, src[i]))
    print('%sGATE band %g-%g A: max z %.1f (state %d, n %d), K*min p %.1e, Bonferroni z threshold %.2f -> %s ; %d states flagged %s | global all-state R %.4f z %.1f'
          % (tag, lo, hi, rows[j][1], ks[j], nk[j], K * p.min(), zthr, 'GO' if go else 'STOP', int((K * p < a.alpha).sum()), [ks[i] for i in range(K) if K * p[i] < a.alpha], Rg, zg))
