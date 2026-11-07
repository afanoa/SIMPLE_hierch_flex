#!/usr/bin/env python3
"""Signedness of a coordinate node's end-to-end difference (problem 1: continuous motion vs discrete composition, truth-free).
A motion moves density (positive and negative lobes of equal mass), an occupancy change adds or removes it (one sign). Per node: the
coordinate leaves' unfiltered half-maps in coordinate order (sphere-masked, low-passed at LP); end difference D = particle-weighted
mean of the last NEND leaves - mean of the first NEND, per half; significant voxels = cross-half product D_E D_O above the (1 - ALPHA)
quantile of its negative side (the noise); signedness s = (M+ - M-) / (M+ + M-) of the significant mass of the half-mean D (0 =
transport, +-1 = occupancy), the same with cross-half energy instead of mass, and the global balance sum(D) / (M+ + M-) over the
sphere. With truth: the composition-weighted truth maps of the same leaves through the same pipeline (support = the voxels holding
90 % of D_true^2) and each end's three largest conformers. Variants: truth with junk as zero density (the dilution the data see) and
'reg' = every end mean regressed per Fourier shell on the node mean first (cross-half coefficient, node_masks2) so global scale /
blur / dilution differences between the ends drop out (data and junk-as-zero truth alike).
usage: p1_sign.py OUT --node NAME:RECDIR:LABELS [--node ...] [--gt GT --truth DIR] [--mskdiam 220 --lp 20,15,10 --nend 2 --alpha 0.001]"""
import os, glob, argparse, warnings, numpy as np
warnings.filterwarnings('ignore')
ap = argparse.ArgumentParser(); ap.add_argument('out'); ap.add_argument('--node', action='append', required=True)
ap.add_argument('--gt', default=''); ap.add_argument('--truth', default=''); ap.add_argument('--mskdiam', type=float, default=220.)
ap.add_argument('--lp', default='20,15,10'); ap.add_argument('--nend', type=int, default=2); ap.add_argument('--alpha', type=float, default=0.001)
a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
def rd(p):
    h = np.fromfile(p, np.int32, 256); n = int(h[0]); return np.fromfile(p, np.float32, offset=1024 + int(h[23]), count=n ** 3).reshape(n, n, n).astype(np.float64)
G = np.loadtxt(a.gt, comments='#') if a.gt else None; tcache = {}
def tmap(c, n):
    if c not in tcache:
        fs = sorted(glob.glob('%s/conf%02d_*m160.mrc' % (a.truth, c))) if a.truth else []; tcache[c] = rd(fs[0]) if fs else None
        if tcache[c] is not None and tcache[c].shape[0] != n: tcache[c] = None
    return tcache[c]
out = ['# node lp variant | sig voxels %% of sphere | data: s_mass s_energy balance | truth (junk excluded): s_mass s_energy balance | truth (junk = 0): same | end difference / mean (cross-half rms ratio)']
for spec in a.node:
    name, rdir, labf = spec.split(':'); lab = np.loadtxt(labf, dtype=int); K = lab.max(); cnt = np.array([(lab == k).sum() for k in range(1, K + 1)], float)
    f0 = '%s/recvol_state01_even_unfil.mrc' % rdir; n = int(np.fromfile(f0, np.int32, 1)[0]); smpd = float(np.fromfile(f0, np.float32, 13)[10]) / n
    g = np.arange(n) - n // 2; R = np.sqrt(g[:, None, None] ** 2 + g[None, :, None] ** 2 + g[None, None, :] ** 2); r0 = a.mskdiam / 2 / smpd
    M = np.where(R <= r0, 1., np.where(R <= r0 + 6, 0.5 + 0.5 * np.cos(np.pi * (R - r0) / 6), 0.)); inm = M > 0.5
    f = np.fft.fftfreq(n); kx, ky, kz = np.meshgrid(f, f, np.fft.rfftfreq(n), indexing='ij'); kr = np.sqrt(kx ** 2 + ky ** 2 + kz ** 2)
    EF = [np.fft.rfftn(M * rd('%s/recvol_state%02d_even_unfil.mrc' % (rdir, k))) for k in range(1, K + 1)]
    OF = [np.fft.rfftn(M * rd('%s/recvol_state%02d_odd_unfil.mrc' % (rdir, k))) for k in range(1, K + 1)]
    lo, hi = list(range(a.nend)), list(range(K - a.nend, K)); wl, wh = cnt[lo] / cnt[lo].sum(), cnt[hi] / cnt[hi].sum(); w = cnt / cnt.sum()
    TF = None; TZ = None; ends = ''
    if G is not None:
        conf = np.full(len(lab), -1); gi = G[:, 0].astype(int) - 1; ok = (gi >= 0) & (gi < len(lab)); conf[gi[ok]] = G[ok, 1].astype(int); TF = []; TZ = []
        for k in range(1, K + 1):
            cs, cn = np.unique(conf[lab == k], return_counts=True); good = [(c, m) for c, m in zip(cs, cn) if c > 0 and tmap(c, n) is not None]; ng = sum(m for _, m in good)
            TF.append(np.fft.rfftn(M * sum(m * tmap(c, n) for c, m in good) / ng) if good else np.zeros_like(EF[0])); TZ.append(TF[-1] * ng / max(cn.sum(), 1))
        for tag, idx in (('low end', lo), ('high end', hi)):
            m = np.isin(lab, [i + 1 for i in idx]); cs, cn = np.unique(conf[m], return_counts=True); o = np.argsort(-cn)[:3]
            ends += '   %s %s (%d ptcls): %s\n' % (name, tag, m.sum(), ' '.join('conf%d %.0f%%' % (cs[j], 100. * cn[j] / m.sum()) for j in o))
        tcache.clear()
    sh = np.minimum(np.round(kr * n).astype(int), n // 2); wg = np.where((kz == 0) | (np.abs(kz) == 0.5), 1., 2.)
    def ssum(z): return np.bincount(sh.ravel(), (wg * z).ravel(), n // 2 + 1)
    def endpair(F, Mx, My, reg):   # high end - low end; reg: each end minus its per-shell regression on the node mean (cross-half coefficient)
        xh, xl = sum(wh[i] * F[j] for i, j in enumerate(hi)), sum(wl[i] * F[j] for i, j in enumerate(lo))
        if not reg: return xh - xl
        den = ssum(np.real(Mx * My.conj())); okd = np.abs(den) > 1e-6 * np.abs(den).max(); den = np.where(okd, den, 1.)
        b = lambda x: np.where(okd, ssum(np.real(x * My.conj())) / den, 0.)[sh]
        return (xh - b(xh) * Mx) - (xl - b(xl) * Mx)
    MEF, MOF = sum(w[j] * EF[j] for j in range(K)), sum(w[j] * OF[j] for j in range(K))
    for lpv, reg in [(float(x), r_) for x in a.lp.split(',') for r_ in (False, True)]:
        lpf = np.exp(-(kr / smpd * lpv) ** 2 * np.log(2)); back = lambda F: np.fft.irfftn(F * lpf, s=(n, n, n))
        DE = back(endpair(EF, MEF, MOF, reg)); DO = back(endpair(OF, MOF, MEF, reg))
        mE, mO = back(MEF), back(MOF)
        v = DE * DO; thr = np.quantile(-v[inm], 1 - a.alpha); S = (v > thr) & inm; Db = (DE + DO) / 2
        def signs(D, sup, en):
            mp, mm = D[sup & (D > 0)].sum(), -D[sup & (D < 0)].sum(); ep, em = en[sup & (D > 0)].sum(), en[sup & (D < 0)].sum()
            return (mp - mm) / max(mp + mm, 1e-300), (ep - em) / max(ep + em, 1e-300), D[inm].sum() / max(mp + mm, 1e-300)
        sd = signs(Db, S, v); rel = np.sqrt(max(v[inm].sum(), 0) / max((mE * mO)[inm].sum(), 1e-300))
        st, sz = (np.nan,) * 3, (np.nan,) * 3
        def tsign(TT):
            MT = sum(w[j] * TT[j] for j in range(K)); DT = back(endpair(TT, MT, MT, reg)); e2 = np.where(inm, DT ** 2, 0.)
            o = np.argsort(-e2.ravel()); cs = np.cumsum(e2.ravel()[o]) / max(e2.sum(), 1e-300); sup = np.zeros(e2.size, bool); sup[o[:np.searchsorted(cs, 0.9) + 1]] = True
            return signs(DT, sup.reshape(e2.shape), e2)
        if TF is not None: st, sz = (tsign(TF) if not reg else (np.nan,) * 3), tsign(TZ)
        out.append('%s %4.0f %-3s | %5.2f | %+.3f %+.3f %+.3f | %+.3f %+.3f %+.3f | %+.3f %+.3f %+.3f | %.3f' % (name, lpv, 'reg' if reg else 'raw', 100. * S.sum() / inm.sum(), sd[0], sd[1], sd[2], st[0], st[1], st[2], sz[0], sz[1], sz[2], rel))
    if ends: out.append(ends.rstrip('\n'))
open(a.out + '/signedness.txt', 'w').write('\n'.join(out) + '\n'); print('\n'.join(out))
