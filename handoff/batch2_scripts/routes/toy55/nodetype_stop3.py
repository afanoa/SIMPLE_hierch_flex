#!/usr/bin/env python3
"""Stop test v3 (traps from the split-decision session): a K-group organisation (REC dir with recvol_stateNN_{even,odd}_unfil.mrc +
per-row LABELS). (1) the test band starts >= 2 shells past the assignment band edge (soft masks couple neighbouring shells): pass the
bands accordingly; (2) groups sorted by the particles' pose errors give translated copies that reproduce in both halves, so with
--shiftfit every group's half-map is first aligned by a sub-pixel 3-D translation fitted IN THE TEST BAND against the OTHER half's
consensus (members-weighted mean of the other half's groups). Then per-shell regression on the mean, R = cross-half reproducibility of
the deviations, SE from 8 Fourier sectors, z = R / SE. Prints R, z without and with the shift fit and the fitted shifts.
usage: nodetype_stop3.py REC LABELS [--bands 13:9.5] [--mskdiam 220] [--tag NAME]"""
import argparse, os, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('rec'); ap.add_argument('labels'); ap.add_argument('--bands', default='13:9.5'); ap.add_argument('--mskdiam', type=float, default=220.)
ap.add_argument('--tag', default=''); a = ap.parse_args()
def rd(p):
    h = np.fromfile(p, np.int32, 256); n = int(h[0]); return np.fromfile(p, np.float32, offset=1024 + int(h[23]), count=n ** 3).reshape(n, n, n).astype(np.float64), float(np.fromfile(p, np.float32, 13)[10]) / n
lab = np.loadtxt(a.labels).astype(int); ks = [k for k in range(1, lab.max() + 1) if (lab == k).any() and os.path.exists(os.path.join(a.rec, 'recvol_state%02d_even_unfil.mrc' % k))]
w = np.array([(lab == k).sum() for k in ks], float); w /= w.sum()
v0, smpd = rd(os.path.join(a.rec, 'recvol_state%02d_even_unfil.mrc' % ks[0])); n = v0.shape[0]; c = n // 2
g = np.arange(n) - c; R3 = np.sqrt(g[:, None, None] ** 2 + g[None, :, None] ** 2 + g[None, None, :] ** 2); r0 = a.mskdiam / 2 / smpd
M = np.where(R3 <= r0, 1., np.where(R3 <= r0 + 6, 0.5 + 0.5 * np.cos(np.pi * (R3 - r0) / 6), 0.))
f = np.fft.fftfreq(n); fz = np.fft.rfftfreq(n); KX, KY, KZ = np.meshgrid(f, f, fz, indexing='ij'); kr = np.sqrt(KX ** 2 + KY ** 2 + KZ ** 2); res = np.where(kr > 0, smpd / np.maximum(kr, 1e-9), 1e9)
wb_full = np.where((KZ == 0) | (np.abs(KZ) == 0.5), 1., 2.); sh_full = np.round(kr * n).astype(int); tag = (a.tag + ' ') if a.tag else ''
E = np.array([np.fft.rfftn(M * rd(os.path.join(a.rec, 'recvol_state%02d_even_unfil.mrc' % k))[0]) for k in ks]); O = np.array([np.fft.rfftn(M * rd(os.path.join(a.rec, 'recvol_state%02d_odd_unfil.mrc' % k))[0]) for k in ks])
def fitshift(X, Ref, band):   # sub-pixel 3-D translation maximising the band-limited cross-correlation of X with Ref (within +-4 px)
    P = np.zeros_like(X); P[band] = X[band] * Ref[band].conj(); cc = np.fft.irfftn(P, s=(n, n, n)); cc = np.roll(cc, (4, 4, 4), axis=(0, 1, 2))[:9, :9, :9]
    i = np.array(np.unravel_index(np.argmax(cc), cc.shape)); t = i - 4.
    for d in range(3):
        if 0 < i[d] < 8:
            lo_ = cc[tuple(np.where(np.arange(3) == d, i - 1, i))]; mi = cc[tuple(i)]; hi_ = cc[tuple(np.where(np.arange(3) == d, i + 1, i))]; den = lo_ - 2 * mi + hi_
            if den < 0: t[d] += 0.5 * (lo_ - hi_) / den
    return t
def shifted(X, t): return X * np.exp(2j * np.pi * (KX * t[0] + KY * t[1] + KZ * t[2]))   # the cross-correlation peak t means x(r) ~ ref(r - t): undo it
def stat(Ex, Ox, band):
    Em, Om = w @ Ex.reshape(len(ks), -1), w @ Ox.reshape(len(ks), -1); Em = Em.reshape(Ex.shape[1:]); Om = Om.reshape(Ox.shape[1:])
    m = band; sh = sh_full[m]; wb = wb_full[m]; nsh = sh_full.max() + 1; ss = lambda z: np.bincount(sh, weights=wb * z, minlength=nsh)
    den = ss(np.real(Em[m] * Om[m].conj())); ok = np.abs(den) > 1e-9 * np.abs(den).max(); den = np.where(ok, den, 1.)
    dE = np.array([e[m] - np.where(ok, ss(np.real(e[m] * Om[m].conj())) / den, 1.)[sh] * Em[m] for e in Ex]); dO = np.array([o[m] - np.where(ok, ss(np.real(o[m] * Em[m].conj())) / den, 1.)[sh] * Om[m] for o in Ox])
    x = (w[:, None] * np.real(dE * dO.conj())).sum(0) * wb; pe = (w[:, None] * np.abs(dE) ** 2).sum(0) * wb; po = (w[:, None] * np.abs(dO) ** 2).sum(0) * wb
    sec = ((KX[m] >= 0).astype(int) + 2 * (KY[m] >= 0) + 4 * (np.abs(KZ[m]) >= np.abs(KZ[m]).max() / 2)).astype(int)
    Rv = x.sum() / np.sqrt(pe.sum() * po.sum()); Rs = np.array([x[sec == s].sum() / np.sqrt(pe[sec == s].sum() * po[sec == s].sum()) for s in range(8)]); se = Rs.std(ddof=1) / np.sqrt(8)
    return Rv, se
for b in a.bands.split(','):
    lo, hi = [float(t) for t in b.split(':')]; band = (res <= lo) & (res >= hi)
    R0, se0 = stat(E, O, band)
    Em = (w[:, None, None, None] * E).sum(0); Om = (w[:, None, None, None] * O).sum(0)
    tE = [fitshift(E[i], Om, band) for i in range(len(ks))]; tO = [fitshift(O[i], Em, band) for i in range(len(ks))]
    Es = np.array([shifted(E[i], tE[i]) for i in range(len(ks))]); Os = np.array([shifted(O[i], tO[i]) for i in range(len(ks))])
    R1, se1 = stat(Es, Os, band)
    sz = np.array([np.linalg.norm(t) for t in tE + tO]); agree = np.median([np.linalg.norm(np.array(tE[i]) - np.array(tO[i])) for i in range(len(ks))])
    print('%sband %g-%g A: R %.4f z %.1f | after per-group 3-D shift fit R %.4f z %.1f | fitted |shift| median %.2f px (max %.2f), even-odd shift disagreement median %.2f px'
          % (tag, lo, hi, R0, R0 / se0, R1, R1 / se1, np.median(sz), sz.max(), agree))
