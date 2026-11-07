#!/usr/bin/env python3
"""Geometry of a node's K hard states (its focused refine3D: recvol_stateNN_{even,odd}_unfil.mrc, any box; truth maps are Fourier-
cropped to it): do the states lie along ONE path (a continuum sampled by states) or spread in several directions (a composition)?
Per band, after per-shell regression of every state on the members-weighted mean: the cross-half Gram matrix of the deviations
G_kl = Re<dE_k, dO_l> (symmetrised; noise-free in expectation because the halves are independent) -> its positive eigenvalue
shares: s1 = 1-D share (one path / a line between two groups), s2 = 2-D share, and R = cross-half reproducibility. The refinement
assigned with lp = LPA (10 A), so only bands finer than LPA are free of selection bias. With --gt the same on the members'
conformer-weighted truth maps (noise-free targets) and the truth composition.
usage: nodetype_states.py RECDIR LABELS [--lpa 10] [--mskdiam 220] [--gt gt --tmaps DIR] [--tag NAME]"""
import argparse, glob, os, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('rec'); ap.add_argument('labels'); ap.add_argument('--lpa', type=float, default=10.); ap.add_argument('--mskdiam', type=float, default=220.)
ap.add_argument('--gt', default=''); ap.add_argument('--tmaps', default=''); ap.add_argument('--tag', default=''); a = ap.parse_args()
def rd(p):
    h = np.fromfile(p, np.int32, 256); n = int(h[0]); return np.fromfile(p, np.float32, offset=1024 + int(h[23]), count=n ** 3).reshape(n, n, n).astype(np.float64), float(np.fromfile(p, np.float32, 13)[10]) / n
lab = np.loadtxt(a.labels).astype(int); ks = [k for k in range(1, lab.max() + 1) if (lab == k).any() and os.path.exists(os.path.join(a.rec, 'recvol_state%02d_even_unfil.mrc' % k))]
w = np.array([(lab == k).sum() for k in ks], float); w /= w.sum()
v0, smpd = rd(os.path.join(a.rec, 'recvol_state%02d_even_unfil.mrc' % ks[0])); n = v0.shape[0]; c = n // 2
g = np.arange(n) - c; R3 = np.sqrt(g[:, None, None] ** 2 + g[None, :, None] ** 2 + g[None, None, :] ** 2); r0 = a.mskdiam / 2 / smpd
M = np.where(R3 <= r0, 1., np.where(R3 <= r0 + 6, 0.5 + 0.5 * np.cos(np.pi * (R3 - r0) / 6), 0.))
f = np.fft.fftfreq(n); kx, ky, kz = np.meshgrid(f, f, np.fft.rfftfreq(n), indexing='ij'); kr = np.sqrt(kx ** 2 + ky ** 2 + kz ** 2)
res = np.where(kr > 0, smpd / np.maximum(kr, 1e-9), 1e9); keep = res >= max(2.05 * smpd, 5.9); sh = np.round(kr * n).astype(int)[keep]; nsh = sh.max() + 1
wb = np.where((kz == 0) | (np.abs(kz) == 0.5), 1., 2.)[keep]; resk = res[keep]
bands = {'assign <=%g A (biased)' % a.lpa: (117., a.lpa), 'free %g-6 A' % a.lpa: (a.lpa, 6.), 'free %g-8 A' % a.lpa: (a.lpa, 8.)}
def F(x): return np.fft.rfftn(M * x)[keep]
def crop(v):
    m = v.shape[0]
    if m == n: return v
    Fv = np.fft.fftshift(np.fft.fftn(np.fft.ifftshift(v))); o = m // 2 - n // 2; return np.real(np.fft.fftshift(np.fft.ifftn(np.fft.ifftshift(Fv[o:o + n, o:o + n, o:o + n]))))
def regress(E, O):
    Em, Om = w @ E, w @ O; ss = lambda z: np.bincount(sh, weights=wb * z, minlength=nsh)
    den = ss(np.real(Em * Om.conj())); ok = np.abs(den) > 1e-9 * np.abs(den).max(); den = np.where(ok, den, 1.)
    return (np.array([e - np.where(ok, ss(np.real(e * Om.conj())) / den, 1.)[sh] * Em for e in E]), np.array([o - np.where(ok, ss(np.real(o * Em.conj())) / den, 1.)[sh] * Om for o in O]))
def geo(dE, dO, lo, hi):
    m = (resk <= lo) & (resk >= hi); ww = wb[m]; A = dE[:, m] * np.sqrt(w)[:, None]; B = dO[:, m] * np.sqrt(w)[:, None]   # member-weighted
    G = np.real((A * ww) @ B.conj().T); G = 0.5 * (G + G.T); lam = np.sort(np.linalg.eigvalsh(G))[::-1]; pos = np.maximum(lam, 0); tot = max(pos.sum(), 1e-300)
    Rv = np.trace(G) / np.sqrt(np.real((np.abs(A) ** 2) @ ww).sum() * np.real((np.abs(B) ** 2) @ ww).sum() + 1e-300)
    return Rv, pos[0] / tot, pos[:2].sum() / tot, lam[-1] / max(lam[0], 1e-300)
E = np.array([F(rd(os.path.join(a.rec, 'recvol_state%02d_even_unfil.mrc' % k))[0]) for k in ks]); O = np.array([F(rd(os.path.join(a.rec, 'recvol_state%02d_odd_unfil.mrc' % k))[0]) for k in ks])
dE, dO = regress(E, O); tag = (a.tag + ' ') if a.tag else ''
print('%s%d states (members %s), box %d @ %.2f A; assignment lp %g A' % (tag, len(ks), ' '.join(str((lab == k).sum()) for k in ks), n, smpd, a.lpa))
for bn, (lo, hi) in bands.items():
    Rv, s1, s2, neg = geo(dE, dO, lo, hi); print('%s  %-24s R %.3f  1-D share %.3f  2-D share %.3f  (most negative eigenvalue / top %.3f)' % (tag, bn, Rv, s1, s2, neg))
if a.gt and a.tmaps:
    Gt = np.loadtxt(a.gt, comments='#', ndmin=2); conf = np.zeros(len(lab), int); conf[Gt[:, 0].astype(int) - 1] = Gt[:, 1].astype(int); cache = {}; X = []
    for k in ks:
        cs, cn = np.unique(conf[(lab == k) & (conf > 0)], return_counts=True); acc, ws = 0, 0
        for ci, ni in zip(cs, cn):
            if ci not in cache:
                fl = glob.glob(os.path.join(a.tmaps, ('conf%02d_*_m160.mrc' % ci) if ci < 100 else ('conf%d_*_m160.mrc' % ci))); cache[ci] = F(crop(rd(fl[0])[0])) if fl else None
            if cache[ci] is not None: acc = acc + ni * cache[ci]; ws += ni
        X.append(acc / max(ws, 1))
    dX, _ = regress(np.array(X), np.array(X))
    for bn, (lo, hi) in bands.items():
        Rv, s1, s2, neg = geo(dX, dX, lo, hi); print('%s  TRUTH-maps %-13s          1-D share %.3f  2-D share %.3f' % (tag, bn, s1, s2))
    sw = (conf >= 200) & (conf <= 230)
    print('%s  TRUTH per state: swing share / mean swing / top conformer share: %s' % (tag, ' '.join('%.0f%%/%s/%.2f' % (100 * sw[lab == k].mean(), ('%.1f' % (conf[(lab == k) & sw] - 200).mean()) if (sw & (lab == k)).any() else '-',
          np.bincount(conf[(lab == k) & (conf > 0)]).max() / max(((lab == k) & (conf > 0)).sum(), 1)) for k in ks)))
