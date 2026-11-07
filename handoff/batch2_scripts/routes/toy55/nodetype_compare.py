#!/usr/bin/env python3
"""Hard states vs coordinate bins at one node: which organisation captures more REPRODUCIBLE heterogeneity? For each organisation
(a K-state reconstruction dir + per-row labels): every state's even/odd maps, per-shell regression on the members-weighted mean, then
H = sum_k w_k Re<dE_k, dO_k> / Re<Ebar, Obar> in a band neither assignment used (default 10-8 A: the hard states were assigned at
lp 10, the cheap children at <= 15 A) = member-weighted between-group variance relative to the mean's power, noise-free in expectation.
Also label agreement between the organisations (NMI, and the share of the bins' coordinate variance the hard states explain).
With --gt: the same H on the members' conformer-weighted truth maps, and NMI(conformer) / swing eta2 per organisation.
usage: nodetype_compare.py --a NAME:RECDIR:LABELS --b NAME:RECDIR:LABELS [--coords coords.txt] [--band 10,8] [--mskdiam 220] [--gt gt --tmaps DIR]"""
import argparse, glob, os, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('--a', required=True); ap.add_argument('--b', required=True); ap.add_argument('--coords', default='')
ap.add_argument('--band', default='10,8'); ap.add_argument('--mskdiam', type=float, default=220.); ap.add_argument('--gt', default=''); ap.add_argument('--tmaps', default=''); ap.add_argument('--tag', default='')
a = ap.parse_args(); lo, hi = [float(x) for x in a.band.split(',')]; tag = (a.tag + ' ') if a.tag else ''
def rd(p):
    h = np.fromfile(p, np.int32, 256); n = int(h[0]); return np.fromfile(p, np.float32, offset=1024 + int(h[23]), count=n ** 3).reshape(n, n, n).astype(np.float64), float(np.fromfile(p, np.float32, 13)[10]) / n
conf = None
if a.gt:
    Gt = np.loadtxt(a.gt, comments='#', ndmin=2); NRg = int(Gt[:, 0].max()); conf = np.zeros(NRg, int); conf[Gt[:, 0].astype(int) - 1] = Gt[:, 1].astype(int)
def nmi(x, y):
    _, xi = np.unique(x, return_inverse=True); _, yi = np.unique(y, return_inverse=True); J = np.zeros((xi.max() + 1, yi.max() + 1)); np.add.at(J, (xi, yi), 1); J /= J.sum()
    px, py = J.sum(1), J.sum(0); nz = J > 0; I = (J[nz] * np.log(J[nz] / np.outer(px, py)[nz])).sum(); h = lambda p: -(p[p > 0] * np.log(p[p > 0])).sum(); return I / max(np.sqrt(h(px) * h(py)), 1e-12)
def organisation(spec):
    nm, rec, labf = spec.split(':'); lab = np.loadtxt(labf).astype(int); ks = [k for k in range(1, lab.max() + 1) if (lab == k).any() and os.path.exists(os.path.join(rec, 'recvol_state%02d_even_unfil.mrc' % k))]
    w = np.array([(lab == k).sum() for k in ks], float); w /= w.sum(); v0, smpd = rd(os.path.join(rec, 'recvol_state%02d_even_unfil.mrc' % ks[0])); n = v0.shape[0]; c = n // 2
    g = np.arange(n) - c; R3 = np.sqrt(g[:, None, None] ** 2 + g[None, :, None] ** 2 + g[None, None, :] ** 2); r0 = a.mskdiam / 2 / smpd
    M = np.where(R3 <= r0, 1., np.where(R3 <= r0 + 6, 0.5 + 0.5 * np.cos(np.pi * (R3 - r0) / 6), 0.))
    f = np.fft.fftfreq(n); kx, ky, kz = np.meshgrid(f, f, np.fft.rfftfreq(n), indexing='ij'); kr = np.sqrt(kx ** 2 + ky ** 2 + kz ** 2); res = np.where(kr > 0, smpd / np.maximum(kr, 1e-9), 1e9)
    keep = res >= hi - 1e-6; sh = np.round(kr * n).astype(int)[keep]; nsh = sh.max() + 1; wb = np.where((kz == 0) | (np.abs(kz) == 0.5), 1., 2.)[keep]; band = (res[keep] <= lo) & (res[keep] >= hi)
    F = lambda x: np.fft.rfftn(M * x)[keep]
    def crop(v):
        m = v.shape[0]
        if m == n: return v
        Fv = np.fft.fftshift(np.fft.fftn(np.fft.ifftshift(v))); o = m // 2 - n // 2; return np.real(np.fft.fftshift(np.fft.ifftn(np.fft.ifftshift(Fv[o:o + n, o:o + n, o:o + n]))))
    def H(E, O):
        Em, Om = w @ E, w @ O; ss = lambda z: np.bincount(sh, weights=wb * z, minlength=nsh); den = ss(np.real(Em * Om.conj())); ok = np.abs(den) > 1e-9 * np.abs(den).max(); den = np.where(ok, den, 1.)
        dE = np.array([e - np.where(ok, ss(np.real(e * Om.conj())) / den, 1.)[sh] * Em for e in E]); dO = np.array([o - np.where(ok, ss(np.real(o * Em.conj())) / den, 1.)[sh] * Om for o in O])
        num = sum(wk * np.real((dE[i][band] * dO[i][band].conj() * wb[band]).sum()) for i, wk in enumerate(w)); return num / np.real((Em[band] * Om[band].conj() * wb[band]).sum())
    E = np.array([F(rd(os.path.join(rec, 'recvol_state%02d_even_unfil.mrc' % k))[0]) for k in ks]); O = np.array([F(rd(os.path.join(rec, 'recvol_state%02d_odd_unfil.mrc' % k))[0]) for k in ks])
    out = dict(name=nm, lab=lab, K=len(ks), H=H(E, O), box=n)
    if conf is not None and a.tmaps:
        cache = {}; X = []
        for k in ks:
            cs, cn = np.unique(conf[(lab[:len(conf)] == k) & (conf > 0)], return_counts=True); acc, ws = 0, 0
            for ci, ni in zip(cs, cn):
                if ci not in cache:
                    fl = glob.glob(os.path.join(a.tmaps, ('conf%02d_*_m160.mrc' % ci) if ci < 100 else ('conf%d_*_m160.mrc' % ci))); cache[ci] = F(crop(rd(fl[0])[0])) if fl else None
                if cache[ci] is not None: acc = acc + ni * cache[ci]; ws += ni
            X.append(acc / max(ws, 1))
        out['Htrue'] = H(np.array(X), np.array(X)); m = (lab[:len(conf)] > 0) & (conf > 0); out['nmi'] = nmi(conf[m], lab[:len(conf)][m])
        sw = (lab[:len(conf)] > 0) & (conf >= 200) & (conf <= 230)
        if sw.sum() > 100:
            v = conf[sw] - 200.; u, inv = np.unique(lab[:len(conf)][sw], return_inverse=True); mm = np.bincount(inv, v) / np.bincount(inv); out['eta2'] = 1 - ((v - mm[inv]) ** 2).sum() / ((v - v.mean()) ** 2).sum()
    return out
A, B = organisation(a.a), organisation(a.b); m = (A['lab'] > 0) & (B['lab'] > 0)
line = '%s%s (K %d) H %.4f | %s (K %d) H %.4f | ratio %s/%s %.2f | label NMI between them %.3f' % (tag, A['name'], A['K'], A['H'], B['name'], B['K'], B['H'], A['name'], B['name'], A['H'] / B['H'] if B['H'] != 0 else np.nan, nmi(A['lab'][m], B['lab'][m]))
if a.coords:
    y = np.loadtxt(a.coords)[:, 0]; mm = m & np.isfinite(y); u, inv = np.unique(A['lab'][mm], return_inverse=True); mu = np.bincount(inv, y[mm]) / np.bincount(inv)
    line += ' | share of the PC1 coordinate variance the %s explain %.3f' % (A['name'], 1 - ((y[mm] - mu[inv]) ** 2).sum() / ((y[mm] - y[mm].mean()) ** 2).sum())
print(line + '   [band %g-%g A]' % (lo, hi))
if 'Htrue' in A:
    print('%s  TRUTH: H on truth maps %s %.4f | %s %.4f | ratio %.2f ; NMI(conformer) %s %.3f | %s %.3f ; swing eta2 %s %s | %s %s' % (tag, A['name'], A['Htrue'], B['name'], B['Htrue'], A['Htrue'] / B['Htrue'] if B['Htrue'] else np.nan,
          A['name'], A['nmi'], B['name'], B['nmi'], A['name'], '%.3f' % A['eta2'] if 'eta2' in A else '-', B['name'], '%.3f' % B['eta2'] if 'eta2' in B else '-'))
