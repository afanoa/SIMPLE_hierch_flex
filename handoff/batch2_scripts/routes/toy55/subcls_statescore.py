#!/usr/bin/env python3
"""Per-state map truth score of a K-state reconstruction (recvol_stateNN_{even,odd}_unfil.mrc, any box: truth maps are Fourier-cropped
to it): per state its members, half-map FSC 0.143 resolution and the FIDELITY of its distinct features = per band corr between
(state map - members-weighted mean of the state maps) and the same for the members' conformer-weighted truth maps (junk left out);
summary = member-weighted mean fidelity per band 117-30 / 30-20 / 20-12 / 12-8 A, and the split-half reproducibility of the deviations.
usage: subcls_statescore.py RECDIR LABELS GT TRUTHMAPDIR [--mskdiam 220] [--tag NAME]"""
import sys, glob, os, argparse, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('rec'); ap.add_argument('labels'); ap.add_argument('gt'); ap.add_argument('tmaps')
ap.add_argument('--mskdiam', type=float, default=220.); ap.add_argument('--tag', default=''); a = ap.parse_args()
def rd(p):
    h = np.fromfile(p, np.int32, 256); n = int(h[0]); return np.fromfile(p, np.float32, offset=1024 + int(h[23]), count=n ** 3).reshape(n, n, n).astype(np.float64), float(np.fromfile(p, np.float32, 13)[10]) / n
lab = np.loadtxt(a.labels).astype(int); G = np.loadtxt(a.gt, comments='#', ndmin=2); conf = np.zeros(len(lab), int); conf[G[:, 0].astype(int) - 1] = G[:, 1].astype(int)
ks = [k for k in range(1, lab.max() + 1) if os.path.exists(os.path.join(a.rec, 'recvol_state%02d_even_unfil.mrc' % k)) and (lab == k).any()]
v0, smpd = rd(os.path.join(a.rec, 'recvol_state%02d_even_unfil.mrc' % ks[0])); n = v0.shape[0]; c = n // 2
zz, yy, xx = np.mgrid[:n, :n, :n]; rr = np.sqrt((xx - c) ** 2. + (yy - c) ** 2. + (zz - c) ** 2.); R = a.mskdiam / 2. / smpd; msk = np.clip((R + 3 - rr) / 6., 0, 1)
k = np.fft.fftfreq(n) * n; KZ, KY, KX = np.meshgrid(k, k, k, indexing='ij'); sh = np.sqrt(KX ** 2 + KY ** 2 + KZ ** 2).round().astype(int).ravel(); nsh = n // 2
def ft(v): return np.fft.fftn(v * msk).ravel()
def crop(v):  # centred Fourier crop of a truth map to the reconstruction box (same physical extent)
    m = v.shape[0]
    if m == n: return v
    F = np.fft.fftshift(np.fft.fftn(np.fft.ifftshift(v))); o = m // 2 - n // 2; return np.real(np.fft.fftshift(np.fft.ifftn(np.fft.ifftshift(F[o:o + n, o:o + n, o:o + n]))))
bands = [(117., 30.), (30., 20.), (20., 12.), (12., 8.)]; bm = [(sh >= n * smpd / lo) & (sh <= n * smpd / hi) for lo, hi in bands]
def bc(A, B): return np.array([(A[m] * B[m].conj()).real.sum() / np.sqrt((np.abs(A[m]) ** 2).sum() * (np.abs(B[m]) ** 2).sum() + 1e-30) for m in bm])
def fsc143(A, B):
    m = sh < nsh; num = np.bincount(sh[m], (A[m] * B[m].conj()).real, nsh); d1 = np.bincount(sh[m], np.abs(A[m]) ** 2, nsh); d2 = np.bincount(sh[m], np.abs(B[m]) ** 2, nsh)
    f = num / np.sqrt(np.maximum(d1 * d2, 1e-30)); below = np.where(f[1:] < 0.143)[0]; return n * smpd / max((below[0] + 1) if len(below) else nsh - 1, 1)
tm = {}
def tmap(cid):
    if cid not in tm:
        fl = glob.glob(os.path.join(a.tmaps, ('conf%02d_*_m160.mrc' % cid) if cid < 100 else ('conf%d_*_m160.mrc' % cid))); tm[cid] = ft(crop(rd(fl[0])[0])) if fl else None
    return tm[cid]
E, O, X, W, res = {}, {}, {}, {}, {}
for s in ks:
    m = lab == s; E[s] = ft(rd(os.path.join(a.rec, 'recvol_state%02d_even_unfil.mrc' % s))[0]); O[s] = ft(rd(os.path.join(a.rec, 'recvol_state%02d_odd_unfil.mrc' % s))[0])
    cs, cn = np.unique(conf[m & (conf > 0)], return_counts=True); acc, wsum = 0, 0
    for ci, ni in zip(cs, cn):
        t = tmap(int(ci))
        if t is not None: acc = acc + ni * t; wsum += ni
    X[s] = acc / max(wsum, 1); W[s] = m.sum(); res[s] = fsc143(E[s], O[s])
w = np.array([W[s] for s in ks], float); w /= w.sum()
Em = sum(wi * E[s] for wi, s in zip(w, ks)); Om = sum(wi * O[s] for wi, s in zip(w, ks)); Xm = sum(wi * X[s] for wi, s in zip(w, ks))
fid = np.array([bc(E[s] + O[s] - Em - Om, X[s] - Xm) for s in ks]); rep = np.array([bc(E[s] - Em, O[s] - Om) for s in ks])
print('%s%d states at box %d (%.2f A/px): state members / FSC0.143 A / deviation fidelity vs truth per band (117-30 30-20 20-12 12-8 A)' % ((a.tag + ': ') if a.tag else '', len(ks), n, smpd))
for i, s in enumerate(ks): print('  %2d %7d %5.1f   %s' % (s, W[s], res[s], ' '.join('%5.2f' % x for x in fid[i])))
print('%smember-weighted mean: deviation fidelity %s | deviation split-half reproducibility %s | median FSC0.143 %.1f A'
      % ((a.tag + ' ') if a.tag else '', ' '.join('%.2f' % x for x in w @ fid), ' '.join('%.2f' % x for x in w @ rep), np.median([res[s] for s in ks])))
