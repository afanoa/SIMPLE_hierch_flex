#!/usr/bin/env python3
"""Truth-free node-type statistics from a node's cheap view-balanced children (cheapkids output: children.txt + rec/*_reconstruct3D
recvol_stateNN_{even,odd}_unfil.mrc; children ordered along PC1 of the node's flex latent). The children are ASSIGNED with the flex
latent (band <= LPA, 15 A); every statistic below comes from cross-half products of the children's half-maps in bands the assignment
never saw (finer than LPA), after per-shell regression of every child on the children's weighted mean (removes scale, blur and
alignment-quality differences). The even and odd half-maps are independent, so neither the particles' noise nor the selection can
create them. Per band:
  R   cross-half reproducibility of the children's deviations, sum_k <dE_k,dO_k> / sqrt(sum_k|dE_k|^2 sum_k|dO_k|^2): heterogeneity
  L   linearity of the children path: top-eigenvalue share of the cross-half Gram matrix G_kl = Re<dE_k,dO_l> (symmetrised)
  B2  bend: share of the path's positive variance beyond its top component
  U   neighbour-difference profile p_k = Re<dE_k+1 - dE_k, dO_k+1 - dO_k>: ends / middle ratio (U-shape = unimodal density along a
      path, middle peak = a gap between clusters)
With --gt, the same statistics on the members' conformer-weighted TRUTH maps per child (noise-free geometry of the same children).
usage: nodetype_stats.py KIDSDIR [--lpa 15] [--mskdiam 220] [--gt gt --tmaps DIR] [--tag NAME]"""
import argparse, glob, os, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('kids'); ap.add_argument('--lpa', type=float, default=15.); ap.add_argument('--mskdiam', type=float, default=220.)
ap.add_argument('--gt', default=''); ap.add_argument('--tmaps', default=''); ap.add_argument('--tag', default=''); a = ap.parse_args()
def rd(p):
    h = np.fromfile(p, np.int32, 256); n = int(h[0]); return np.fromfile(p, np.float32, offset=1024 + int(h[23]), count=n ** 3).reshape(n, n, n).astype(np.float64), float(np.fromfile(p, np.float32, 13)[10]) / n
lab = np.loadtxt(os.path.join(a.kids, 'children.txt')).astype(int); RR = sorted(glob.glob(os.path.join(a.kids, 'rec', '*_reconstruct3D')))[-1]; K = lab.max()
w = np.array([(lab == k).sum() for k in range(1, K + 1)], float); w /= w.sum()
v0, smpd = rd(os.path.join(RR, 'recvol_state01_even_unfil.mrc')); n = v0.shape[0]; c = n // 2
g = np.arange(n) - c; R3 = np.sqrt(g[:, None, None] ** 2 + g[None, :, None] ** 2 + g[None, None, :] ** 2); r0 = a.mskdiam / 2 / smpd
M = np.where(R3 <= r0, 1., np.where(R3 <= r0 + 6, 0.5 + 0.5 * np.cos(np.pi * (R3 - r0) / 6), 0.))
f = np.fft.fftfreq(n); kx, ky, kz = np.meshgrid(f, f, np.fft.rfftfreq(n), indexing='ij'); kr = np.sqrt(kx ** 2 + ky ** 2 + kz ** 2)
res = np.where(kr > 0, smpd / np.maximum(kr, 1e-9), 1e9); keep = res >= 7.9; sh = np.round(kr * n).astype(int)[keep]; nsh = sh.max() + 1
wb = np.where((kz == 0) | (np.abs(kz) == 0.5), 1., 2.)[keep]; resk = res[keep]
bands = {'assign <=%g A (biased)' % a.lpa: (117., a.lpa), 'B %g-8 A' % a.lpa: (a.lpa, 8.), 'B1 %g-11 A' % a.lpa: (a.lpa, 11.), 'B2 11-8 A': (11., 8.)}
def F(x): return np.fft.rfftn(M * x)[keep]
def regress(E, O):  # per-shell regression of every child on the weighted mean, cross-half coefficients (focus_mask.py)
    Em, Om = w @ E, w @ O; ss = lambda z: np.bincount(sh, weights=wb * z, minlength=nsh)
    den = ss(np.real(Em * Om.conj())); ok = np.abs(den) > 1e-9 * np.abs(den).max(); den = np.where(ok, den, 1.)
    dE = np.array([e - np.where(ok, ss(np.real(e * Om.conj())) / den, 1.)[sh] * Em for e in E]); dO = np.array([o - np.where(ok, ss(np.real(o * Em.conj())) / den, 1.)[sh] * Om for o in O])
    return dE, dO
def stats(dE, dO, lo, hi):
    m = (resk <= lo) & (resk >= hi); ww = wb[m]; A = dE[:, m]; B = dO[:, m]
    G = np.real((A * ww) @ B.conj().T); G = 0.5 * (G + G.T)
    Rv = np.trace(G) / np.sqrt(np.real((np.abs(A) ** 2) @ ww).sum() * np.real((np.abs(B) ** 2) @ ww).sum() + 1e-300)
    lam = np.sort(np.linalg.eigvalsh(G))[::-1]; pos = np.maximum(lam, 0); L = pos[0] / max(pos.sum(), 1e-300)
    dA = A[1:] - A[:-1]; dB = B[1:] - B[:-1]; p = np.real(np.sum(dA * dB.conj() * ww, 1)); pm = p / max(np.abs(p).mean(), 1e-300)
    kk = len(p); ends = np.r_[pm[:2], pm[-2:]].mean(); mid = pm[kk // 2 - 1: kk // 2 + 2].mean(); U = ends / mid if mid > 0 else np.inf
    return Rv, L, 1 - L, U, pm
E = np.array([F(rd(os.path.join(RR, 'recvol_state%02d_even_unfil.mrc' % k))[0]) for k in range(1, K + 1)])
O = np.array([F(rd(os.path.join(RR, 'recvol_state%02d_odd_unfil.mrc' % k))[0]) for k in range(1, K + 1)])
dE, dO = regress(E, O)
tag = (a.tag + ' ') if a.tag else ''
print('%s%d children (members %d..%d), box %d @ %.3f A; assignment band <= %g A, measurement bands disjoint' % (tag, K, (lab == 1).sum(), (lab == K).sum(), n, smpd, a.lpa))
for bn, (lo, hi) in bands.items():
    Rv, L, Bd, U, pm = stats(dE, dO, lo, hi)
    print('%s  %-24s R %.3f  L %.3f  bend %.3f  U(ends/mid) %.2f  profile %s' % (tag, bn, Rv, L, Bd, U, ' '.join('%.2f' % x for x in pm)))
if a.gt and a.tmaps:
    G_ = np.loadtxt(a.gt, comments='#', ndmin=2); conf = np.zeros(len(lab), int); conf[G_[:, 0].astype(int) - 1] = G_[:, 1].astype(int)
    cache = {}; X = []
    for k in range(1, K + 1):
        cs, cn = np.unique(conf[(lab == k) & (conf > 0)], return_counts=True); acc, ws = 0, 0
        for ci, ni in zip(cs, cn):
            if ci not in cache:
                fl = glob.glob(os.path.join(a.tmaps, ('conf%02d_*_m160.mrc' % ci) if ci < 100 else ('conf%d_*_m160.mrc' % ci))); cache[ci] = F(rd(fl[0])[0]) if fl else None
            if cache[ci] is not None: acc = acc + ni * cache[ci]; ws += ni
        X.append(acc / max(ws, 1))
    X = np.array(X); dX, _ = regress(X, X)
    for bn, (lo, hi) in bands.items():
        Rv, L, Bd, U, pm = stats(dX, dX, lo, hi)
        print('%s  TRUTH-maps %-13s          L %.3f  bend %.3f  U(ends/mid) %.2f  profile %s' % (tag, bn, L, Bd, U, ' '.join('%.2f' % x for x in pm)))
    sw = (conf >= 200) & (conf <= 230); nj = conf > 0
    print('%s  TRUTH: junk %.1f %%, swing %.1f %% of members; per child swing share / mean swing: %s' % (tag, 100 * (conf[lab > 0] == 0).mean(), 100 * sw[lab > 0].mean(),
          ' '.join('%.0f%%/%s' % (100 * sw[lab == k].mean(), ('%.1f' % (conf[(lab == k) & sw] - 200).mean()) if (sw & (lab == k)).any() else '-') for k in range(1, K + 1))))
    top = np.unique(conf[nj & (lab > 0)], return_counts=True); o = np.argsort(-top[1])[:4]
    print('%s  TRUTH: top conformers %s; their share per child: %s' % (tag, ' '.join('%d' % top[0][i] for i in o),
          ' | '.join(' '.join('%.2f' % (conf[lab == k] == top[0][i]).mean() for i in o) for k in range(1, K + 1))))
