#!/usr/bin/env python3
"""land_rank.py: TWO-STATE vs SEQUENTIAL test along an ordered series of maps (windows of a coordinate). The arm's closed-site region and
open-site region (arm_readout.py: deposited closed map vs consensus, lp 10, standardised) are cut into radial shells around the box centre;
each sub-region's occupancy is read per map. If particles only exchange between two states, every sub-region is the same linear function
of the closed fraction: the centred (sub-regions x maps) matrix is rank 1. A second component above the half-map noise (singular value
of the even + odd average vs the even - odd difference, both centred) means the sub-regions change at different points of the path:
a sequence (intermediates). Also printed: per sub-region the map index where its occupancy crosses half way between its end values.
usage: land_rank.py CONS CLOSED OUT.txt --maps M1 .. Mn [--halves E1,O1 .. En,On] [--bins 40,60,75,90,100,115] | --mix (principle test on
       the truth swing maps: --truth DIR, exact two-state mixtures of t=0 and t=30 vs single maps along t)"""
import argparse, glob, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('cons'); ap.add_argument('closed'); ap.add_argument('out'); ap.add_argument('--maps', nargs='*', default=[]); ap.add_argument('--halves', nargs='*', default=[])
ap.add_argument('--bins', default='40,60,75,90,100,118'); ap.add_argument('--truth', default=''); ap.add_argument('--smpd', type=float, default=1.464); a = ap.parse_args()
def rd(p):
    h = np.fromfile(p, np.int32, 256); n = int(h[0]); return np.fromfile(p, np.float32, offset=1024 + int(h[23]), count=n ** 3).reshape(n, n, n).astype(np.float64)
c0 = rd(a.cons); n = c0.shape[0]; f = np.fft.fftfreq(n, d=a.smpd); kx, ky, kz = np.meshgrid(f, f, f, indexing='ij'); G = np.exp(-0.5 * (kx ** 2 + ky ** 2 + kz ** 2) * (10. / 1.2) ** 2)
lp = lambda v: np.real(np.fft.ifftn(np.fft.fftn(v) * G))
cons = lp(c0); supp = cons > 0.2 * cons.max(); outside = cons < 0.02 * cons.max()
def std(v): v = v - v[outside].mean(); return v / max(np.sqrt(np.mean(v[supp] ** 2)), 1e-12)
consS = std(cons); clS = std(lp(rd(a.closed))); d = clS - consS; thr = 0.5 * np.percentile(np.abs(d[supp | (clS > 0.2 * clS.max())]), 99); mC = d > thr; mO = d < -thr
g = (np.arange(n) - n / 2.) * a.smpd; R = np.sqrt(g[:, None, None] ** 2 + g[None, :, None] ** 2 + g[None, None, :] ** 2); edges = [float(x) for x in a.bins.split(',')]
regs = []
for nm, m in (('closed', mC), ('open', mO)):
    for lo, hi in zip(edges[:-1], edges[1:]):
        r = m & (R >= lo) & (R < hi)
        if r.sum() >= 30: regs.append(('%s %3.0f-%3.0f A' % (nm, lo, hi), r, (clS if nm == 'closed' else consS)[r].mean()))
def readout(v): v = std(lp(v)); return np.array([v[r].mean() / ref for _, r, ref in regs])
def rank_test(Mavg, Mdif=None):
    X = Mavg - Mavg.mean(1, keepdims=True); s = np.linalg.svd(X, compute_uv=False)
    if Mdif is None: return s, None
    E = Mdif - Mdif.mean(1, keepdims=True); se = np.linalg.svd(E, compute_uv=False); return s, se
def crossings(Mavg):
    out = []
    for i in range(Mavg.shape[0]):
        y = Mavg[i]; y0, y1 = y[0], y[-1]; h = 0.5 * (y0 + y1); s = np.sign(y - h); k = np.where(s[:-1] != s[1:])[0]
        out.append(float(k[0] + (h - y[k[0]]) / (y[k[0] + 1] - y[k[0]])) if len(k) else np.nan)
    return np.array(out)
lines = ['# sub-regions: ' + '; '.join('%s (%d vox)' % (nm, r.sum()) for nm, r, _ in regs)]
if a.truth:
    tr = sorted(glob.glob(a.truth + '/conf2??_*.mrc')); V = [rd(p) for p in tr]; T0, T1 = V[0], V[-1]
    A2 = np.array([readout((1 - x) * T0 + x * T1) for x in np.linspace(0, 1, 11)]).T; As = np.array([readout(v) for v in V[::3]]).T
    for lab, M in (('two-state mixtures t0/t30', A2), ('single maps along t (rigid swing)', As)):
        s, _ = rank_test(M); lines.append('%-36s singular values %s | second/first %.3f | half-way crossings %s' % (lab, ' '.join('%.3f' % x for x in s[:3]), s[1] / s[0],
            ' '.join('%.1f' % x for x in crossings(M))))
else:
    M = np.array([readout(rd(p)) for p in a.maps]).T
    if a.halves:
        HE = np.array([readout(rd(h.split(',')[0])) for h in a.halves]).T; HO = np.array([readout(rd(h.split(',')[1])) for h in a.halves]).T
        cl = np.array([nm.startswith('closed') for nm, _, _ in regs])
        for lab, sel in (('all sub-regions', np.ones(len(regs), bool)), ('closed-site sub-regions', cl)):
            s, se = rank_test(0.5 * (HE + HO)[sel], 0.5 * (HE - HO)[sel])
            lines.append('%s: singular values (even+odd)/2 %s | noise (even-odd)/2 %s | second / noise-first %.2f (two-state if <~1)' % (lab, ' '.join('%.3f' % x for x in s[:4]),
                         ' '.join('%.3f' % x for x in se[:4]), s[1] / max(se[0], 1e-12)))
    else:
        s, _ = rank_test(M); lines.append('singular values %s (no halves: no noise reference)' % ' '.join('%.3f' % x for x in s[:4]))
    lines.append('half-way crossing (map index, 0 = first map): ' + ' | '.join('%s %.2f' % (nm, x) for (nm, _, _), x in zip(regs, crossings(M))))
    for (nm, _, _), row in zip(regs, M): lines.append('  %-18s ' % nm + ' '.join('%.3f' % x for x in row))
open(a.out, 'w').write('\n'.join(lines) + '\n'); print('\n'.join(lines))
