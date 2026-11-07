#!/usr/bin/env python3
"""land_truthocc.py: closed-site occupancy of every truth map with land_occ.py's crop-48 readout (same masks, lp, standardisation, docking
omitted: truth maps share the consensus frame). Writes one line of NG values. usage: land_truthocc.py LEAVES CONS CLOSED OUT [--crop 48 --mskdiam 230]"""
import sys, argparse, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('leaves'); ap.add_argument('cons'); ap.add_argument('closed'); ap.add_argument('out'); ap.add_argument('--crop', type=int, default=48)
ap.add_argument('--mskdiam', type=float, default=230.); ap.add_argument('--smpd', type=float, default=1.464); a = ap.parse_args()
def rd(p):
    h = np.fromfile(p, np.int32, 256); n_ = int(h[0]); return np.fromfile(p, np.float32, offset=1024 + int(h[23]), count=n_ ** 3).reshape(n_, n_, n_).astype(np.float64)
n = rd(a.cons).shape[0]; c = a.crop; R0 = a.mskdiam / 2 / a.smpd
g3 = np.arange(c) - c // 2; r3 = np.sqrt(g3[:, None, None] ** 2 + g3[None, :, None] ** 2 + g3[None, None, :] ** 2) * n / c; msk3 = np.clip((R0 + 4 - r3) / 8., 0, 1).astype(np.float32)
def crop_four(v):
    nb = v.shape[0]; sb = np.r_[0:c // 2, nb - c // 2:nb]; Fv = np.fft.fftn(np.fft.ifftshift(v)); vc = np.fft.fftshift(np.real(np.fft.ifftn(Fv[np.ix_(sb, sb, sb)]))) * msk3
    return np.fft.fftn(np.fft.ifftshift(vc)).ravel()
gq = np.fft.fftfreq(c) * c / n / a.smpd; qx, qy, qz = np.meshgrid(gq, gq, gq, indexing='ij'); GL = np.exp(-0.5 * (qx ** 2 + qy ** 2 + qz ** 2) * (10. / 1.2) ** 2).ravel()
def real_lp(Fm): return np.fft.fftshift(np.real(np.fft.ifftn((Fm * GL).reshape(c, c, c))))
cons = real_lp(crop_four(rd(a.cons))); supp = cons > 0.2 * cons.max(); outside = cons < 0.02 * cons.max()
def std(v): v = v - v[outside].mean(); return v / max(np.sqrt(np.mean(v[supp] ** 2)), 1e-12)
consS = std(cons); clS = std(real_lp(crop_four(rd(a.closed)))); dd = clS - consS; thr = 0.5 * np.percentile(np.abs(dd[supp | (clS > 0.2 * clS.max())]), 99)
mC = dd > thr; refC = clS[mC].mean()
RT = [std(real_lp(crop_four(rd(l.split()[1]))))[mC].mean() / refC for l in open(a.leaves) if l.strip()]
open(a.out, 'w').write(' '.join('%.4f' % x for x in RT) + '\n'); print('truth closed-site occupancy: ' + ' '.join('%.2f' % x for x in RT))
