#!/usr/bin/env python3
"""Is the mobile region truncated by masking? For the closed-arm region (arm_readout.py's definition: where the deposited closed map has
density the consensus lacks) and the open-arm region, per radial bin (distance from the box centre, BIN A): the regions' voxel counts,
each MAP's density there relative to the deposited closed map's (lp 10, standardised as arm_readout), and each MASK's coverage of the
region (mean mask value); plus the spherical mask's soft edge for mskdiam MSKDIAM (SIMPLE-like cosine edge is not modelled: a hard
radius is printed). usage: radial_extent.py CONS CLOSED --maps name:file ... --masks name:file ... [--bin 5 --mskdiam 232]"""
import argparse, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('cons'); ap.add_argument('closed'); ap.add_argument('--maps', nargs='*', default=[]); ap.add_argument('--masks', nargs='*', default=[])
ap.add_argument('--bin', type=float, default=5.); ap.add_argument('--mskdiam', type=float, default=232.); a = ap.parse_args()
def rd(p):
    h = np.fromfile(p, np.int32, 256); n = int(h[0]); return np.fromfile(p, np.float32, offset=1024 + int(h[23]), count=n ** 3).reshape(n, n, n).astype(np.float64)
smpd = 1.464; LP = 10.; cons0 = rd(a.cons); n = cons0.shape[0]; f = np.fft.fftfreq(n, d=smpd)
kx, ky, kz = np.meshgrid(f, f, np.fft.rfftfreq(n, d=smpd), indexing='ij'); G = np.exp(-0.5 * (kx ** 2 + ky ** 2 + kz ** 2) * (LP / 1.2) ** 2)
lp = lambda v: np.fft.irfftn(np.fft.rfftn(v) * G, s=v.shape)
cons = lp(cons0); supp = cons > 0.2 * cons.max(); outside = cons < 0.02 * cons.max()
def std(v): v = v - v[outside].mean(); return v / max(np.sqrt(np.mean(v[supp] ** 2)), 1e-12)
cons = std(cons); closed = std(lp(rd(a.closed))); d = closed - cons; thr = 0.5 * np.percentile(np.abs(d[supp | (closed > 0.2 * closed.max())]), 99)
mC, mO = d > thr, d < -thr; g = np.arange(n) - n // 2; R = np.sqrt(g[:, None, None] ** 2 + g[None, :, None] ** 2 + g[None, None, :] ** 2) * smpd
edges = np.arange(0, R.max() + a.bin, a.bin); rb = np.digitize(R, edges) - 1
print('box %d x %.3f A (half-box %.0f A); spherical mask radius %.0f A; closed-arm region %d voxels, open-arm %d' % (n, smpd, n * smpd / 2, a.mskdiam / 2, mC.sum(), mO.sum()))
for nmR, reg in (('closed-arm region', mC), ('open-arm region', mO)):
    rs = [i for i in range(len(edges) - 1) if (reg & (rb == i)).sum() >= 5]
    if not rs: continue
    print('== %s: radius bins (A) %s' % (nmR, ' '.join('%3.0f-%-3.0f' % (edges[i], edges[i + 1]) for i in rs)))
    print('   %-26s %s' % ('voxels', ' '.join('%7d' % (reg & (rb == i)).sum() for i in rs)))
    ref = closed if nmR.startswith('closed') else cons
    for spec in a.maps:
        nm, fn = spec.split(':', 1); v = std(lp(rd(fn)))
        print('   %-26s %s' % (nm[:26], ' '.join('%7.2f' % (v[reg & (rb == i)].mean() / max(ref[reg & (rb == i)].mean(), 1e-6)) for i in rs)))
    for spec in a.masks:
        nm, fn = spec.split(':', 1); m = rd(fn)
        print('   %-26s %s' % (('mask ' + nm)[:26], ' '.join('%7.2f' % m[reg & (rb == i)].mean() for i in rs)))
    print('   %-26s %s' % ('sphere (hard) inside', ' '.join('%7.2f' % (R[reg & (rb == i)] <= a.mskdiam / 2).mean() for i in rs)))
