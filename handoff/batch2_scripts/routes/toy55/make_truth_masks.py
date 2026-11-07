#!/usr/bin/env python3
"""Tight 'truth' site masks for scoring site discovery.
arms MAPDIR SECTORDIR OUT [--nmaps 60 --frac 0.1]: IgM-like set - the cross-conformer variance of the true maps (lp 15) above FRAC of its
   maximum, assigned to the arm sectors (SECTORDIR/arm_k.mrc), dilated 2 px, soft: OUT/arm_k_tight.mrc (the arm's swept density).
rnasep CONS CLOSED OUT: real RNase P - arm_readout.py's closed-arm region (density the deposited closed map has and the consensus lacks)
   and open-arm region (the reverse), lp 10 standardised, dilated 2 px, soft: OUT/closed_arm_region.mrc, OUT/open_arm_region.mrc."""
import sys, os, glob, numpy as np
from scipy import ndimage
def rd(p):
    h = np.fromfile(p, np.int32, 256); n = int(h[0]); return np.fromfile(p, np.float32, offset=1024 + int(h[23]), count=n ** 3).reshape(n, n, n).astype(np.float64)
def wr(p, x, smpd=1.464):
    x = np.ascontiguousarray(x, np.float32); h = np.zeros(256, np.int32); h[0:3] = x.shape; h[3] = 2; h[7:10] = x.shape; hf = h.view(np.float32); hf[10:13] = np.array(x.shape) * smpd; hf[13:16] = 90.
    h[16:19] = [1, 2, 3]; hf[19], hf[20], hf[21] = x.min(), x.max(), x.mean(); h[52] = int.from_bytes(b'MAP ', 'little'); h[53] = 0x4444
    with open(p, 'wb') as f: f.write(h.tobytes()); f.write(x.tobytes())
def lpg(v, res, smpd=1.464):
    n = v.shape[0]; f = np.fft.fftfreq(n, d=smpd); kx, ky, kz = np.meshgrid(f, f, np.fft.rfftfreq(n, d=smpd), indexing='ij')
    return np.fft.irfftn(np.fft.rfftn(v) * np.exp(-0.5 * (kx ** 2 + ky ** 2 + kz ** 2) * (res / 1.2) ** 2), s=v.shape)
soft = lambda b: np.clip(ndimage.gaussian_filter(ndimage.binary_dilation(b, iterations=2).astype(float), 1.5) * 1.5, 0, 1)
mode = sys.argv[1]
if mode == 'arms':
    md, sd, out = sys.argv[2:5]; nm = int(sys.argv[5]) if len(sys.argv) > 5 else 60; frac = float(sys.argv[6]) if len(sys.argv) > 6 else 0.1; os.makedirs(out, exist_ok=True)
    fs = sorted(glob.glob(md + '/conf*_synth.mrc'))[:nm]; S = None; S2 = None
    for f in fs:
        v = lpg(rd(f), 15.); S = v if S is None else S + v; S2 = v * v if S2 is None else S2 + v * v
    var = S2 / len(fs) - (S / len(fs)) ** 2; core = var > frac * var.max(); secs = sorted(glob.glob(sd + '/arm_*.mrc'))
    for p in secs:
        sec = rd(p) > 0.5; m = soft(core & sec); k = os.path.basename(p).replace('.mrc', ''); wr('%s/%s_tight.mrc' % (out, k), m)
        print('%s: variance core in the sector %d voxels -> mask > 0.5 %d voxels' % (k, (core & sec).sum(), (m > 0.5).sum()))
    print('variance core %d voxels from %d maps; %d of them in no sector' % (core.sum(), len(fs), (core & ~sum(rd(p) > 0.5 for p in secs).astype(bool)).sum()))
else:
    cons_f, closed_f, out = sys.argv[2:5]; os.makedirs(out, exist_ok=True)
    cs = lpg(rd(cons_f), 10.); supp = cs > 0.2 * cs.max(); outside = cs < 0.02 * cs.max()
    def std(v): v = v - v[outside].mean(); return v / max(np.sqrt(np.mean(v[supp] ** 2)), 1e-12)
    c_ = std(cs); cl = std(lpg(rd(closed_f), 10.)); d = cl - c_; thr = 0.5 * np.percentile(np.abs(d[supp | (cl > 0.2 * cl.max())]), 99)
    for nm, b in (('closed_arm_region', d > thr), ('open_arm_region', d < -thr)): m = soft(b); wr('%s/%s.mrc' % (out, nm), m); print('%s: %d voxels -> mask > 0.5 %d' % (nm, b.sum(), (m > 0.5).sum()))
