#!/usr/bin/env python3
"""Soft mask A without B: out = A x (1 - B) (both soft masks on the same box), re-softened (Gaussian SOFT px); optional readout of the
closed-arm / open-arm regions it covers (arm_readout.py's definition from CONS and CLOSED, used for verification only).
usage: mask_minus.py A.mrc B.mrc OUT.mrc [--soft 1.5 --cons consensus.mrc --closed closed.mrc]"""
import argparse, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('a'); ap.add_argument('b'); ap.add_argument('out'); ap.add_argument('--soft', type=float, default=1.5)
ap.add_argument('--cons', default=''); ap.add_argument('--closed', default=''); a = ap.parse_args()
def rd(p):
    h = np.fromfile(p, np.int32, 256); n = int(h[0]); return np.fromfile(p, np.float32, offset=1024 + int(h[23]), count=n ** 3).reshape(n, n, n).astype(np.float64)
def wr(p, x, smpd):
    h = np.zeros(256, np.int32); h[0:3] = x.shape; h[3] = 2; h[7:10] = x.shape; hf = h.view(np.float32); hf[10:13] = np.array(x.shape) * smpd; hf[13:16] = 90.
    h[16:19] = [1, 2, 3]; hf[19], hf[20], hf[21] = x.min(), x.max(), x.mean(); h[52] = int.from_bytes(b'MAP ', 'little'); h[53] = 0x4444
    with open(p, 'wb') as f: f.write(h.tobytes()); f.write(x.astype(np.float32).tobytes())
A, B = rd(a.a), rd(a.b); n = A.shape[0]; smpd = float(np.fromfile(a.a, np.float32, 13)[10]) / n
gw = np.minimum(np.arange(n), n - np.arange(n)); Rw = np.sqrt(gw[:, None, None] ** 2 + gw[None, :, None] ** 2 + gw[None, None, :] ** 2)
conv = lambda x, kk: np.fft.irfftn(np.fft.rfftn(x) * np.fft.rfftn(kk / kk.sum()), s=x.shape)
m = np.clip(conv(np.clip(A * (1. - np.clip(B, 0, 1)), 0, 1), np.exp(-Rw ** 2 / (2 * a.soft ** 2))), 0, 1); wr(a.out, m, smpd)
line = '%s = %s without %s: %.0f voxels above 0.5 (A %.0f, B %.0f)' % (a.out, a.a, a.b, (m > 0.5).sum(), (A > 0.5).sum(), (B > 0.5).sum())
if a.cons and a.closed:
    f = np.fft.fftfreq(n, d=smpd); kx, ky, kz = np.meshgrid(f, f, np.fft.rfftfreq(n, d=smpd), indexing='ij'); k = np.sqrt(kx ** 2 + ky ** 2 + kz ** 2)
    G10 = np.exp(-0.5 * (k * (10. / 1.2)) ** 2); l10 = lambda x: np.fft.irfftn(np.fft.rfftn(x) * G10, s=x.shape)
    cs = l10(rd(a.cons)); cl = l10(rd(a.closed)); s2 = cs > 0.2 * cs.max(); out_ = cs < 0.02 * cs.max()
    std = lambda x: (x - x[out_].mean()) / max(np.sqrt(np.mean((x - x[out_].mean())[s2] ** 2)), 1e-12); d = std(cl) - std(cs)
    thr = 0.5 * np.percentile(np.abs(d[s2 | (cl > 0.2 * cl.max())]), 99); mC, mO = d > thr, d < -thr
    for nm, x in (('A', A), ('B', B), ('out', m)):
        mb = x > 0.5; line += '\n   %s: closed-arm region covered %.2f, open-arm region covered %.2f' % (nm, (mb & mC).sum() / max(mC.sum(), 1), (mb & mO).sum() / max(mO.sum(), 1))
print(line)
