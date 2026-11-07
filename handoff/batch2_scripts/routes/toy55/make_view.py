#!/usr/bin/env python3
"""Viewing copies of state maps: the average of the unfiltered even/odd halves, low-passed to LP A (raised-cosine edge LP+1 -> LP-0.5 A),
unsharpened, written as OUTDIR/<name>.mrc (full maps only). usage: make_view.py OUTDIR SPEC [--lp 10]; SPEC lines: name even.mrc odd.mrc"""
import os, sys, argparse, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('out'); ap.add_argument('spec'); ap.add_argument('--lp', type=float, default=10.); ap.add_argument('--box', type=int, default=0); a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
def rd(p):
    h = np.fromfile(p, np.int32, 256); n = int(h[0]); return np.fromfile(p, np.float32, offset=1024 + int(h[23]), count=n ** 3).reshape(n, n, n).astype(np.float64), float(np.fromfile(p, np.float32, 13)[10]) / n
def wr(p, x, smpd):
    h = np.zeros(256, np.int32); h[0:3] = x.shape; h[3] = 2; h[7:10] = x.shape; hf = h.view(np.float32); hf[10:13] = np.array(x.shape) * smpd; hf[13:16] = 90.
    h[16:19] = [1, 2, 3]; hf[19], hf[20], hf[21] = x.min(), x.max(), x.mean(); h[52] = int.from_bytes(b'MAP ', 'little'); h[53] = 0x4444
    with open(p, 'wb') as f: f.write(h.tobytes()); f.write(x.astype(np.float32).tobytes())
for l in open(a.spec):
    if not l.strip() or l.startswith('#'): continue
    nm, e, o = l.split()[:3]; E, smpd = rd(e); O, _ = rd(o); n = E.shape[0]
    f = np.fft.fftfreq(n); kx, ky, kz = np.meshgrid(f, f, np.fft.rfftfreq(n), indexing='ij'); k = np.sqrt(kx ** 2 + ky ** 2 + kz ** 2) / smpd
    k1, k2 = 1. / (a.lp + 1.), 1. / (a.lp - 0.5); filt = np.where(k <= k1, 1., np.where(k >= k2, 0., 0.5 + 0.5 * np.cos(np.pi * (k - k1) / (k2 - k1))))
    V = np.fft.irfftn(np.fft.rfftn((E + O) / 2) * filt, s=E.shape)
    if a.box and a.box > n:   # Fourier-pad to the target box (same physical size, finer sampling)
        F = np.fft.fftshift(np.fft.fftn(V)); P = np.zeros((a.box,) * 3, complex); lo = a.box // 2 - n // 2; P[lo:lo + n, lo:lo + n, lo:lo + n] = F
        V = np.real(np.fft.ifftn(np.fft.ifftshift(P))) * (a.box / n) ** 3; smpd = smpd * n / a.box
    wr('%s/%s.mrc' % (a.out, nm), V, smpd)
print('%s: %d maps' % (a.out, sum(1 for l in open(a.spec) if l.strip() and not l.startswith('#'))))
